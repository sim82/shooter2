/*
 * Copyright (C) 2012 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */

#include <ClanLib/display.h>
#include <ClanLib/application.h>
#include <ClanLib/gl.h>
#include <ClanLib/core.h>
#include <GL/gl.h>
#include <GL/glx.h>
#include <CL/cl_gl.h>

#define __CL_ENABLE_EXCEPTIONS
#include "cl.hpp"
#include "cycle.h"
#include "cl_error_codes.h"
#include <execinfo.h>
#include <fstream>
#include <array>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <thread>
#include <atomic>
#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "aligned_buffer.h"
#include "vec_unit.h"

namespace ublas = boost::numeric::ublas;


class spinlock_mutex {
public:
    spinlock_mutex() : flag_(ATOMIC_FLAG_INIT) {}
    
    void lock() {
        while( flag_.test_and_set( std::memory_order_acquire ));
    }
    
    void unlock() {
        flag_.clear();
    }
    
private:
    std::atomic_flag flag_;
};

typedef CL_Vec3i vec3i;
typedef CL_Vec3f vec3f;

struct col3f_sse {
    float r;
    float g;
    float b;
    float x;

    col3f_sse( float r_, float g_, float b_ ) : r(r_), g(g_), b(b_) {}
    col3f_sse( const vec3f & v ) : r(v.r), g(v.g), b(v.b), x(0) {} // { std::cout << "assign: " << r << " " << g << " " << b << "\n";}
    const col3f_sse &operator=( const col3f_sse &other ) {
        r = other.r;
        g = other.g;
        b = other.b;
        x = 0;
        return *this;
    }

    col3f_sse() {}

//  inline operator vec3f() {
//      return vec3f(r, g, b );
//  }

};

template<typename vec_t>
typename vec_t::datatype dist_sqr( const vec_t &v1, const vec_t &v2 ) {
    vec_t vd = v1 - v2;

    return vd.dot(vd);
}



template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}



class bitmap3d : private boost::dynamic_bitset<> {

public:
    typedef boost::dynamic_bitset<>::reference reference;


    bitmap3d( size_t x, size_t y, size_t z )
            : boost::dynamic_bitset<>(x * y * z, false),
            x_(x), y_(y), z_(z),
            slice_(x * z),
            stride_(x)
    {


    }

    bitmap3d()
            :
            x_(0), y_(0), z_(0),
            slice_(0),
            stride_(0)
    {


    }


//  bitmap3d() : bitmap3d(0, 0, 0) // aaahhhhrgggg waiting for gcc 4.7....
//  {
//
//  }


    reference operator()( int x, int y, int z ) {
        if ( !inside(x, y, z)) {
            throw std::runtime_error( "index out of bounds");
        }

        return operator[]( addr(x, y, z) );
    }

    bool operator()( int x, int y, int z ) const {
        if ( !inside( x, y, z )) {
            return false;
        }

        return operator[]( addr(x, y, z) );
    }


    inline size_t x() const {
        return x_;
    }
    inline size_t y() const {
        return y_;
    }
    inline size_t z() const {
        return z_;
    }
    
private:
    inline bool inside( int x, int y, int z ) const {
        return !(x < 0 || y < 0 || z < 0 || x >= int(x_) || y >= int(y_) || z >= int(z_) );
    }

    inline size_t addr( size_t x, size_t y, size_t z ) const {
        return  x + stride_ * z + slice_ * y;
    }

    size_t x_;
    size_t y_;
    size_t z_;

    size_t slice_;
    size_t stride_;

};


class util {
public:


    static bool occluded2(vec3i p0, vec3i p1, const bitmap3d &solid)
    {
        // 3d bresenham, ripped from http://www.cobrabytes.com/index.php?topic=1150.0



        int x0 = p0.x;
        int y0 = p0.y;
        int z0 = p0.z;

        int x1 = p1.x;
        int y1 = p1.y;
        int z1 = p1.z;



        //'steep' xy Line, make longest delta x plane
        const bool swap_xy = abs(y1 - y0) > abs(x1 - x0);
        if ( swap_xy ) {
            std::swap(x0, y0);
            std::swap(x1, y1);
        }

        //do same for xz
        const bool swap_xz = abs(z1 - z0) > abs(x1 - x0);
        if ( swap_xz ) {
            std::swap(x0, z0);
            std::swap(x1, z1);
        }

        //delta is Length in each plane
        int delta_x = abs(x1 - x0);
        int delta_y = abs(y1 - y0);
        int delta_z = abs(z1 - z0);

        //drift controls when to step in 'shallow' planes
        //starting value keeps Line centred
        int drift_xy  = (delta_x / 2);
        int drift_xz  = (delta_x / 2);

        //direction of line
        const int step_x = (x0 > x1) ? -1 : 1;
        const int step_y = (y0 > y1) ? -1 : 1;
        const int step_z = (z0 > z1) ? -1 : 1;

        //starting point
        int y = y0;
        int z = z0;

        //step through longest delta (which we have swapped to x)


        if ( swap_xz ) {
            for ( int x = x0; x != x1; x += step_x ) {

                //copy position
                int cx = z;
                int cy = y;
                int cz = x;

                //unswap (in reverse)

                //std::swap(cx, cz);


                //passes through this point
                //debugmsg(":" + cx + ", " + cy + ", " + cz)
                if ( solid( cx, cy, cz ) ) {
                    return true;
                }
                //update progress in other planes
                drift_xy = drift_xy - delta_y;
                drift_xz = drift_xz - delta_z;

                //step in y plane
                if ( drift_xy < 0 ) {
                    y = y + step_y;
                    drift_xy = drift_xy + delta_x;
                }


                //same in z
                if ( drift_xz < 0 ) {
                    z = z + step_z;
                    drift_xz = drift_xz + delta_x;
                }
            }
        } else if ( swap_xz && swap_xy ) {
            for ( int x = x0; x != x1; x += step_x ) {

                //copy position
                int cx = y;
                int cy = z;
                int cz = x;

                //unswap (in reverse)

//              std::swap(cx, cz);
//              std::swap(cx, cy);

                //passes through this point
                //debugmsg(":" + cx + ", " + cy + ", " + cz)
                if ( solid( cx, cy, cz ) ) {
                    return true;
                }
                //update progress in other planes
                drift_xy = drift_xy - delta_y;
                drift_xz = drift_xz - delta_z;

                //step in y plane
                if ( drift_xy < 0 ) {
                    y = y + step_y;
                    drift_xy = drift_xy + delta_x;
                }


                //same in z
                if ( drift_xz < 0 ) {
                    z = z + step_z;
                    drift_xz = drift_xz + delta_x;
                }
            }

        } else if ( swap_xy ) {
            for ( int x = x0; x != x1; x += step_x ) {

                //copy position
                int cx = y;
                int cy = x;
                int cz = z;

                //unswap (in reverse)


//              std::swap(cx, cy);

                //passes through this point
                //debugmsg(":" + cx + ", " + cy + ", " + cz)
                if ( solid( cx, cy, cz ) ) {
                    return true;
                }
                //update progress in other planes
                drift_xy = drift_xy - delta_y;
                drift_xz = drift_xz - delta_z;

                //step in y plane
                if ( drift_xy < 0 ) {
                    y = y + step_y;
                    drift_xy = drift_xy + delta_x;
                }


                //same in z
                if ( drift_xz < 0 ) {
                    z = z + step_z;
                    drift_xz = drift_xz + delta_x;
                }
            }

        } else {
            for ( int x = x0; x != x1; x += step_x ) {

                //copy position
                int cx = x;
                int cy = y;
                int cz = z;

                //passes through this point
                //debugmsg(":" + cx + ", " + cy + ", " + cz)
                if ( solid( cx, cy, cz ) ) {
                    return true;
                }
                //update progress in other planes
                drift_xy = drift_xy - delta_y;
                drift_xz = drift_xz - delta_z;

                //step in y plane
                if ( drift_xy < 0 ) {
                    y = y + step_y;
                    drift_xy = drift_xy + delta_x;
                }


                //same in z
                if ( drift_xz < 0 ) {
                    z = z + step_z;
                    drift_xz = drift_xz + delta_x;
                }
            }

        }


        return false;

    }
    static bool occluded(vec3i p0, vec3i p1, const bitmap3d &solid)
    {
        // 3d bresenham, ripped from http://www.cobrabytes.com/index.php?topic=1150.0



        int x0 = p0.x;
        int y0 = p0.y;
        int z0 = p0.z;

        int x1 = p1.x;
        int y1 = p1.y;
        int z1 = p1.z;



        //'steep' xy Line, make longest delta x plane
        const bool swap_xy = abs(y1 - y0) > abs(x1 - x0);
        if ( swap_xy ) {
            std::swap(x0, y0);
            std::swap(x1, y1);
        }

        //do same for xz
        const bool swap_xz = abs(z1 - z0) > abs(x1 - x0);
        if ( swap_xz ) {
            std::swap(x0, z0);
            std::swap(x1, z1);
        }

        //delta is Length in each plane
        int delta_x = abs(x1 - x0);
        int delta_y = abs(y1 - y0);
        int delta_z = abs(z1 - z0);

        //drift controls when to step in 'shallow' planes
        //starting value keeps Line centred
        int drift_xy  = (delta_x / 2);
        int drift_xz  = (delta_x / 2);

        //direction of line
        const int step_x = (x0 > x1) ? -1 : 1;
        const int step_y = (y0 > y1) ? -1 : 1;
        const int step_z = (z0 > z1) ? -1 : 1;

        //starting point
        int y = y0;
        int z = z0;

        //step through longest delta (which we have swapped to x)


        for ( int x = x0; x != x1; x += step_x ) {

            //copy position
            int cx = x;
            int cy = y;
            int cz = z;

            //unswap (in reverse)
            if ( swap_xz ) {
                std::swap(cx, cz);
            }

            if ( swap_xy ) {
                std::swap(cx, cy);
            }

            //passes through this point
            //debugmsg(":" + cx + ", " + cy + ", " + cz)
            if ( solid( cx, cy, cz ) ) {
                return true;
            }
            //update progress in other planes
            drift_xy = drift_xy - delta_y;
            drift_xz = drift_xz - delta_z;

            //step in y plane
            if ( drift_xy < 0 ) {
                y = y + step_y;
                drift_xy = drift_xy + delta_x;
            }


            //same in z
            if ( drift_xz < 0 ) {
                z = z + step_z;
                drift_xz = drift_xz + delta_x;
            }
        }

        return false;

    }
};

struct plane {
public:
    enum dir_type {
        dir_xy_p,
        dir_xy_n,
        dir_yz_p,
        dir_yz_n,
        dir_zx_p,
        dir_zx_n,
    };


    bool normal_cull( const plane &other ) const {
//      return (dir_ == dir_xy_p && other.dir_ == dir_xy_n) ||
//              (dir_ == dir_xy_n && other.dir_ == dir_xy_p) ||
//              (dir_ == dir_yz_p && other.dir_ == dir_yz_n) ||
//              (dir_ == dir_yz_n && other.dir_ == dir_yz_p) ||
//              (dir_ == dir_zx_p && other.dir_ == dir_zx_n) ||
//              (dir_ == dir_zx_n && other.dir_ == dir_zx_p);

        return dir_ == other.dir_;
    }


    static vec3f normal( dir_type dt ) {
        switch ( dt ) {
        case dir_xy_p:
            return vec3f( 0.0, 0.0, 1.0 );
        case dir_xy_n:
            return vec3f( 0.0, 0.0, -1.0 );
        case dir_yz_p:
            return vec3f( 1.0, 0.0, 0.0 );
        case dir_yz_n:
            return vec3f( -1.0, 0.0, 0.0 );
        case dir_zx_p:
            return vec3f( 0.0, 1.0, 0.0 );
        case dir_zx_n:
        default:
            return vec3f( 0.0, -1.0, 0.0 );
        }
    }

    static vec3f primary( dir_type dt ) {
        switch ( dt ) {
        case dir_xy_p:
        case dir_xy_n:
            return vec3f( 1.0, 1.0, 0.0 );
        case dir_yz_p:
        case dir_yz_n:
            return vec3f( 0.0, 1.0, 1.0 );
        case dir_zx_p:
        case dir_zx_n:
            return vec3f( 1.0, 0.0, 1.0 );
        }
    }

    static vec3f primary0( dir_type dt ) {
        switch ( dt ) {
        case dir_xy_p:
        case dir_xy_n:
            return vec3f( 1.0, 0.0, 0.0 );
        case dir_yz_p:
        case dir_yz_n:
            return vec3f( 0.0, 1.0, 0.0 );
        case dir_zx_p:
        case dir_zx_n:
        default:
            return vec3f( 0.0, 0.0, 1.0 );
        }
    }

    static vec3f primary1( dir_type dt ) {
        switch ( dt ) {
        case dir_xy_p:
        case dir_xy_n:
            return vec3f( 0.0, 1.0, 0.0 );
        case dir_yz_p:
        case dir_yz_n:
            return vec3f( 0.0, 0.0, 1.0 );
        case dir_zx_p:
        case dir_zx_n:
        default:
            return vec3f( 1.0, 0.0, 0.0 );
        }
    }

    static std::array<float,4> vgen0( dir_type dt ) {
//         std::array<float,4> v;

        switch ( dt ) {
        case dir_xy_p:
        case dir_yz_p:
        case dir_zx_p:
            return std::array<float,4> {{-0.5, 0.5, 0.5, -0.5 }};
//          return v;
//          break;
        case dir_xy_n:
        case dir_yz_n:
        case dir_zx_n:
        default:
            return std::array<float,4> {{-0.5, -0.5, 0.5, 0.5 }};
//          break;
        }

//      return v;
    }

    static std::array<float,4> vgen1( dir_type dt ) {
//      std::array<float,4> v;

        switch ( dt ) {
        case dir_xy_p:
        case dir_yz_p:
        case dir_zx_p:
            return std::array<float,4> {{-0.5, -0.5, 0.5, 0.5 }};

        case dir_xy_n:
        case dir_yz_n:
        case dir_zx_n:
        default:
            return std::array<float,4> {{-0.5, 0.5, 0.5, -0.5 }};
        }
    }

//  static vec3f col_diff( dir_type dt ) {
//      switch( dt ) {
//      case dir_xy_p:
//          return vec3f(1.0, 0.5, 0.0 );
//
//      case dir_yz_p:
//          return vec3f(0.0, 1.0, 0.0 );
//
//      case dir_xy_n:
//      case dir_yz_n:
//      case dir_zx_p:
//      case dir_zx_n:
//      default:
//          return vec3f(0.8,0.8,0.8);
//      }
//  }

    static vec3f col_diff( dir_type dt ) {
        switch ( dt ) {
        case dir_xy_p:
            return vec3f(1.0, 0.5, 0.0 );

        case dir_xy_n:
            return vec3f(0.0, 1.0, 0.0 );

        case dir_yz_p:
        case dir_yz_n:
        case dir_zx_p:
        case dir_zx_n:
        default:
            return vec3f(0.8,0.8,0.8);
        }
    }


    plane( dir_type d, const vec3f &base_pos, const vec3i &pos, float scale, float energy )
            : dir_(d),
            pos_(pos),
            energy_(energy),
            col_diff_(col_diff(d))
    {



        vec3f p0 = primary0(d);
        vec3f p1 = primary1(d);

        std::array<float,4> vg0 = vgen0(d);
        std::array<float,4> vg1 = vgen1(d);

        vec3f norm2 = normal(d) * 0.5;

        vec3f trans_pos = (base_pos + pos);

        for ( size_t i = 0; i < 4; ++i ) {
            verts_[i] = (trans_pos + p0 * vg0[i] + p1 * vg1[i] + norm2) * scale;

//          std::cout << "vert: " << i << " " << verts_[i] << " " << vg0[i] << "\n";
        }

        if ( (pos_.y / 2) % 2 == 0 ) {
            col_diff_ = vec3f(1.0, 1.0, 1.0 );
        }


    }

    void render() const {

        for ( size_t i = 0; i < 4; ++i ) {
            //glColor3f( energy_, energy_, energy_ );
            glColor3f( energy_rgb_.r, energy_rgb_.g, energy_rgb_.b );
            glVertex3fv( verts_[i] );

        }
    }

    template<typename oiter>
    void render_vertex( oiter first ) const {
        for ( size_t i = 0; i < 4; ++i ) {
            *(first++) = verts_[i].x;
            *(first++) = verts_[i].y;
            *(first++) = verts_[i].z;
            
            
        }
        
    }
    template<typename oiter>
    oiter render_color( oiter first ) const {
        for ( size_t i = 0; i < 4; ++i ) {
            *(first++) = energy_rgb_.r;
            *(first++) = energy_rgb_.g;
            *(first++) = energy_rgb_.b;
        }
        
        return first;
    }
    
    const vec3i &pos() const {
        return pos_;
    }
    vec3f norm() const {
        return normal(dir_);
    }

    dir_type dir() const {
        return dir_;
    }

    void energy( float e ) {
        energy_ = e;
    }
    void energy_rgb( const vec3f &e ) {
        energy_rgb_ = e;
    }

    const vec3f &col_diff() const {
        return col_diff_;
    }

private:
    dir_type dir_;
    vec3i pos_;

    std::array<vec3f,4> verts_;

    float energy_;
    vec3f energy_rgb_;
    vec3f col_diff_;

};


class rad_core {
public:
    virtual void set_emit( const std::vector<vec3f> &emit ) = 0;
//     virtual bool update() = 0;
    virtual void copy( std::vector<vec3f> *out ) = 0;
};


class rad_core_sse : public rad_core {
public:
    rad_core_sse( const std::vector<plane> &planes, const std::vector<std::vector<float> > &ffs, const std::vector<std::vector<int> > &ff_target )
            : emit_( planes.size() ), rad_( planes.size() ), rad2_( planes.size() ),
            ffs_(ffs), ff_target_(ff_target),
            planes_(planes)
    {


    }

    virtual void set_emit( const std::vector<vec3f> &emit ) {
        assert( emit_.size() == emit.size() );

        emit_.assign( std::begin(emit), std::end(emit));
    }

    virtual bool update() {
        do_radiosity_sse();

        return true;
    }

    virtual void copy( std::vector<vec3f> *out ) {
        for ( size_t i = 0; i < rad_.size(); ++i ) {
            (*out)[i].r = rad_[i].r;
            (*out)[i].g = rad_[i].g;
            (*out)[i].b = rad_[i].b;
        }
    }


private:

    void do_radiosity_sse() {
        const int steps = 1;
        const float min_ff = 0;
//      std::fill(e_rad_sse_.begin(), e_rad_sse_.end(), vec3f(0.0, 0.0, 0.0));
//      std::fill(e_rad2_sse_.begin(), e_rad2_sse_.end(), vec3f(0.0, 0.0, 0.0));

        //std::copy( emit_sse_.begin(), emit_sse_.end(), e_rad_sse_.begin() );

        typedef vector_unit<float,4> vu;
        typedef vu::vec_t vec_t;
        //steps = 0;

        vec_t reflex_factor = vu::set1(1.0);
        for ( int i = 0; i < steps; ++i ) {
            //for( auto it = pairs_.begin(); it != pairs_.end(); ++it, ++ff_it ) {


            for ( size_t j = 0; j < ffs_.size(); ++j ) {

                const size_t s = ffs_[j].size();
                vec_t rad = vu::set1(0);

                const vec3f cd = planes_[j].col_diff();
                const vec_t col_diff = vu::set( 0, cd.b, cd.g, cd.r );



                for ( size_t k = 0; k < s; ++k ) {
                    size_t target = ff_target_[j][k];
                    //rad_rgb += (col_diff * e_rad_rgb_[target]) * ff2s_[j][k];

                    if ( false && ffs_[j][k] < min_ff ) {
                        continue;
                    }

                    const vec_t ff = vu::set1( ffs_[j][k] );

                    rad = vu::add( rad, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target))), ff ));

                }

                vu::store( vu::add( vu::load((float*)emit_(j)), vu::mul(rad, reflex_factor)), (float*)rad2_(j));
//              std::cout << "col: " << rad2_[j].r << " " << cd.r <<  "\n";
                //e_rad2_rgb_[j] = emit_rgb_[j] + rad_rgb;// * reflex_factor;

            }


//            rad_.swap(rad2_);
            rad_ = rad2_;


        }

        //e_rad_rgb_.assign( e_rad_sse_.begin(), e_rad_sse_.end() );



    }


    aligned_buffer<col3f_sse> emit_;
    aligned_buffer<col3f_sse> rad_;
    aligned_buffer<col3f_sse> rad2_;
    const std::vector<std::vector<float> > &ffs_;
    const std::vector<std::vector<int> > &ff_target_;
    const std::vector<plane> &planes_;

};

class tick_timer {
public:
    tick_timer() :t1_(getticks()  ) {}

    double elapsed() {
        return ::elapsed( getticks(), t1_ );
    }

private:
    ticks t1_;
};

class rad_core_threaded: public rad_core {
    //typedef std::mutex lock_type;
    typedef spinlock_mutex lock_type;
public:
    class worker {

    };

    
    std::vector<std::pair<size_t,size_t>> calc_plane_distribution( const size_t num_partition ) {
        size_t num_ints = 0;
        std::vector<std::pair<size_t,size_t>> parts;
        for( auto &ff : ffs_ ) {
            num_ints += ff.size();
        }
        
        std::cout << "num_ints: " << num_ints << "\n";
        
        size_t ints_per_part = num_ints / num_partition;
        
        size_t first = 0;
        size_t last = 0;
        size_t acc = 0;
        for( size_t i = 0; i < num_partition; ++i ) {
            while( last < ffs_.size() && acc < ints_per_part ) {
                acc += ffs_[last].size();
                ++last;
            }
            acc = 0;
            parts.emplace_back( first, last );
            first = last;
        }
        
        parts.back().second = ffs_.size();
        
        
        for( auto &p : parts ) { 
            std::cout << "part: " << p.first << " " << p.second << "\n";
        }
        return parts;
    }

    rad_core_threaded( const std::vector<plane> &planes, const std::vector<std::vector<float> > &ffs, const std::vector<std::vector<int> > &ff_target )
            : rad_is_new_(false),
            emit_is_new_(false),
            emit_new_(planes.size()), emit_( planes.size() ), rad_( planes.size() ), rad2_( planes.size() ),
            ffs_(ffs), ff_target_(ff_target),
            planes_(planes),
            pints_(0),
            pints_last_(0),
            pints_last_time_(0)
    {

        
        
        if ( !true ) {
            threads_.push_back( std::thread( [&]() {
                work(0, planes_.size(), 0);
            }) );

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
        } else {
            const size_t num_planes = planes_.size();
            const size_t num_threads = 3;
            
            auto part = calc_plane_distribution(num_threads);
            for ( size_t i = 0; i < num_threads; ++i ) {
                
                
                
                //const size_t first = num_planes / num_threads * i;
                //const size_t last = num_planes / num_threads * (i+1);
                const size_t first = part.at(i).first;
                const size_t last = part.at(i).second;
            
                std::cout << "thread: " << i << " " << first << " " << last << " " << num_planes << "\n";
                
                threads_.push_back( std::thread( [=]() {
                    work( first, last, i );
                }));
            }

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
//          thread1_ = std::thread( [&]() { work(planes_.size()/4, planes_.size()/4 * 2);
//          });
//          thread2_ = std::thread( [&]() { work(planes_.size()/4 * 2, planes_.size()/4 * 3);
//          });
//          thread3_ = std::thread( [&]() { work(planes_.size()/4 * 3, planes_.size()/4 * 4);
//          });
        }
    }

    virtual void set_emit( const std::vector<vec3f> &emit ) {
        std::lock_guard<lock_type>lock(mtx_);

        assert( emit_new_.size() == emit.size() );

        emit_new_.assign( std::begin(emit), std::end(emit));
        emit_is_new_ = true;
    }

    virtual bool update() {
        return rad_is_new_;
//         std::lock_guard<std::mutex> lock(mtx_);
// 
//         if ( rad_is_new_ ) {
//             rad_is_new_ = false;
//             return true;
//         } else {
//             return false;
//         }


    }

    virtual void copy( std::vector<vec3f> *out ) {
        std::lock_guard<lock_type> lock(mtx_);
        for ( size_t i = 0; i < rad_.size(); ++i ) {
            (*out)[i].r = rad_[i].r;
            (*out)[i].g = rad_[i].g;
            (*out)[i].b = rad_[i].b;
        }

        rad_is_new_ = false;
        cl_ubyte64 time = CL_System::get_microseconds();
        
        cl_ubyte64 dt = time - pints_last_time_;
        
        if( dt >= 1e6 ) {
            const size_t dp = pints_ - pints_last_;
            pints_last_ = pints_;
            
            pints_last_time_ = time;
            
            std::cout << "pint/s: " << dp / (dt / 1.0e6) << "\n";
        }
    }




private:

    void work( size_t first, size_t last, size_t rank ) {
        while (true) {
            if( rank == 0 )
            {
                std::lock_guard<lock_type> lock(mtx_);

                if ( emit_is_new_ ) {
                    //emit_.swap(emit_new_);
                    //emit_ = std::move(emit_new_);
                    assert( emit_new_.size() == emit_.size());
                    std::copy( emit_new_.begin(), emit_new_.end(), emit_.begin() );
                    emit_is_new_ = false;
                }

            }


//             tick_timer tt;
            do_radiosity_sse(first, last);

//             std::cout << "elapsed: " << tt.elapsed() << "\n";

            if( rank == 0 )
            {
                std::lock_guard<lock_type> lock(mtx_);
                rad_is_new_ = true;
            }

        }
    }

    void do_radiosity_sse( size_t first, size_t last ) {
        const int steps = 1;
        const float min_ff = 0;
//      std::fill(e_rad_sse_.begin(), e_rad_sse_.end(), vec3f(0.0, 0.0, 0.0));
//      std::fill(e_rad2_sse_.begin(), e_rad2_sse_.end(), vec3f(0.0, 0.0, 0.0));

        //std::copy( emit_sse_.begin(), emit_sse_.end(), e_rad_sse_.begin() );

        typedef vector_unit<float,4> vu;
        typedef vu::vec_t vec_t;
        //steps = 0;
//         std::cout << "do: " << first << " " << last << "\n";
        vec_t reflex_factor = vu::set1(1.0);
        
        size_t pints = 0;
        for ( int i = 0; i < steps; ++i ) {
            //for( auto it = pairs_.begin(); it != pairs_.end(); ++it, ++ff_it ) {


            for ( size_t j = first; j < last; ++j ) {

                const size_t s = ffs_[j].size();
                vec_t rad = vu::set1(0);

                const vec3f cd = planes_[j].col_diff();
                const vec_t col_diff = vu::set( 0, cd.b, cd.g, cd.r );



                for ( size_t k = 0; k < s; ++k ) {
                    size_t target = ff_target_[j][k];
                    //rad_rgb += (col_diff * e_rad_rgb_[target]) * ff2s_[j][k];

                    if ( false && ffs_[j][k] < min_ff ) {
                        continue;
                    }

                    const vec_t ff = vu::set1( ffs_[j][k] );

                    rad = vu::add( rad, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target))), ff ));

                }
                pints += s;
                vu::store( vu::add( vu::load((float*)emit_(j)), vu::mul(rad, reflex_factor)), (float*)rad2_(j));
//              std::cout << "col: " << rad2_[j].r << " " << cd.r <<  "\n";
                //e_rad2_rgb_[j] = emit_rgb_[j] + rad_rgb;// * reflex_factor;

            }


            {
                std::lock_guard<lock_type>lock(mtx_);
                //rad_.swap(rad2_);
                std::copy( &rad2_[first], &rad2_[last], &rad_[first]);
                pints_ += pints;
            }
            //rad_ = rad2__;


        }

        //e_rad_rgb_.assign( e_rad_sse_.begin(), e_rad_sse_.end() );



    }

    lock_type mtx_;
//  std::thread thread0_;
//  std::thread thread1_;
//  std::thread thread2_;
//  std::thread thread3_;

    std::vector<std::thread> threads_;

    bool rad_is_new_;
    bool emit_is_new_;

    aligned_buffer<col3f_sse> emit_new_;
    aligned_buffer<col3f_sse> emit_;
    aligned_buffer<col3f_sse> rad_;
    aligned_buffer<col3f_sse> rad2_;
    const std::vector<std::vector<float> > &ffs_;
    const std::vector<std::vector<int> > &ff_target_;
    const std::vector<plane> &planes_;
    size_t pints_;
    size_t pints_last_;
    cl_ubyte64 pints_last_time_;
};



class light_scene {
public:
    light_scene( const std::vector<plane> &planes, const bitmap3d &solid )
            : planes_(planes), solid_(solid)
    {


        emit_rgb_.resize( planes_.size() );
        e_rad_rgb_.resize( planes_.size() );
        e_rad2_rgb_.resize( planes_.size() );

        emit_sse_.resize( planes_.size() );
        e_rad_sse_.resize( planes_.size() );
        e_rad2_sse_.resize( planes_.size() );

        glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if ( false ) {
            setup_formfactors();
            {
                std::ofstream os( "ff.bin" );


                write_formfactors(os);

            }
        } else {
            try
            {


                std::ifstream is( "ff.bin" );
                if ( !is.good() ) {
                    throw std::runtime_error( "cannot open ff.bin");
                }
                read_formfactors(is);

            } catch ( std::runtime_error x ) {
                std::cerr << x.what() << "\n";
                std::cerr << "error while reading ff.bin. regenerating\n";
                setup_formfactors();
            }
        }


        rad_core_ = make_unique<rad_core_threaded>(planes_, ff2s_, ff2_target_);

    }


    void reset_emit() {
        std::fill( emit_rgb_.begin(), emit_rgb_.end(), vec3f(0.0, 0.0, 0.0));
    }

    void render_light( const vec3f &light_pos, const vec3f &light_color ) {
        //std::fill( emit_rgb_.begin(), emit_rgb_.end(), vec3f(0.2, 0.2, 0.2 ));
        std::fill( emit_rgb_.begin(), emit_rgb_.end(), vec3f(0.0, 0.0, 0.0 ));
        for ( size_t i = 0; i < planes_.size(); ++i ) {

            auto &p = planes_[i];

            vec3f trace_pos = p.pos() + p.norm();

            const bool occ = true && util::occluded( light_pos, trace_pos, solid_ );


            if ( !occ ) {
                vec3f d = light_pos - p.pos();
                d.normalize();
                float len = d.length();
                d /= len;
                float dot = d.dot( p.norm() );

                if ( dot > 0 ) {
                    emit_rgb_[i] += (p.col_diff() * light_color) * dot * (5/(2*3.1415*len*len));
                }

            }
        }

        //emit_sse_.assign( emit_rgb_.begin(), emit_rgb_.end() );


    }

    void post() {
        std::copy( emit_rgb_.begin(), emit_rgb_.end(), emit_sse_.begin() );
    }

    float randf() const {
        return std::rand() / float(RAND_MAX);
    }

    void render_emit_patches() {
        static int x = 0;
        static int xd = 1;
        for ( size_t i = 0; i < planes_.size(); ++i ) {



//          if( randf() > 0.9 ) {
//              emit_rgb_[i] += vec3f(randf(), randf(), randf());
//          }
            if ( (planes_[i].pos().y / 2) == x ) {
                //if( planes_[i].dir() == plane::dir_zx_p ) {
                emit_rgb_[i] += planes_[i].col_diff();
            }


        }
        x += xd;
        if ( x > 15 || x == 0) {
            xd = -xd;
        }

    }


    void do_radiosity( int steps = 10,  float min_ff = 0.0 ) {
        // TODO: rename and/or remove parameters
        rad_core_->set_emit( emit_rgb_ );
//         rad_core_->update();
        rad_core_->copy( &e_rad_rgb_ );
        return;

    }



    const vec3f &rad_rgb( size_t i ) {
        return e_rad_rgb_[i];
    }

    const std::vector<vec3f> &rad_rgb() const {
        return e_rad_rgb_;
    }
     
    
private:
    inline bool normal_cull( const plane &pl1, const plane &pl2 ) {
        const plane::dir_type d1 = pl1.dir();
        const plane::dir_type d2 = pl2.dir();

        const vec3i &p1 = pl1.pos();
        const vec3i &p2 = pl2.pos();

        return p1 == p2 ||
               d1 == d2 ||
               (d1 == plane::dir_xy_n && d2 == plane::dir_xy_p && p1.z < p2.z) ||
               (d1 == plane::dir_xy_p && d2 == plane::dir_xy_n && p1.z > p2.z) ||
               (d1 == plane::dir_yz_n && d2 == plane::dir_yz_p && p1.x < p2.x) ||
               (d1 == plane::dir_yz_p && d2 == plane::dir_yz_n && p1.x > p2.x) ||
               (d1 == plane::dir_zx_n && d2 == plane::dir_zx_p && p1.y < p2.y) ||
               (d1 == plane::dir_zx_p && d2 == plane::dir_zx_n && p1.y > p2.y);

    }

    void write_formfactors( std::ostream &os ) {
        size_t size1 = ff2s_.size();

        os.write((char*) &size1, sizeof(size_t));
        for ( size_t i = 0; i < size1; ++i ) {
            size_t size2 = ff2s_[i].size();
            os.write((char*) &size2, sizeof(size_t));

            os.write( (char*) ff2s_[i].data(), size2 * sizeof(float));
            os.write( (char*) ff2_target_[i].data(), size2 * sizeof(int));
        }

    }

    void read_formfactors( std::istream &is ) {
        size_t size1;
        is.read( (char *) &size1, sizeof( size_t ));

        if ( size1 != planes_.size() ) {
            throw std::runtime_error( "cannot read form factors: size1 != planes_.size()");
        }


        std::vector<std::vector<float> > ff2s(size1);
        std::vector<std::vector<int> > ff2_target(size1);


        for ( size_t i = 0; i < size1; ++i ) {
            size_t size2;

            is.read( (char *) &size2, sizeof( size_t ));

            if ( size2 > planes_.size() ) {
                throw std::runtime_error( "cannot read form factors: size2 != planes_.size()");
            }

            ff2s[i].resize(size2);
            ff2_target[i].resize(size2);

            is.read( (char*) ff2s[i].data(), size2 * sizeof(float));
            is.read( (char*) ff2_target[i].data(), size2 * sizeof(int));

            std::for_each(ff2_target[i].begin(), ff2_target[i].end(), [&](int t) {
                if ( size_t(t) >= planes_.size() ) {
                    throw std::runtime_error( "cannot read form factors: target out of range");
                }
            });
        }

        // exception safe!
        ff2s_.swap( ff2s );
        ff2_target_.swap( ff2_target );

    }

    void setup_formfactors() {

//      std::ofstream os( "ff.txt" );

        std::ofstream os;//( "matrix.pnm");
        os << "P1\n";
        os << planes_.size() << " " << planes_.size() << "\n";

        ff2s_.resize(planes_.size());
        ff2_target_.resize(planes_.size());

        std::vector<float> ff_tmp;
        std::vector<int> target_tmp;
        size_t num_ff = 0;

        for ( size_t i = 0; i < planes_.size(); ++i ) {
            vec3i p1 = planes_.at(i).pos();


            vec3f norm1 = planes_[i].norm();
            vec3f p1f = p1;//(p1 + norm1 * 0.5);

//             size_t minj = size_t(-1);
//             size_t maxj = 0;



            ff_tmp.clear();
            target_tmp.clear();

            for ( size_t j = 0; j < i /*planes_.size()*/; ++j ) {

                //                  std::cerr << "i: " << i << " " << j << "\n";


                if ( normal_cull( planes_[i], planes_[j] )) {
//                  os << "0 ";
                    continue;
                }

//              if( planes_[i].normal_cull(planes_[j])) {
//                  os << "0 ";
//                  continue;
//              }

                const vec3i &p2 = planes_[j].pos();


                const vec3f &norm2 = planes_[j].norm();
                vec3f p2f = p2;// + (norm2 * 0.5);

                float d2 = dist_sqr( p1f, p2f );



                bool dist_cull = false;


                float ff = 0;
                {


                    vec3f dn = (p1f - p2f).normalize();
                    //std::cout << p1_3d << " " << p2_3d << " " << dn << "\n";
                    //.norm();
                    float ff1 = std::max( 0.0f, norm1.dot(vec3f(0.0, 0.0, 0.0)-dn));
                    float ff2 = std::max( 0.0f, norm2.dot(dn));

                    ff = ff1 * ff2;
//                  ff = std::max( 0.0f, ff );
                    //                      std::cout << "ff: " << ff << "\n";
                    //dist_cull = ff < 0.01;
                }



                //dist_cull = false;



                ff /=  (3.1415 * d2);

//              os << ff << "\n";

                dist_cull = ff < 5e-5;

                if ( !dist_cull && i != j ) {

                    if ( util::occluded( p1 + norm1, p2 + norm2, solid_ )) {
                        //                  os << "0 ";
                        continue;
                    }

//                  pairs_.push_back(std::make_pair(i,j));
//
//
//                  ffs_.push_back(ff);// / (3.1415 * d2));
                    ff_tmp.push_back(ff);
                    target_tmp.push_back(j);

                    // j < i => ff2s_[j] is already initialized.
                    ff2s_[j].push_back(ff);
                    ff2_target_[j].push_back(i);

//                  minj = std::min( minj, j );
//                  maxj = std::max( maxj, j );
                    ++num_ff;
//                  os << "1 ";
                } else {
                    os << "0 ";
                }

            }

            ff2s_[i].assign(ff_tmp.begin(), ff_tmp.end());
            ff2_target_[i].assign(target_tmp.begin(), target_tmp.end());
            os << "\n";

            //std::cout << "i: " << i << " " << num_ff << " " << minj << " - " << maxj << "\n";
        }






        std::cout << "num interactions (half): " << num_ff << "\n";

//      throw std::runtime_error("xxx");

//      e_emit.resize(patches_.size());
//      e_rad.resize(patches_.size());
//
//      e_emit_rgb.resize(patches_.size());
//
//
//      col_diff.resize(patches_.size());

    }

    const std::vector<plane> planes_;
    const bitmap3d solid_;


    std::vector<vec3f> emit_rgb_;
    std::vector<vec3f> e_rad_rgb_;
    std::vector<vec3f> e_rad2_rgb_;


    aligned_buffer<col3f_sse> emit_sse_;
    aligned_buffer<col3f_sse> e_rad_sse_;
    aligned_buffer<col3f_sse> e_rad2_sse_;



    std::vector<std::vector<float> > ff2s_;
    std::vector<std::vector<int> > ff2_target_;

    std::unique_ptr<rad_core> rad_core_;

};


class input_mapper {
public:

    void add_mapping( int keycode, bool *indicator ) {
        mappings_.push_back(mapping( keycode, indicator ));
    }

    void input( const CL_InputDevice &dev ) {
        for ( const mapping &m : mappings_ ) {
            *m.indicator_ = dev.get_keycode(m.keycode_);
        }
    }

private:

    struct mapping {
        int keycode_;
        bool *indicator_;

        mapping( int keycode, bool *indicator ) : keycode_(keycode), indicator_(indicator) {
            *indicator_ = false;
        }

    };

    std::vector<mapping> mappings_;

};

class mouse_mapper {
public:
    mouse_mapper() : valid_(false), delta_x_(0), delta_y_(0) {}

    void input( const CL_InputDevice &mouse ) {
        int new_x = mouse.get_x();
        int new_y = mouse.get_y();

        if ( valid_ ) {
            delta_x_ = new_x - old_x_;
            delta_y_ = new_y - old_y_;
        }

        old_x_ = new_x;
        old_y_ = new_y;
        valid_ = true;
    }

    float delta_x() const {
        return delta_x_;
    }

    float delta_y() const {
        return delta_y_;
    }

private:

    bool valid_;
    int old_x_;
    int old_y_;

    float delta_x_;
    float delta_y_;
};

class player {
public:

    player() {
        input_mapper_.add_mapping( CL_KEY_W, &i_forward_);
        input_mapper_.add_mapping( CL_KEY_S, &i_backward_);
        input_mapper_.add_mapping( CL_KEY_A, &i_left_);
        input_mapper_.add_mapping( CL_KEY_D, &i_right_);

        pos_.x = 0;
        pos_.y = 0;
        pos_.z = 0;

        rot_y_ = 0.0;
        rot_x_ = 0.0;
    }

    void input( const CL_InputDevice &keyboard, const CL_InputDevice &mouse ) {
        input_mapper_.input(keyboard);

        mouse_mapper_.input(mouse);
    }

    void frame( double time, double dt ) {

        
        // set up player head to world rotation from x/y rotations
        rot_y_ -= mouse_mapper_.delta_x(); // NOTE: rotations are CCW, so clockwise head rotation (to the right) means negative mouse delta
        rot_x_ -= mouse_mapper_.delta_y();

//          CL_Mat4f rot = CL_Mat4f::rotate(CL_Angle(rot_x_, cl_degrees), CL_Angle(rot_y_, cl_degrees), CL_Angle(), cl_YXZ );

        CL_Quaternionf rot_quat(CL_Angle(rot_x_, cl_degrees), CL_Angle(rot_y_, cl_degrees), CL_Angle(), cl_YXZ );

        

        CL_Vec4f trans_vec(0.0, 0.0, 0.0, 1.0);

        const float move_speed = 4.0; // 4 m/s

        if ( i_forward_ ) {
            trans_vec.z -= move_speed * dt; // NOTE: translation is in 'player head coordinate system' so forward is -z
        }
        if ( i_backward_ ) {
            trans_vec.z += move_speed * dt;
        }
        if ( i_left_ ) {
            trans_vec.x -= move_speed * dt;
        }
        if ( i_right_ ) {
            trans_vec.x += move_speed * dt;
        }

        // rotate translation vector from 'player head coordinates' into world coordinates
        
        trans_vec = rot_quat.rotate_vector(trans_vec);
        
        pos_.x += trans_vec.x;
        pos_.y += trans_vec.y;
        pos_.z += trans_vec.z;
        
        

    }
    const vec3f &pos() const {
        return pos_;
    }
    float rot_x() const {
        return rot_x_;
    }
    float rot_y() const {
        return rot_y_;
    }

private:
    bool i_forward_;
    bool i_backward_;
    bool i_left_;
    bool i_right_;

    input_mapper input_mapper_;
    mouse_mapper mouse_mapper_;

    vec3f pos_;
    float rot_y_;
    float rot_x_;

};

class vbo_builder {
  
    
public:
    vbo_builder( size_t num_planes ) : num_planes_(num_planes) {
        GLuint b;
        glGenBuffers( 2, buffers_ );
        glGenBuffers( 1, &index_buffer_ );
        
        const size_t vertex_size = num_planes_ * (4 * 3) * sizeof(float);
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] );
        glBufferData( GL_ARRAY_BUFFER, vertex_size, 0, GL_STATIC_DRAW );
        
        color_size = num_planes_ * (4 * 4) * sizeof(GLubyte);
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] );
        glBufferData( GL_ARRAY_BUFFER, color_size, 0, GL_DYNAMIC_DRAW );
        
        glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_buffer_ );
        glBufferData( GL_ELEMENT_ARRAY_BUFFER, num_planes_ * 4 * sizeof(int), 0, GL_STATIC_DRAW );
        
        
    }
    
    void update_index_buffer( size_t num_planes ) {
        glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_buffer_ );
        GLuint *b = (GLuint *) glMapBuffer( GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY );
        assert( b != nullptr );
        for( size_t i = 0; i < num_planes * 4; ++i ) {
            b[i] = i;
        }
        
        glUnmapBuffer( GL_ELEMENT_ARRAY_BUFFER );
        
    }
    
    template<typename iiter>
    void update_vertices( iiter first, iiter last ) {
        //const size_t num_planes_ = std::distance( first, last );
        
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] );
        float *b = (float*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY );
        assert( b != nullptr );
        for( ; first != last; ++first ) {
         
            const plane &p = *first;
            //size_t s1 = vertex_.size();
            
            p.render_vertex( b );
            b += 4 * 3;
            
        }
        
        glUnmapBuffer( GL_ARRAY_BUFFER );
    }
    
    template<const int minv, const int maxv >
    inline GLubyte clamp( int v ) {
        //return std::max( minv, std::min( maxv, v ));
        return std::min( maxv, v );
        //return v;
    }
    
    template<typename iiter>
    void update_color( iiter first, iiter last ) {
//         const size_t num_planes_ = std::distance( first, last );
        
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] );
        GLubyte *b = (GLubyte*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY );
        assert( b != nullptr );
        
        //b += 4 * 3 * num_planes_;
        
        for( ; first != last; ++first ) {
            for( size_t i = 0; i < 4; ++i ) {
                *(b++) = clamp<0,255>(255 * first->r);
                *(b++) = clamp<0,255>(255 * first->g);
                *(b++) = clamp<0,255>(255 * first->b);
                *(b++) = 255;
            }   
            
        }
        
        glUnmapBuffer( GL_ARRAY_BUFFER );
    }
    
    void draw_arrays() {
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] );
        glVertexPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL));
       // glColorPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL+ 4 * 3 * num_planes_ * sizeof(GLfloat) ));
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] );
        glColorPointer(4, GL_UNSIGNED_BYTE, 0, (GLvoid*)((char*)NULL));
        
//         glBindBuffer(GL_ARRAY_BUFFER, buffer_);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
// This is the actual draw command
        glDrawElements(GL_QUADS, num_planes_ * 4, GL_UNSIGNED_INT, (GLvoid*)((char*)NULL));
    }
//private:
    size_t color_size;
    
    GLuint buffers_[2];
    GLuint index_buffer_;
    
    const size_t num_planes_;
};

class vertex_array_builder {
public:
    vertex_array_builder() {
        
    }

    template<typename iiter>
    void render( iiter first, iiter last ) {
        num_planes_ = std::distance( first, last );
        
        for( ; first != last; ++first ) {
         
            const plane &p = *first;
            size_t s1 = vertex_.size();
            
            p.render_color( std::back_inserter(color_));
            p.render_vertex( std::back_inserter(vertex_));
            
            size_t s2 =vertex_.size();
            
            for( ;s1 != s2; ++s1 ){
                index_.push_back(s1);
            }
        }
        
        std::cout << "vab: " << vertex_.size() << " " << index_.size() << "\n";
    }
    
    template<typename iiter>
    void update_color( iiter first, iiter last ) {
        auto cfirst = color_.begin();
        
        assert( std::distance( first, last ) == num_planes_ );
        
        for( ;first != last; ++first ) {

            for( size_t i = 0; i < 4; ++i ) {
                *(cfirst++) = first->r;
                *(cfirst++) = first->g;
                *(cfirst++) = first->b;
            }
            
//             const plane &p = *first;
//  
//             cfirst = p.render_color( cfirst );

        
        
            assert( cfirst <= color_.end() );
        }
        
        
        
    }
    
    void setup_gl_pointers() {
        glVertexPointer(3, GL_FLOAT, 0, vertex_.data());
        glColorPointer(3, GL_FLOAT, 0, color_.data() );
        glIndexPointer(GL_INT, 0, index_.data());
        
        glEnableClientState( GL_VERTEX_ARRAY );
        glEnableClientState( GL_COLOR_ARRAY );
        glEnableClientState( GL_INDEX_ARRAY );
    }
    
    
    void draw_arrays() {
        glDrawArrays(GL_QUADS, 0, num_planes_ * 4);
    }
    
private:
    std::vector<float> vertex_;
    std::vector<float> color_;
    std::vector<int> index_;
    
    size_t num_planes_;
};

class ortho {

public:

    ublas::matrix<int> pump( const ublas::matrix<int> &in, const size_t factor ) {
        ublas::matrix<int> out( in.size1() * factor, in.size2() * factor );
        
        for( size_t row = 0; row != in.size1(); ++row ) {
            ublas::matrix_row<ublas::matrix<int>> r( out, row * factor );
            
           
            
            for( size_t col = 0; col != in.size2(); ++col ) {
                for( size_t i = 0; i < factor; ++i ) {
                    r[col*factor+i] = in(row, col);
                }
            }
            for( size_t i = 0; i < factor; ++i ) {
                ublas::matrix_row<ublas::matrix<int>> r2( out, row * factor + i);
                std::copy( r.begin(), r.end(), r2.begin() );
            }
            
        }
        
        return out;
    }
    



    std::vector<std::vector<int>> matrix_to_intvec2d( const ublas::matrix<int> &in ) {
        std::vector<std::vector<int>> out;
        out.reserve(in.size1());
        
        for( auto it1 = in.begin1(); it1 != in.end1(); ++it1 ) {
            out.emplace_back( it1.begin(), it1.end() );
            
//             std::cout << "len2: " << out.back().size() << "\n";
        }
//         std::cout << "len1: " << out.size() << "\n";
        return out;
    }
    
    ublas::matrix<int> load_crystal_slice( std::istream &is, size_t width, size_t height ) {
        ublas::matrix<int> slice( height, width );
        
        auto mapc = [](char c) {
                c = std::tolower(c);

                if ( c == ' ' ) {
                    return int(0);
                } else if ( c >= 'a' && c <= 'z' ) {
                    return int(1 + c - 'a');
                } else if ( c >= '0' && c <= '9' ) {
                    return int(2 + 'z' - 'a' + c - '0');
                } else {
                    std::cerr << "bad: " << int(c) << "\n";
                    throw std::runtime_error( "bad character in map");
                }
        };
        
        
        auto it1 = slice.begin1();
        for( size_t i = 0; i < height; ++i, ++it1 ) {
            auto it2 = it1.begin();
            for( size_t j = 0; j < width; ++j, ++it2 ) {
                *it2 = mapc( is.get() );
            }
            
            int nl = is.get();
            assert( nl == '\n' );
        }
        
        return slice;
    }
    
    
    std::vector<ublas::matrix<int> > load_crystal( std::istream &is ) {
        size_t width;
        size_t height;
        size_t num;
        
        is >> width;
        is >> height;
        is >> num;
        
        while( is.get() != '\n' ) {}
        
        std::cout << "size: " << width << " " << height << " " << num << "\n";
        
        std::vector<ublas::matrix<int> > out;
        
        for( size_t i = 0; i < num; ++i ) {
            
            out.emplace_back(pump(load_crystal_slice( is, width, height ), pump_factor_));

            std::cout << "map: " << i << "\n";
            
            for( auto it1 = out.back().begin1(); it1 != out.back().end1(); ++it1 ) {
                std::copy( it1.begin(), it1.end(), std::ostream_iterator<int>(std::cout, " " ));
                std::cout << "\n";
            }

            
        }
        
        
        std::cout << "size: " << out.size() << "\n";
        
        return out;
        
        //throw "exit";
        
//         size_t len = size_t(-1);
//         while ( !is.eof() ) {
//             std::string line;
// 
//             std::getline(is, line);
// 
// //          while( !is.eof() ) {
// //              char c = is.get();
// //              if( c == '\n' ) {
// //                  break;
// //              }
// //              line.push_back(c);
// //          }
// 
//             if ( line.empty()) {
//                 break;
//             }
// 
//             std::cout << "len: " << line.size() << "'" << std::string(line.begin(), line.end()) << "'\n";
// 
// 
//             if ( len == size_t(-1)) {
//                 len = line.size();
//             } else {
//                 assert( len == line.size() );
//             }
// 
//             ret.push_back(std::vector<int>(len));
// 
//             std::transform( line.begin(), line.end(), ret.back().begin(), [](char c) {
//                 c = std::tolower(c);
// 
//                 if ( c == ' ' ) {
//                     return int(0);
//                 } else if ( c >= 'a' && c <= 'z' ) {
//                     return int(1 + c - 'a');
//                 } else if ( c >= '0' && c <= '9' ) {
//                     return int(2 + 'z' - 'a' + c - '0');
//                 } else {
//                     std::cerr << "bad: " << int(c) << "\n";
//                     throw std::runtime_error( "bad character in map");
//                 }
//             });
// 
//         }



        
    }

    ortho() :
      pump_factor_(4) 
    {

//         glVertexPointer(4, GL_FLOAT, 0, ptr);
        //glEnableClientState( GL_VERTEX_ARRAY );
//         glColorPointer(

        std::ifstream is( "cryistal-castle-hidden-ramp.txt" );
//         std::ifstream is( "house1.txt" );
        //std::ifstream is( "cryistal-castle-tree-wave.txt" );

        assert( is.good() );
        height_fields_ = load_crystal(is);
        std::cout << "hf: " << height_fields_.size() << "\n";
        
        init_solid(height_fields_);
        
        

        init_planes();

        CL_OpenGLWindowDescription desc;
        desc.set_size( CL_Size( 1024, 768 ), true );

        desc.set_depth_size(16);
        //std::cout << "depth: " << desc.get_depth_size();
        
        wnd_ = CL_DisplayWindow(desc);

        CL_GraphicContext_GL gc = wnd_.get_gc();
//      //CL_Mat4f proj = CL_Mat4f::ortho( 0, 1024, 0, 768, 100, -100 );


        gc.set_active();
        
        try {
            init_opencl();
        } catch( cl::Error x ) {
            
//             std::array<void*, 256> bt;
//             //void *bt[256];
//             
//             size_t size = backtrace( bt.data(), bt.size() );
//             char **strings = backtrace_symbols( bt.data(), size );
//             std::cout << "backtrace: " << size << "\n";
//             for( size_t i = 0; i < size; ++i ) {
//                 std::cout << i << " " << strings[i] << "\n";
//             }
//             free( strings );
            
            std::cerr << "opencl initialization failed\ncall: " << x.what() << "\nerror code: " << cl_str_error( x.err() ) << "\n";            
            throw;
        }
        
      //  throw 0;
        
        glMatrixMode(GL_PROJECTION);                        //hello

        
        

        CL_Mat4f proj = CL_Mat4f::perspective( 60, 1.5, 2, 200 );
//      CL_Mat4f proj = CL_Mat4f::ortho( -20.0 * pump_factor_, 20.0 * pump_factor_, -15.0 * pump_factor_, 15.0 * pump_factor_, 0, 200 );
        //CL_Mat4f proj = CL_Mat4f::ortho( -40, 40, -30, 30, 0, 200 );


        glLoadMatrixf( proj.matrix );


        //CL_Texture tex(gc, 64, 64 );


        struct texel {
            GLubyte col[3];
            GLubyte alpha;

            texel() {
                col[0] = 128;
                col[1] = 128;
                col[2] = 128;
                alpha = 255;

            }
        };

        std::array<texel,64 * 64> tex_data;

        glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 0.0);
//      glGenTextures(1, &texName);
//      glBindTexture(GL_TEXTURE_2D, texName);
//      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 64, 64, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex_data.data());


//      gc.set_map_mode(cl_user_projection);
//      gc.set_projection(proj);
//
//      gc.flush_batcher();
//      glMatrixMode(GL_PROJECTION);


        glMatrixMode(GL_MODELVIEW);


        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDepthFunc(GL_LESS);
    }

    void setup_perspective( const player &camera ) {
        glMatrixMode(GL_PROJECTION);                        //hello

        {
            CL_Mat4f proj = CL_Mat4f::perspective( 60, 1.5, 0.2, 50 );
            //      CL_Mat4f proj = CL_Mat4f::ortho( -20.0 * pump_factor_, 20.0 * pump_factor_, -15.0 * pump_factor_, 15.0 * pump_factor_, 0, 200 );
            //CL_Mat4f proj = CL_Mat4f::ortho( -40, 40, -30, 30, 0, 200 );


            glLoadMatrixf( proj.matrix );
        }

        //          std::cout << "pos: " << player_pos << "\n";

        glMatrixMode( GL_MODELVIEW );
        {
            const vec3f &player_pos = camera.pos();
            CL_Mat4f proj = CL_Mat4f::translate(-player_pos.x, -player_pos.y, -player_pos.z) * CL_Mat4f::rotate(CL_Angle(-camera.rot_x(), cl_degrees),CL_Angle(-camera.rot_y(),cl_degrees),CL_Angle(), cl_XYZ);
            glLoadMatrixf(proj.matrix);
        }
    }

    void setup_ortho() {
        glMatrixMode(GL_PROJECTION);                        //hello

        {

            CL_Mat4f proj = CL_Mat4f::ortho( -20.0, 20.0, -15.0, 15.0, 0, 200 );
            //CL_Mat4f proj = CL_Mat4f::ortho( -40, 40, -30, 30, 0, 200 );


            glLoadMatrixf( proj.matrix );
        }

        //          std::cout << "pos: " << player_pos << "\n";

        glMatrixMode( GL_MODELVIEW );
        {

            CL_Mat4f proj = CL_Mat4f::look_at( 10, 10, 10, 0, 0, 0, 0.0, 1.0, 0.0 );
//          CL_Mat4f proj = CL_Mat4f::translate(-player_pos.x, -player_pos.y, -player_pos.z) * CL_Mat4f::rotate(CL_Angle(-p1.rot_x(), cl_degrees),CL_Angle(-p1.rot_y(),cl_degrees),CL_Angle(), cl_XYZ);
            glLoadMatrixf(proj.matrix);
        }

    }

    void render_quads() {
        assert(false && "out of order" );
        glBegin(GL_QUADS);
        for ( size_t i = 0; i < planes_.size(); ++i ) {
            plane &p = planes_[i];

//          p.energy_rgb(ls.rad_rgb(i));

            p.render();
        }
        glEnd();

    }

//     vertex_array_builder setup_quad_buffers() {
//         vertex_array_builder vab;
//         
//         
//         
//         for ( size_t i = 0; i < planes_.size(); ++i ) {
//            
//             vab.render( planes_[i] );
//             
//         }
//         
//         return vab;
//     }
    
    void start() {



        float x1 = -10;
        float y1 = -10;




//      glBindTexture(GL_TEXTURE_2D, texName);
//         float roty = -45;
        int light_x = 0;
        int light_xd = 1;

        vec3i light_pos( 0, 40, 0 );

        int steps = 1;

        light_scene ls( planes_, solid_ );


        auto t_old = CL_System::get_microseconds();
        bool light_on = true;
        bool light_button_down = false;
        double delta_t = 0.01;
        player p1;

//      float min_ff = 5e-5;

        vertex_array_builder vab;
//         vab.render(planes_.begin(), planes_.end() );
        
        vbo_builder vbob(planes_.size());
        vbob.update_index_buffer(planes_.size());
        vbob.update_vertices( planes_.begin(), planes_.end());
        
        cl::Buffer buf;
        try {
            cl_int cl_err;
            buf = cl::Buffer( clCreateFromGLBuffer( cl_context_(), CL_MEM_WRITE_ONLY, vbob.buffers_[1], &cl_err ));
            assert( cl_err == CL_SUCCESS );
        
        
            cl_fcolor_ = cl::Buffer( cl_context_, CL_MEM_READ_ONLY, planes_.size() * 3 * sizeof(float) );
            
            cl_kernel_.setArg(0, buf() );
            cl_kernel_.setArg(1, cl_fcolor_ );
            cl_uint cl_color_size = planes_.size();
            cl_kernel_.setArg(2, cl_color_size );
            
            
            
        } catch( cl::Error x ) {
            std::cerr << "cl error during gl buffer setup\ncall: " << x.what() << "\nerror code: " << cl_str_error( x.err() ) << "\n";            
            throw;
        }
        bool light_changed = true;
        
        
        
        while ( true ) {

            //cube c(x1, 0, y1);

            CL_GraphicContext gc = wnd_.get_gc();
            CL_InputDevice &keyboard = wnd_.get_ic().get_keyboard();
            CL_InputDevice &mouse = wnd_.get_ic().get_mouse();

//          if( keyboard.get_keycode(CL_KEY_I) ) {
//          //  roty += 1;
//
//              min_ff += 5e-6;
//          }
//          if(  keyboard.get_keycode(CL_KEY_K) ) {
//              //roty -= 1;
//              min_ff -= 5e-6;
//          }
//
//          std::cout << "min ff: " << min_ff << "\n";

            p1.input(keyboard, mouse);
            p1.frame(t_old * 1.0e-3, delta_t);

            if ( keyboard.get_keycode(CL_KEY_LEFT) ) {
                light_pos.x += 1;
                light_changed = true;
            }
            if (  keyboard.get_keycode(CL_KEY_RIGHT) ) {
                light_pos.x -= 1;
                light_changed = true;
            }
            if ( keyboard.get_keycode(CL_KEY_UP) ) {
                light_pos.z += 1;
                light_changed = true;
            }
            if (  keyboard.get_keycode(CL_KEY_DOWN) ) {
                light_pos.z -= 1;
                light_changed = true;
            }
            if ( keyboard.get_keycode(CL_KEY_L )) {
                if ( !light_button_down ) {
                    light_on = !light_on;
                    light_changed = true;
                }
                light_button_down = true;
            } else {
                light_button_down = false;
            }



            glEnable(GL_CULL_FACE);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


//          glTranslatef(0,0,-100);
//          glRotatef(45,1,0,0);
//          glRotatef(roty,0,1,0);


            if ( light_x > 60 || light_x < -60 ) {
                light_xd = -light_xd;
//              light_x = -20;
            }

            //          light_planes(vec3i(light_x, 10, 10 ));
            
            if( light_changed ) {
                ls.reset_emit();
                //ls.render_light(vec3f( 40, 50.0, light_x ), vec3f(1.0, 1.0, 1.0 ));
                if ( light_on ) {
                    //ls.render_light(light_pos, vec3f(1.0, 0.8, 0.6 ));
                    
                    //vec3f light_pos( p1.pos)
                    
                    //              vec3f light_pos = (p1.pos()* pump_factor_) - base_pos_;
                    vec3f light_weird = (light_pos * pump_factor_) - base_pos_;
                    ls.render_light( light_weird, vec3f(1.0, 0.8, 0.6 ));
                }
                ls.post();
                light_changed = false;
            }
            

            //ls.render_emit_patches();

            steps = 1;
            ls.do_radiosity( steps );

             // stupid: transfer rgb energy fomr light scene to planes
            for ( size_t i = 0; i < planes_.size(); ++i ) {
                plane &p = planes_[i];
                p.energy_rgb(ls.rad_rgb(i));
            }
            
            light_x += light_xd;


//          glPushMatrix();
            
            //CL_Mat4f proj = CL_Mat4f::look_at( 20, 20, 20, 0, 0, 0, 0.0, 1.0, 0.0 );

            //vab.update_color( planes_.begin(), planes_.end() );
//             vab.update_color( ls.rad_rgb().begin(), ls.rad_rgb().end() );
            
//             vab.setup_gl_pointers();
     //       vbob.update_color(ls.rad_rgb().begin(), ls.rad_rgb().end());
#if 1
            try
            {
                glFinish();
                const auto &rad_colors = ls.rad_rgb();
                
                cl_cqueue_.enqueueWriteBuffer( cl_fcolor_, false, 0, rad_colors.size() * sizeof(vec3f), rad_colors.data() );
                
                cl_int err = clEnqueueAcquireGLObjects( cl_cqueue_(), 1, &(buf()), 0, 0, 0 );
                assert( err == CL_SUCCESS );


                // Set arg 3 and execute the kernel
                
                cl_cqueue_.enqueueNDRangeKernel( cl_kernel_, 0, cl::NDRange(planes_.size()) );
                
                               

                // unmap buffer object
                err = clEnqueueReleaseGLObjects( cl_cqueue_(), 1, &(buf()), 0,0,0);
                assert( err == CL_SUCCESS );
                
                cl_cqueue_.finish();
                
                
                
            } catch( cl::Error x ) {
                std::cerr << "cl error during kernel execution\ncall: " << x.what() << "\nerror code: " << cl_str_error( x.err() ) << "\n";               
                throw;
            }
#endif       
            
            //setup_perspective( p1 );
            setup_ortho();

            int ortho_width = 320;
            int ortho_heigth = 200;
            glViewport( gc.get_width() - ortho_width, gc.get_height() - ortho_heigth, ortho_width, ortho_heigth );
       //     vab.draw_arrays();
            //  render_quads();
//          glPopMatrix();
            setup_perspective(p1);
            glViewport( 0, 0, gc.get_width(), gc.get_height());
//             render_quads();
//            vab.setup_gl_pointers();
//             vab.draw_arrays();
            vbob.draw_arrays();

            wnd_.flip(1);



            auto t = CL_System::get_microseconds();

//            std::cout << "fps: " << 1e6 / (t - t_old) << "\n";
            delta_t = (t - t_old) * 1.0e-6;
            t_old = t;



//           std::cout << "delta: " << delta_t << "\n";

//             CL_System::sleep( 1000 / 60. );
//          getchar();
            CL_KeepAlive::process();


            x1 += 1;
            if ( x1 > 10 ) {
                x1 = -10;
                y1 += 1;
            }
            if ( y1 > 10 ) {
                y1 = -10;
            }

        }

    }

    static int main( const std::vector<CL_String> &args ) {
	 _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
//      plane p( plane::dir_zx_p, vec3f( 0.5, 0.5, 0.5 ));

//      return 0;
        ortho o;
        o.start();

        return 0;
    }

private:

    void init_solid( const std::vector<ublas::matrix<int>> &slices ) {
        const auto &slice0 = slices.at(0);
        size_t size_z = slice0.size1();
        size_t size_x = slice0.size2();
        
        size_t size_y = *std::max_element( slice0.data().begin(), slice0.data().end() ) + 1;
        
        std::cout << "size: " << size_x << " " << size_z << " " << size_y << "\n";
        
        solid_ = bitmap3d( size_x, size_y, size_z );
        
        bool add = true;
        for( auto &slice : slices ) {
        
            int starty = 0;
            
            // kind of hack: do not touch blocks below the 1-level in subtractive passes,
            // to prevent drilling the ground plane...
            if( !add ) {
                starty = 1;
            }
            
            for ( int y = starty; y < int(size_y); ++y ) {
                for ( size_t z = 0; z < size_z; ++z ) {
                    for ( size_t x = 0; x < size_x; ++x ) {
                        int h = slice(z,x);

                        
                        if ( h >= y ) {
                            solid_(x, y, z) = add;
                        }    
//                         solid_(x, y, z) = true;
                        
                        
                    }
                }
            }
            
            add = !add;
        }
    }
        
        


    void init_planes() {

//         size_t size_z = height_.size();
//         size_t size_x = height_.front().size();
//         int size_y = 0;
// 
//         for ( auto it = height_.begin(); it != height_.end(); ++it ) {
//             size_y = std::max( size_y, (*std::max_element(it->begin(), it->end())) + 1 );
//         }
// 
//         std::cout << "solid: " << size_x << " " << size_y << " " << size_z << "\n";
// 
//         solid_ = bitmap3d( size_x, size_y, size_z );
// 
// 
//         //vec3f base_pos;
//         
// 
// 
//         for ( int y = 0; y < size_y; ++y ) {
//             for ( size_t z = 0; z < height_.size(); ++z ) {
//                 for ( size_t x = 0; x < height_[z].size(); ++x ) {
//                     int h = height_[z][x];
// 
// 
//                     if ( h >= y ) {
//                         solid_(x, y, z) = true;
//                     }
// 
// 
//                 }
//             }
//         }
        
        //vec3f base_pos( -10.5, -10.5, -10.5);
        vec3f base_pos( -(solid_.x() / 2.0 + 0.5), -10.5, -(solid_.z() / 2.0 + 0.5));
        base_pos_ = base_pos;
        
//         std::cout << "base pos: " << base_pos_ << " " << solid_.x() << "\n";
        
        const auto &solidc = solid_;

        vec3i light_pos( 10, 10, 10 );
        float scale = 1.0 / pump_factor_;

        for ( int y = 0; y < int(solid_.y()); ++y ) {
            for ( size_t z = 0; z < solid_.z(); ++z ) {
                for ( size_t x = 0; x < solid_.x(); ++x ) {
                    if ( solid_(x, y, z) ) {

                        const bool occ = util::occluded( light_pos, vec3i(x,y,z), solid_ );

                        float energy = occ ? 0.2 : 1.0;



                        if ( !solidc(x,y,z+1)) {
                            planes_.push_back( plane( plane::dir_xy_p, base_pos, vec3i( x, y, z ), scale, energy));
                        }
                        if ( !solidc(x,y,z-1)) {
                            planes_.push_back( plane( plane::dir_xy_n, base_pos, vec3i( x, y, z ), scale, energy));
                        }
                        if ( !solidc(x+1,y,z)) {
                            planes_.push_back( plane( plane::dir_yz_p, base_pos, vec3i( x, y, z ), scale, energy));
                        }
                        if ( !solidc(x-1,y,z)) {
                            planes_.push_back( plane( plane::dir_yz_n, base_pos, vec3i( x, y, z ), scale, energy));
                        }
                        if ( !solidc(x,y+1,z)) {
                            planes_.push_back( plane( plane::dir_zx_p, base_pos, vec3i( x, y, z ), scale, energy));
                        }
                        if ( y > 0 && !solidc(x,y-1,z)) {
                            planes_.push_back( plane( plane::dir_zx_n, base_pos, vec3i( x, y, z ), scale, energy));
                        }


                    }
                }
            }
        }
//      for( size_t z = 0; z < height_.size(); ++z ) {
//          for( size_t x = 0; x < height_[z].size(); ++x ) {
//              planes_.push_back( plane( plane::dir_zx_n, base_pos, vec3i( x, 20, z ), 0.0));
//          }
//      }


        std::cout << "planes: " << planes_.size() << "\n";




    }


    void light_planes( const vec3i &light_pos ) {
        for ( plane &p : planes_ ) {

            vec3f trace_pos = p.pos() + p.norm();

            const bool occ = true && util::occluded( light_pos, trace_pos, solid_ );


            if ( !occ ) {
                vec3f d = light_pos - p.pos();
                d.normalize();
                float len = d.length();
                d /= len;
                float dot = d.dot( p.norm() );

                if ( dot > 0 ) {
                    p.energy( dot * (10/(2*3.1415*len*len)));
                } else {
                    p.energy(0.1);
                }

            } else {
                p.energy( 0.1 );
            }


        }
    }

    
    void init_opencl() {
        
        // get opencl platform list.
        std::vector<cl::Platform> platforms;
        cl_int err = cl::Platform::get( &platforms );
        
        assert( err == CL_SUCCESS );
        
        // find NVIDIA platform
        
        for( cl::Platform & p : platforms ) {
            std::string vendor = p.getInfo<CL_PLATFORM_VENDOR>();
            std::cout << "vendor: " << vendor << "\n";
            
            if( vendor.find( "NVIDIA" ) == 0 ) {
                cl_platform_ = p;
                std::cout << "selected\n";
            }
            
        }
        
        assert( cl_platform_() != 0 );

        
        // get device list for platform
        std::vector<cl::Device> devices;
        err = cl_platform_.getDevices( CL_DEVICE_TYPE_GPU, &devices );
        
        assert( err == CL_SUCCESS );
        
        // choose first device that supports gl sharing
        for( cl::Device & d : devices ) {
            auto name = d.getInfo<CL_DEVICE_NAME>();
            std::cout << "device: " << name << "\n";
            std::string extensions = d.getInfo<CL_DEVICE_EXTENSIONS>();
            const char * gl_sharing_ext = "cl_khr_gl_sharing";
            
            if( cl_used_devices_.empty() && extensions.find( gl_sharing_ext ) != extensions.npos ) {
                std::cout << "using\n";
                cl_used_devices_.push_back( d );
            }
        }
        
        assert( cl_used_devices_.size() == 1 );

        // create context
        cl_context_properties props[] =
        {
            CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
            CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
            CL_CONTEXT_PLATFORM, (cl_context_properties)cl_platform_(),
            0
        };
        
        cl_context_ = cl::Context( cl_used_devices_, props );
        assert( cl_context_() != 0 );
        auto num_devs = cl_context_.getInfo<CL_CONTEXT_NUM_DEVICES>();
        
        std::cout << "num devs: " << num_devs << "\n";
        
        cl_cqueue_ = cl::CommandQueue( cl_context_, cl_used_devices_.front(), 0, &err );
        assert( err == CL_SUCCESS );
        
                //cxGPUContext = clCreateContext(props, 1, &cdDevices[uiDeviceUsed], NULL, NULL, &ciErrNum);

        std::string source;
        {
            std::ifstream is( "test.cl" );
            
            // I hate myself for doing it this way, but it's the easiest...
            
            while( is.good() ) {
                
                char c = is.get();
                if( is.good() ) {
                    source.push_back(c);
                }
            }
        }
        
       // std::cout << "source >>>>" << source << "<<<<\n";
        assert( !source.empty() );
        
        try {
            cl_program_ = cl::Program( cl_context_, source, false, &err );
            
            cl_program_.build();
        } catch( cl::Error x ) {
            std::cout << "cl error during build: " << x.what() << " " << cl_str_error(x.err()) << "\n"; 
            
            if( x.err() == CL_BUILD_PROGRAM_FAILURE ) {
                auto log = cl_program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>( cl_used_devices_.front() );
                
                std::cout << "build error. log: \n" << log << "end of log.\n";
            }
            
            throw;
        }
        assert( err == CL_SUCCESS );
        
        cl_kernel_ = cl::Kernel( cl_program_, "convert_colors", &err );
        assert( err == CL_SUCCESS );
        
        
        
        
    }
    
    CL_SetupCore setup_core_;




    CL_SetupDisplay display_;
    //CL_SetupSWRender swrender;

    CL_SetupGL setup_gl_;
    CL_DisplayWindow wnd_;
    GLuint texName;

    
    std::vector<ublas::matrix<int>> height_fields_;
    
    std::vector<plane> planes_;
    bitmap3d solid_;

    const size_t pump_factor_;
    vec3f base_pos_;
    
    cl::Platform cl_platform_;
    std::vector<cl::Device> cl_used_devices_;
    cl::Context cl_context_;
    
    cl::Program cl_program_;
    cl::Kernel cl_kernel_;
    
    cl::CommandQueue cl_cqueue_;
    cl::Buffer cl_fcolor_;
};



CL_ClanApplication app(&ortho::main);
