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

#include <fstream>
#include <iomanip>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "misc_utils.h"
#include "scene_bits.h"

const uint32_t scene_static::restart_idx = 0x1FFFFFFF;

bool util::occluded2(vec3i p0, vec3i p1, const bitmap3d& solid) {
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
bool util::occluded(vec3i p0, vec3i p1, const bitmap3d& solid) {
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
vec3i plane::normali(plane::dir_type dt) {
   switch ( dt ) {
    case dir_xy_p:
        return vec3i( 0, 0, 1 );
    case dir_xy_n:
        return vec3i( 0, 0,-1 );
    case dir_yz_p:
        return vec3i( 1, 0, 0 );
    case dir_yz_n:
        return vec3i(-1, 0, 0 );
    case dir_zx_p:
        return vec3i( 0, 1, 0 );
    case dir_zx_n:
    default:
        return vec3i( 0,-1, 0 );
    }
    
}

vec3f plane::normal(plane::dir_type dt) {
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
vec3f plane::primary(plane::dir_type dt) {
    switch ( dt ) {
    case dir_xy_p:
    case dir_xy_n:
        return vec3f( 1.0, 1.0, 0.0 );
    case dir_yz_p:
    case dir_yz_n:
        return vec3f( 0.0, 1.0, 1.0 );
    case dir_zx_p:
    case dir_zx_n:
    default:
        return vec3f( 1.0, 0.0, 1.0 );
    }
}
vec3f plane::primary0(plane::dir_type dt) {
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
vec3f plane::primary1(plane::dir_type dt) {
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
std::array< float, 4 > plane::vgen0(plane::dir_type dt) {
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
std::array< float, 4 > plane::vgen1(plane::dir_type dt) {
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
vec3f plane::col_diff(plane::dir_type dt) {
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

plane::plane(plane::dir_type d, const vec3f& base_pos, const vec3i& pos, float scale, float energy) : dir_(d),
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

class face_iterator {
public:
    face_iterator( const bitmap3d &bm, plane::dir_type dir ) 
    : 
    x_(0), y_(0), z_(0),
    bitmap_(bm), dir_(dir), norm_( plane::normali(dir)), last_dim_changed_(false)
    {
        
    }
    
    bool is_face() {
        
        
        return bitmap_(x_, y_, z_) && !bitmap_(x_ + norm_.x, y_ + norm_.y, z_ + norm_.z );
        
//         switch( dir_ ) {
//         case plane::dir_zx_p:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_, y_ + 1, z_);
//         case plane::dir_zx_n:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_, y_ - 1, z_);
//         case plane::dir_yz_p:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_ + 1, y_, z_);
//         case plane::dir_yz_n:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_ - 1, y_, z_);
//         case plane::dir_xy_p:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_, y_, z_ + 1);
//         case plane::dir_xy_p:
//             return bitmap_(x_, y_, z_) && !bitmap_(x_, y_, z_ + 1);
//             
//         }
    }
    
    bool inc() {
        switch( dir_ ) {
        case plane::dir_zx_p:
        case plane::dir_zx_n:
            ++x_;
            if( x_ >= bitmap_.x() ) {
                ++z_;
                x_ = 0;
            }
            
            last_dim_changed_ = z_ >= bitmap_.z();
            if( last_dim_changed_ ) {
                ++y_;
                z_ = 0;
                last_dim_changed_ = true;
            } 
            
            return y_ >= bitmap_.y();
            
            
        case plane::dir_yz_p:
        case plane::dir_yz_n:
            ++z_;
            if( z_ >= bitmap_.z() ) {
                ++y_;
                z_ = 0;
            }
            
            last_dim_changed_ = y_ >= bitmap_.y();
            if( last_dim_changed_ ) {
                ++x_;
                y_ = 0;
            }
            
            return x_ >= bitmap_.x();
            
        case plane::dir_xy_p:
        case plane::dir_xy_n:
            ++y_;
            if( y_ >= bitmap_.y() ) {
                ++x_;
                y_ = 0;
            }
            last_dim_changed_ = x_ >= bitmap_.x();
            if( last_dim_changed_ ) {
                ++z_;
                x_ = 0;
            }
            
            return z_ >= bitmap_.z();
        default:
            throw std::runtime_error( "unknown plane::dir_type" );    
        }
        
    }
    vec3i pos() const {
        return vec3i( x_, y_, z_ );
    }
    
    bool last_dim_changed() const {
        return last_dim_changed_;
    }
    vec2i pos2d() const {
        switch( dir_ ) {
        case plane::dir_zx_p:
        case plane::dir_zx_n:
            return vec2i(z_,x_);
        case plane::dir_yz_p:
        case plane::dir_yz_n:
            return vec2i(y_,z_);
        case plane::dir_xy_p:
        case plane::dir_xy_n:
        default: 
            return vec2i(x_,y_);
        }
    }
    
    int last_dim_pos() const {
        switch( dir_ ) {
        case plane::dir_zx_p:
        case plane::dir_zx_n:
            return y_;
        case plane::dir_yz_p:
        case plane::dir_yz_n:
            return x_;
        case plane::dir_xy_p:
        case plane::dir_xy_n:
        default: 
            return z_;
        }
    }
    
private:
    size_t x_, y_, z_;
    
    const bitmap3d &bitmap_;
    plane::dir_type dir_;
    vec3i norm_;
    
    bool last_dim_changed_;
};



class lightmap_atlas_slice {
public:
    lightmap_atlas_slice( plane::dir_type dir, int last_dim ) 
     : 
     dir_(dir),
     last_dim_(last_dim),
     closed_(false)
    {
        
    }
    
    ~lightmap_atlas_slice() {
        if( false ) {
            std::stringstream name;
            
            name << "la_" << plane::dir_to_string(dir_) << "_" << std::setw(4) << std::setfill('0') << last_dim_ << ".pnm";
            
            std::cout << "~lightmap_atlas() writing: " << name.str() << std::endl;
            write_pnm(name.str().c_str());
        }
    }
    
    
    
    void alloc( vec2i p, size_t plane ) {
        
        assert( !closed_ );
        assert( plane != size_t(-1) );
        alloc_tex_.push_back(p);
        alloc_tex_plane_.push_back(plane);
    }
    
    void close() {
        max_size();
        closed_ = true;
    }
    
    vec2i max_size () {
        if( !closed_ ) {
            
            vec2i smax(0,0);
        
            std::for_each( alloc_tex_.begin(), alloc_tex_.end(), [&]( vec2i p ) {
                smax.x = std::max( smax.x, p.x );
                smax.y = std::max( smax.y, p.y );
            });
            
            max_size_ = smax;

        }
        
        
        return max_size_;
    }
    
    int width() const {
        assert( closed_ );
        
        return max_size_.x;
    }
    
    int height() const {
        assert( closed_ );
        
        return max_size_.y;
    }
    
    size_t get_plane( vec2i v ) const {
        // yes I know, this makes me look like a blithering idiot...
        auto it = std::find( alloc_tex_.begin(), alloc_tex_.end(), v );
        
        if( it != alloc_tex_.end() ) {
            size_t idx = std::distance( alloc_tex_.begin(), it );
            assert( idx < alloc_tex_plane_.size() );
            return alloc_tex_plane_[idx];
        } else {
            return size_t(-1);
        }
            
    }
    
    bool is_alloc( vec2i v ) const {
        
        //return std::find( alloc_tex_.begin(), alloc_tex_.end(), v ) != alloc_tex_.end();
        
        return get_plane( v ) != size_t(-1);
    }
    
    size_t num_alloc() const {
        return alloc_tex_.size();
    }
    
    const vec2i &texel_at( size_t idx ) const {
        assert( idx < alloc_tex_.size() );
        return alloc_tex_[idx];
    }
    
    size_t texel_plane_at( size_t idx ) const {
        assert( idx < alloc_tex_plane_.size() );
        return alloc_tex_plane_[idx];
    }
private:
    
    
    std::vector<vec2i> alloc_tex_;
    std::vector<size_t> alloc_tex_plane_;
    plane::dir_type dir_;
    int last_dim_;
    mutable bool closed_;
    mutable vec2i max_size_;
    
    void write_pnm(const char* c_str) {
        vec2i smax = max_size();
        
        if( smax.x <= 0 && smax.y <= 0 ) {
            return;
            
        }
        std::ofstream os(c_str);
        
        os << "P1\n";
        os << smax.x + 1 << " " << smax.y + 1<< "\n";
        for( int i = 0; i <= smax.y; ++i ) {
            for( int j = 0; j <= smax.x; ++j ) {
                if( std::find( alloc_tex_.begin(), alloc_tex_.end(), vec2i( i, j )) != alloc_tex_.end()) {
                    os << "1 ";
                } else {
                    os << "0 ";
                }
            }
            os << "\n";
        }
    }
    
};


class binpacker {

public:
    binpacker( vec2i size ) : size_(size) {}
    
    
    typedef std::tuple<size_t,int,int,const lightmap_atlas_slice *> bin_mapping;
    std::vector<bin_mapping> realize() {
        std::vector<bin_mapping> out;
        
        for( size_t i_bin = 0; i_bin < bins_.size(); ++i_bin ) {
            const auto &bin = bins_.at(i_bin);
            
            for( const auto & level : bin.levels() ) {
                int y_ptr = level.first;
                
                for( const auto & spair : level.second.slices() ) {
                    int x_ptr = spair.first;
                    
//                     std::cout << i_bin << " " << y_ptr << " " << x_ptr << ": " << spair.second << "\n";
                    
                    out.emplace_back(i_bin, y_ptr, x_ptr, spair.second);
                }
            }
        }
        
        return out;
    }
    ~binpacker() {
        std::cout << "binpacker(): " << bins_.size() << "\n";
    }
    
    void insert( lightmap_atlas_slice *slice ) {
        auto size = slice->max_size();
        
        
        // first test if existing bin/level can take the slice
        for( auto & b : bins_ ) {
            auto l = b.fit(size.y, size.x);
            
            if( l != nullptr ) {
                l->insert( size.x, slice );
                return;
            }
        }
        
        // secondly create new level in first bin with enough space
        for( auto & b : bins_ ) {
            auto l = b.create_level( size.y );
            
            if( l != nullptr ) {
                assert( l->fit(size.x ) && "bin width to small");
                l->insert(size.x, slice );
                
                return;
            }
        }
        
        // thirdly create new bin and level
        bins_.emplace_back( size_ );
        auto l = bins_.back().create_level( size.y );
        
        assert( l != nullptr && "could not create level in new bin. impossible to fit" );
        assert( l->fit(size.x ) && "bin width to small in new bin. impossible to fit");
        l->insert(size.x, slice );
    }
    
    
private:
    class bin {
    public:
        
        
        bin( bin && ) = default;
        bin & operator=(bin &&) = default;
        
        bin( vec2i size ) : width_(size.x), height_(size.y), ptr_(0) {}
    
        
        class level {
        public:
            
            typedef std::pair<int,const lightmap_atlas_slice *> ps_pair;
            
            level( int height, int width ) : height_(height), width_(width), ptr_(0) {
                
            }
            
            bool fit( int width ) {
                return ptr_ + width <= width_;
            }
            
            void insert( int width, const lightmap_atlas_slice *slice ) {
                assert( ptr_ <= width_ );
                slices_.emplace_back( ptr_, slice );
                ptr_ += width;
                
                assert( ptr_ <= width_ );
            }
            
            int height() const {
                return height_;
            }
            
            const std::vector<ps_pair> &slices() const {
                return slices_;
            }
            
        private:
            int height_;
            int width_;
            int ptr_;
            
            
            std::vector<ps_pair> slices_;
        };
        
        typedef std::pair<int,level> pl_pair;
        
        level *fit( int height, int width ) {
            for( auto & l : levels_ ) {
                if( l.second.height() >= height ) {
                    if( l.second.fit( width ) ) {
                        return &(l.second);
                    }
                }
             
            }
            return nullptr;
        }
        
        level *create_level( int height ) {
            if( ptr_ + height <= height_ ) {
                levels_.emplace_back(ptr_, level(height, width_));
                ptr_ += height;
                return &(levels_.back().second);
            }
            return nullptr;
        }
        const std::vector<pl_pair> &levels() const {
            return levels_;
        }
      
        
    private:
        
        int width_;
        int height_;
        std::vector<pl_pair> levels_;
        int ptr_;
        
    };
    
    const vec2i size_;
    std::vector<bin> bins_;
   
    
    
};


void scene_static::init_solid(const std::vector< crystal_bits::matrix_ptr >& slices) {
    const auto &slice0 = *slices.at(0);
    size_t size_z = slice0.size1();
    size_t size_x = slice0.size2();

    size_t size_y = *std::max_element( slice0.data().begin(), slice0.data().end() ) + 1;

    std::cout << "size: " << size_x << " " << size_z << " " << size_y << "\n";

    solid_ = bitmap3d( size_x, size_y, size_z );

    bool add = true;
    for( auto &slice_ptr : slices ) {
        const auto &slice = *slice_ptr;
        
        
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
void scene_static::init_planes() {

   // vec3f base_pos( -(solid_.x() / 2.0 + 0.5), -10.5, -(solid_.z() / 2.0 + 0.5));
    //base_pos_ = base_pos;

//         std::cout << "base pos: " << base_pos_ << " " << solid_.x() << "\n";

    const auto &solidc = solid_;

    vec3i light_pos( 10, 10, 10 );

    int pump_factor_ = 2;
    float scale = 1.0 / pump_factor_;

    
#if 0
    for ( int y = 0; y < int(solid_.y()); ++y ) {
        for ( size_t z = 0; z < solid_.z(); ++z ) {
            for ( size_t x = 0; x < solid_.x(); ++x ) {
                if ( solid_(x, y, z) ) {

                    const bool occ = util::occluded( light_pos, vec3i(x,y,z), solid_ );

                    float energy = occ ? 0.2 : 1.0;



                    if ( !solidc(x,y,z+1)) {
                        planes_.push_back( plane( plane::dir_xy_p, base_pos_, vec3i( x, y, z ), scale, energy));
                    }
                    if ( !solidc(x,y,z-1)) {
                        planes_.push_back( plane( plane::dir_xy_n, base_pos_, vec3i( x, y, z ), scale, energy));
                    }
                    if ( !solidc(x+1,y,z)) {
                        planes_.push_back( plane( plane::dir_yz_p, base_pos_, vec3i( x, y, z ), scale, energy));
                    }
                    if ( !solidc(x-1,y,z)) {
                        planes_.push_back( plane( plane::dir_yz_n, base_pos_, vec3i( x, y, z ), scale, energy));
                    }
                    if ( !solidc(x,y+1,z)) {
                        planes_.push_back( plane( plane::dir_zx_p, base_pos_, vec3i( x, y, z ), scale, energy));
                    }
                    if ( y > 0 && !solidc(x,y-1,z)) {
                        planes_.push_back( plane( plane::dir_zx_n, base_pos_, vec3i( x, y, z ), scale, energy));
                    }


                }
            }
        }
    }
#else
    for( auto dir : {plane::dir_zx_p, plane::dir_zx_n, plane::dir_yz_p, plane::dir_yz_n, plane::dir_xy_p, plane::dir_xy_n} ) {
        face_iterator it(solidc, dir);
        
        do {
            if( it.is_face() ) {
                vec3i pos = it.pos();
                
                if( dir == plane::dir_zx_n && pos.y == 0 ) {
                    continue; // skip faces on the underside of the level
                }
                planes_.push_back( plane( dir, base_pos_, pos, scale, 1.0));
            }
            
        } while( !it.inc() );
        
    }
#endif
//      for( size_t z = 0; z < height_.size(); ++z ) {
//          for( size_t x = 0; x < height_[z].size(); ++x ) {
//              planes_.push_back( plane( plane::dir_zx_n, base_pos, vec3i( x, 20, z ), 0.0));
//          }
//      }


    std::cout << "planes: " << planes_.size() << "\n";
    planes_.shrink_to_fit();



}
const static std::array<vec3f, 4>reorder_strip( const std::array<vec3f, 4> &in, plane::dir_type dir ) {
    switch( dir ) {
    case plane::dir_zx_p:
    case plane::dir_yz_p:
    case plane::dir_xy_p:
    default:
        return std::array<vec3f, 4>{{in[0], in[1], in[3], in[2]}};
        
        
        
    case plane::dir_zx_n:
    case plane::dir_yz_n:
    case plane::dir_xy_n:
        //return std::array<vec3f, 4>{in[2], in[1], in[3], in[0]};
        return std::array<vec3f, 4>{{in[3], in[0], in[2], in[1]}};
    };
    
}

static void write_rgb_pnm( std::ostream &os, int width, int height, const std::vector<vec3i> &bm, std::function<size_t (int,int)> coord ) {
    os << "P3\n";
    os << width << " " << height << "\n";
    os << 255 << "\n";
    
    for( int i = 0; i < height; ++i ) {
        for( int j = 0; j < width; ++j ) {
            vec3i col = bm[coord(j, i)];
            
            os << col.r << " " << col.g << " " << col.b << " ";
            
        }
        
        os << "\n";
    }
    
}

// static void fill_rect( int x, int y, int width, int height, std::vector<vec3i> *bm, std::function<size_t (int,int)> coord, vec3i color ) {
//     for( int i = y; i < y + height; ++i ) {
//         for( int j = x; j < x + width; ++j ) {
//             (*bm)[coord(j, i)] = color;
//         }
//     }
// }

static void blit_slice( int x, int y, const lightmap_atlas_slice *slice, std::vector<vec3i> *bm, std::function<size_t (int,int)> coord, vec3i color ) {
    for( int i = y; i < y + slice->height(); ++i ) {
        for( int j = x; j < x + slice->width(); ++j ) {
            
            if( !true || slice->is_alloc(vec2i(j - x, i - y)) ) {
                (*bm)[coord(j, i)] = color;
            }
        }
    }
}


static void visualize( int width, int height, const std::vector< binpacker::bin_mapping > &mapping) {
    
    auto coord = [=]( int x, int y ) { assert( x >= 0 && x < width ); assert( y >= 0 && y < height ); return x + y * width; };
    std::vector<std::vector<vec3i>> binmaps;
    
    vec3i next_col( 128, 128, 128 );
    
    for( const auto &bm : mapping ) {
        size_t bin;
        int y;
        int x;
        const lightmap_atlas_slice *slice;

        std::tie(bin, y, x, slice) = bm; // c++11 is just ridiculously good...
        
        std::cout << bin << " " << y << " " << x << ": " << slice << "\n";
        
        while( bin >= binmaps.size() ) {
            binmaps.emplace_back( width * height, vec3i(0,0,0));
        }
        
        
        //fill_rect( x, y, slice->width(), slice->height(), &binmaps.at(bin), coord, next_col );
        blit_slice( x, y, slice, &binmaps.at(bin), coord, next_col );
        
        // get some variation into the colors
        next_col.r = (next_col.r + 37) % 256;
        next_col.g = (next_col.g + 53) % 256;
        next_col.b = (next_col.b + 79) % 256;
    }
    
    for( size_t i = 0; i < binmaps.size(); ++i ) {
        std::stringstream ss;
        ss << "bin_" << std::setw(4) << std::setfill('0') << i << ".pnm";
        
        std::ofstream os( ss.str().c_str() );
        assert( os.good() );
        write_rgb_pnm(os, width, height, binmaps[i], coord);
        
    }
    
}

// TODO: review: make this smaller e.g. char*small*small should be enough
typedef std::tuple<size_t,int,int> texel_address;

std::vector<texel_address> plane_to_texel_map( const std::vector< binpacker::bin_mapping > &mapping ) {
    std::vector<texel_address> out;
    
    for( const auto &bm : mapping ) {
        size_t bin;
        int y;
        int x;
        const lightmap_atlas_slice *slice;

        std::tie(bin, y, x, slice) = bm;
        
        const size_t num_alloc = slice->num_alloc();
        
        for( size_t i = 0; i < num_alloc; ++i ) {
            vec2i pos = slice->texel_at(i);
            size_t plane = slice->texel_plane_at(i);
            
            
//             std::cout << "plane: " << plane << "\n";
            if( out.size() <= plane ) {
                out.resize( plane + 1, std::make_tuple<size_t,int,int>(-1,-1,-1));
            }
            
            out.at(plane) = std::make_tuple(bin,x + pos.x, y + pos.y );
        }
    }
    
    std::for_each( out.begin(), out.end(), [] (const texel_address &t ) {
        if( std::get<0>(t) == size_t(-1) ) {
            std::cout << "meeeeep: bad texel\n";
        }
    } );
    return out;
}

void scene_static::init_strips() {

   
    const auto &solidc = solid_;

    vec3i light_pos( 10, 10, 10 );

    int pump_factor_ = 2;
    float scale = 1.0 / pump_factor_;

    std::vector<lightmap_atlas_slice> atlases;
    
    
    bool restart = false;
    for( plane::dir_type dir : {plane::dir_zx_p, plane::dir_zx_n, plane::dir_yz_p, plane::dir_yz_n, plane::dir_xy_p, plane::dir_xy_n} ) {
//     for( auto dir : {plane::dir_zx_n, plane::dir_zx_n} ) {
        face_iterator it(solidc, dir);
        atlases.push_back(lightmap_atlas_slice( dir, it.last_dim_pos()));
        do {
            if( it.last_dim_changed() ) {
                atlases.push_back(lightmap_atlas_slice(dir, it.last_dim_pos()));
            }
            
            auto &atlas = atlases.back();
            
            if( it.is_face() ) {
                vec3i pos = it.pos();
                
                if( dir == plane::dir_zx_n && pos.y == 0 ) {
                    continue; // skip faces on the underside of the level
                }
                
                
                atlas.alloc(it.pos2d(), planes_.size());
                
                
                planes_.push_back( plane( dir, base_pos_, pos, scale, 1.0));
                
                auto &qverts = planes_.back().verts();
                        
                auto verts = reorder_strip(qverts, dir);
                uint32_t first_idx = strip_vecs_.size();
                if( restart ) {
                    //strip_idx_.push_back(restart_idx);
                    // output degenerated tris to jump to restart position
                    if( true ) {
                        strip_idx_.push_back(strip_vecs_.size()-1);
                        strip_idx_.push_back(strip_vecs_.size());
                    } else {
                        strip_idx_.push_back(strip_vecs_.size());
                        
                        vec3f t = strip_vecs().back();
                        strip_vecs_.push_back(t);
                        strip_idx_.push_back(strip_vecs_.size());
                        strip_vecs_.push_back(verts[0]);
                    }
                    
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[0]);
                    
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[1]);
                    
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[2]);
                    
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[3]);
                    
                    restart = false;
                } else {
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[2]);
                            
                    strip_idx_.push_back(strip_vecs_.size());
                    strip_vecs_.push_back(verts[3]);
                }
                strip_idx_pairs_.emplace_back( first_idx, strip_vecs_.size());
                
            } else {
                restart = true;
            }
            
        } while( !it.inc() );
        restart = true;
        
    }

    
    std::cout << "planes (striped): " << planes_.size() << "\n";
    std::cout << "vecs: " << strip_vecs_.size() << "\n";
    planes_.shrink_to_fit();


    // close all altlas strips, so that the remaining stuff can be done with const refs/ptrs.
    std::for_each( atlases.begin(), atlases.end(), [](lightmap_atlas_slice &s) {s.close();});
    
    const int bin_width = 256;
    const int bin_height = bin_width;
    
    binpacker bp( vec2i( bin_width, bin_height ) );
    std::sort( atlases.begin(), atlases.end(), [] ( const lightmap_atlas_slice &s1, const lightmap_atlas_slice &s2 ) {
        return s1.height() > s2.height();
    });
    
    for( auto & s : atlases ) {
//         std::cout << "height: " << s.height() << "\n";
        
        bp.insert( &s );
    }

   auto mapping = bp.realize();
   visualize( bin_width, bin_height, mapping );
   
   auto ptm = plane_to_texel_map( mapping );
}


void light_utils::render_light(std::vector< vec3f >* emitptr, const scene_static& scene, const vec3f& light_pos, const vec3f& light_color) {
    assert( emitptr != nullptr );

    // std::cerr << "size: " << scene.planes().size() << " " << emitptr->size() << std::endl;
    assert( scene.planes().size() == emitptr->size() );

    auto &emit_rgb_ = *emitptr; // convenience

    //std::fill( emit_rgb_.begin(), emit_rgb_.end(), vec3f(0.2, 0.2, 0.2 ));
    std::fill( emit_rgb_.begin(), emit_rgb_.end(), vec3f(0.0, 0.0, 0.0 ));
    for ( size_t i = 0; i < scene.planes().size(); ++i ) {

        auto &p = scene.planes()[i];

        vec3f trace_pos = p.pos() + p.norm();

        const bool occ = true && util::occluded( light_pos, trace_pos, scene.solid() );


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
light_static::light_static(std::istream& is, uint64_t hash) {


    uint64_t hashf = -1;
    is.read( (char *) &hashf, sizeof( size_t ));

    if( !is.good() ) {
        throw std::runtime_error( "file hande bad (empty file)" );
    }

    if( hashf != hash ) {
        throw std::runtime_error( "different hash" );
    }

    size_t size1;
    is.read( (char *) &size1, sizeof( size_t ));

    if( !is.good() ) {
        throw std::runtime_error( "file hande bad (preliminary end of file)" );
    }

    f_fact_.resize(size1);
    f_target_.resize(size1);


    for ( size_t i = 0; i < size1; ++i ) {
        size_t size2;

        is.read( (char *) &size2, sizeof( size_t ));

        if( !is.good() ) {
            throw std::runtime_error( "file hande bad (preliminary end of file)" );
        }

        f_fact_[i].resize(size2);
        f_target_[i].resize(size2);

        is.read( (char*) f_fact_[i].data(), size2 * sizeof(float));
        is.read( (char*) f_target_[i].data(), size2 * sizeof(int));

        if( !is.good() ) {
            throw std::runtime_error( "file hande bad (preliminary end of file)" );
        }

        std::for_each(f_target_[i].begin(), f_target_[i].end(), [&](int t) {
            if ( size_t(t) >= size1 ) {
                throw std::runtime_error( "cannot read form factors: target out of range");
            }
        });
    }


//         f_pairs_ = init_pairs(f_fact_, f_target_);
}
void light_static::write(std::ostream& os, uint64_t hash) {
    os.write((char*) &hash, sizeof(uint64_t));

    size_t size1 = f_fact_.size();

    os.write((char*) &size1, sizeof(size_t));
    for ( size_t i = 0; i < size1; ++i ) {
        size_t size2 = f_fact_[i].size();
        os.write((char*) &size2, sizeof(size_t));

        os.write( (char*) f_fact_[i].data(), size2 * sizeof(float));
        os.write( (char*) f_target_[i].data(), size2 * sizeof(int));
    }

}


light_static setup_formfactors( const std::vector<plane> &planes_, const bitmap3d &solid_ ) {
    
    //      std::ofstream os( "ff.txt" );
    
    std::ofstream os;//( "matrix.pnm");
    os << "P1\n";
    os << planes_.size() << " " << planes_.size() << "\n";
    
    std::vector<std::vector<float> > ff2s_(planes_.size());
    std::vector<std::vector<int> > ff2_target_(planes_.size());
    
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
            
            
            if ( light_utils::normal_cull( planes_[i], planes_[j] )) {
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
    
    return light_static( std::move( ff2s_ ), std::move( ff2_target_ ));
}
void scene_static::init_solid_from_crystal(std::istream& is, size_t pump) {

    init_solid(crystal_bits::load_crystal(is, pump));
}

