#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "misc_utils.h"
#include "scene_bits.h"

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
void scene_static::init_solid(const std::vector< boost::numeric::ublas::matrix< int > >& slices) {
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
void scene_static::init_planes() {

    vec3f base_pos( -(solid_.x() / 2.0 + 0.5), -10.5, -(solid_.z() / 2.0 + 0.5));
    base_pos_ = base_pos;

//         std::cout << "base pos: " << base_pos_ << " " << solid_.x() << "\n";

    const auto &solidc = solid_;

    vec3i light_pos( 10, 10, 10 );

    int pump_factor_ = 2;
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
    planes_.shrink_to_fit();



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