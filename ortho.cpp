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



#ifdef TEST_OPENCL
#define __CL_ENABLE_EXCEPTIONS
#include "cl.hpp"
#endif

#include "cycle.h"
#include "cl_error_codes.h"
#include <execinfo.h>
#include <fstream>
#include <array>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <sstream>

// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/matrix_proxy.hpp>
// #include <boost/numeric/ublas/fwd.hpp>

#include "aligned_buffer.h"
#include "vec_unit.h"
#include "misc_utils.h"
#include "scene_bits.h"
#include "player_bits.h"
#include "crystal_bits.h"
#include "gl_bits.h"
#include "rad_core.h"

namespace ublas = boost::numeric::ublas;







// class rad_core_sse : public rad_core {
// public:
//     rad_core_sse( const std::vector<plane> &planes, const std::vector<std::vector<float> > &ffs, const std::vector<std::vector<int> > &ff_target )
//             : emit_( planes.size() ), rad_( planes.size() ), rad2_( planes.size() ),
//             ffs_(ffs), ff_target_(ff_target),
//             planes_(planes)
//     {
// 
// 
//     }
// 
//     virtual void set_emit( const std::vector<vec3f> &emit ) {
//         assert( emit_.size() == emit.size() );
// 
//         emit_.assign( std::begin(emit), std::end(emit));
//     }
// 
//     virtual bool update() {
//         do_radiosity_sse();
// 
//         return true;
//     }
// 
//     virtual void copy( std::vector<vec3f> *out ) {
//         for ( size_t i = 0; i < rad_.size(); ++i ) {
//             (*out)[i].r = rad_[i].r;
//             (*out)[i].g = rad_[i].g;
//             (*out)[i].b = rad_[i].b;
//         }
//     }
// 
// 
// private:
// 
//     void do_radiosity_sse() {
//         const int steps = 1;
//         const float min_ff = 0;
// //      std::fill(e_rad_sse_.begin(), e_rad_sse_.end(), vec3f(0.0, 0.0, 0.0));
// //      std::fill(e_rad2_sse_.begin(), e_rad2_sse_.end(), vec3f(0.0, 0.0, 0.0));
// 
//         //std::copy( emit_sse_.begin(), emit_sse_.end(), e_rad_sse_.begin() );
// 
//         typedef vector_unit<float,4> vu;
//         typedef vu::vec_t vec_t;
//         //steps = 0;
// 
//         vec_t reflex_factor = vu::set1(1.0);
//         for ( int i = 0; i < steps; ++i ) {
//             //for( auto it = pairs_.begin(); it != pairs_.end(); ++it, ++ff_it ) {
// 
// 
//             for ( size_t j = 0; j < ffs_.size(); ++j ) {
// 
//                 const size_t s = ffs_[j].size();
//                 vec_t rad = vu::set1(0);
// 
//                 const vec3f cd = planes_[j].col_diff();
//                 const vec_t col_diff = vu::set( 0, cd.b, cd.g, cd.r );
// 
// 
// 
//                 for ( size_t k = 0; k < s; ++k ) {
//                     size_t target = ff_target_[j][k];
//                     //rad_rgb += (col_diff * e_rad_rgb_[target]) * ff2s_[j][k];
// 
//                     if ( false && ffs_[j][k] < min_ff ) {
//                         continue;
//                     }
// 
//                     const vec_t ff = vu::set1( ffs_[j][k] );
// 
//                     rad = vu::add( rad, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target))), ff ));
// 
//                 }
// 
//                 vu::store( vu::add( vu::load((float*)emit_(j)), vu::mul(rad, reflex_factor)), (float*)rad2_(j));
// //              std::cout << "col: " << rad2_[j].r << " " << cd.r <<  "\n";
//                 //e_rad2_rgb_[j] = emit_rgb_[j] + rad_rgb;// * reflex_factor;
// 
//             }
// 
// 
// //            rad_.swap(rad2_);
//             rad_ = rad2_;
// 
// 
//         }
// 
//         //e_rad_rgb_.assign( e_rad_sse_.begin(), e_rad_sse_.end() );
// 
// 
// 
//     }
// 
// 
//     aligned_buffer<col3f_sse> emit_;
//     aligned_buffer<col3f_sse> rad_;
//     aligned_buffer<col3f_sse> rad2_;
//     const std::vector<std::vector<float> > &ffs_;
//     const std::vector<std::vector<int> > &ff_target_;
//     const std::vector<plane> &planes_;
// 
// };



template<typename IDX>
std::vector<IDX> ff_size_sort_permutation( const std::vector<std::vector<float>> &ffs ) 
{
        typedef std::pair<size_t,size_t> ss_pair;
        std::vector<ss_pair> fs;
        for( size_t i = 0; i < ffs.size(); ++i ) {
            //std::cout << "ffs: " << i << " " << ff2s_[i].size() << "\n";
            fs.emplace_back( ffs[i].size(), i );
        }
        
        std::sort( fs.begin(), fs.end(), []( const ss_pair &p1, const ss_pair &p2 ) { return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first; } );
        
        std::vector<IDX> perm_size(fs.size());
        perm_size.resize(fs.size());
        
        std::transform( fs.begin(), fs.end(), perm_size.begin(), [](const ss_pair &p) { return p.second; } );
      
        return perm_size;
}

class light_dynamic {
public:
    light_dynamic() {}
    light_dynamic( size_t num ) : emit_( num ), rad_( num ) {}
    
    void clear_emit() {
        std::fill( emit_.begin(), emit_.end(), vec3f(0.0, 0.0, 0.0));
    }

    std::vector<vec3f> *emit() {
        return &emit_;
    }
    
    std::vector<vec3f> *rad() {
        return &rad_;
    }
    
private:
    std::vector<vec3f> emit_;
    std::vector<vec3f> rad_;
    
    
};

std::string hash_to_filename( uint64_t hash ) {
    std::stringstream ss;
    ss << "baked";
    for( size_t i = 0; i < 8; ++i ) {
        size_t c = hash & 0xff;
        
        
        ss << std::hex << c;
        
        
        
        hash >>= 1;
        
        
    }
    ss << ".bin";
    return ss.str();
}

// class scene_unit {
// public:
//     scene_unit( std::istream &is, const vec3i &base_pos )
//       : 
//       base_pos_(base_pos),
//       scene_static_( base_pos )
//     {
//         //std::ifstream is( "cryistal-castle-hidden-ramp.txt" );
//      //    std::ifstream is( "house1.txt" );
//         //std::ifstream is( "cryistal-castle-tree-wave.txt" );
// 
//         
//     }
//     
//     light_dynamic *get_light_dynamic() {
//         return &light_dynamic_;
//     }
//     
//     void rad_update() {
//         
//     }
//     
//     size_t num_planes() {
//         return scene_static_.planes().size();
//     }
//     
//     const scene_static &get_scene_static() const {
//         return scene_static_;
//     }
// private:
//     
//     vec3i base_pos_;
//     
//     
//     
//     
// };


class render_unit {
public:
    render_unit( std::istream &is, const vec3i &base_pos )
    : base_pos_(base_pos),
    scene_static_(base_pos)
    {
        assert( is.good() );
//         height_fields_ = crystal_bits::load_crystal(is, pump_factor_);
//         std::cout << "hf: " << height_fields_.size() << "\n";
//         
//         
//         
//         scene_static_.init_solid(height_fields_);
//         
        
        const size_t pump_factor = 4;
        
         base_pos_.x *= pump_factor;
         base_pos_.z *= pump_factor;
        scene_static_.init_solid_from_crystal(is, pump_factor);
        

//        scene_static_.init_planes();
        scene_static_.init_binmaps();
        scene_static_.init_strips();
        uint64_t scene_hash = scene_static_.hash();
        auto bin_name = hash_to_filename(scene_hash);
        
        std::cout << "baked name: " << bin_name << "\n";
        try {
            std::ifstream is( bin_name.c_str() );
            
            
            light_static_ = light_static( is, scene_hash );
        } catch( std::runtime_error x ) {
            
            std::cerr << "load failed. recreating. error:\n" << x.what() << std::endl;
            
            light_static_ = setup_formfactors(scene_static_.planes(), scene_static_.solid());    
        }
        
        if( !false ) {
            std::ofstream os( bin_name.c_str() );
            light_static_.write(os, scene_hash);
        }
        
        
        light_dynamic_ = light_dynamic(scene_static_.planes().size() );
        rad_core_ = make_rad_core_threaded(scene_static_, light_static_);
        
        
        
//         vbob_ = vbo_builder(scene_static_.planes().size() );
//         vbob_.update_index_buffer( scene_static_.planes().size());
//         vbob_.update_vertices( scene_static_.planes().begin(), scene_static_.planes().end());
//         
        vbob_ts_ = vbo_builder_tristrip( scene_static_.tristrip_at(0) );
        
    }
    
    void clear_emit() {
        light_dynamic_.clear_emit();
    }
    void render_light( const vec3f &pos, const vec3f &color ) {
        light_utils::render_light(light_dynamic_.emit(), scene_static_, pos - base_pos_, color);
    }
    
    void update() {
        
        rad_core_->set_emit( *light_dynamic_.emit() );
        rad_core_->copy( light_dynamic_.rad() );
        
        
//         vbob_.update_color( light_dynamic_.rad()->begin(), light_dynamic_.rad()->end());
        vbob_ts_.update_color( light_dynamic_.rad()->begin(), light_dynamic_.rad()->end());
        
    }
    void draw() {
//         vbob_.draw_arrays();
        vbob_ts_.draw_arrays();
    }
private:
    vec3i base_pos_;
    
    scene_static scene_static_;
    
    light_static light_static_;
    light_dynamic light_dynamic_;
    
    std::unique_ptr<rad_core> rad_core_;
    
    
    
//     vbo_builder vbob_;
    vbo_builder_tristrip vbob_ts_;
};


#ifdef TEST_OPENCL
class rad_core_opencl : public rad_core {
public:
    
    rad_core_opencl( cl::Context ctx, cl::CommandQueue queue, const scene_static &scene_static, const light_static &light_static ) 
    : ctx_(ctx), queue_(queue), num_planes_(scene_static.planes().size()), buf_dir_(true),
    pints_(0), pints_last_(0), pints_last_time_(0)
    {
        
        size_t num_facts = 0;
        std::vector<cl_int> stage_offsets;
        std::vector<cl_int> stage_num;
        
        std::vector<cl_float> stage_fact;
        std::vector<cl_int> stage_target;

        
        auto & fact = light_static.f_fact();
        auto & target = light_static.f_target();
        
        for( size_t i = 0; i < fact.size(); ++i ) {
        //std::for_each(light_static.f_fact().begin(), light_static.f_fact().end(), [&]( const std::vector<float> & f) {
            stage_offsets.push_back(num_facts);
            stage_num.push_back(fact[i].size());
            num_facts += fact[i].size();
            
            stage_fact.insert( stage_fact.end(), fact[i].begin(), fact[i].end() );
            stage_target.insert( stage_target.end(), target[i].begin(), target[i].end() );
        }
        
        assert( num_facts < size_t(std::numeric_limits<cl_int>::max()) );
        assert( stage_fact.size() == num_facts );
        assert( stage_target.size() == num_facts );
        
        std::vector<cl_int> stage_perm = ff_size_sort_permutation<cl_int>( light_static.f_fact());
        assert( stage_perm.size() == scene_static.planes().size());
        pints_iter_ = num_facts;
        
        //pints_iter_ = *std::max_element(stage_num.begin(), stage_num.end()) * stage_num.size();
        
        std::vector<cl_float4> stage_col_diff;//( num_planes_ );
        for( auto & p : scene_static.planes() ) {
            vec3f ci = p.col_diff();
            
            cl_float4 co;
            co.s[0] = ci.r;
            co.s[1] = ci.g;
            co.s[2] = ci.b;
            co.s[3] = 1.0;
            stage_col_diff.push_back( co );
        }
        
        
        f_fact_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, num_facts * sizeof(cl_float) );
        f_target_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, num_facts * sizeof(cl_int) );
        f_offsets_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, stage_offsets.size() * sizeof(cl_int) );
        f_num_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, stage_num.size() * sizeof(cl_int) );
        
        col_diff_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, stage_col_diff.size() * sizeof(cl_float4) );
        
        f_perm_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, stage_perm.size() * sizeof(cl_int) );
        
        queue_.enqueueWriteBuffer( f_fact_, false, 0, num_facts * sizeof(cl_float), stage_fact.data() );
        queue_.enqueueWriteBuffer( f_target_, false, 0, num_facts * sizeof(cl_int), stage_target.data() );
        queue_.enqueueWriteBuffer( f_offsets_, false, 0, stage_offsets.size() * sizeof(cl_int), stage_offsets.data() );
        queue_.enqueueWriteBuffer( f_num_, false, 0, stage_num.size() * sizeof(cl_int), stage_num.data() );
        queue_.enqueueWriteBuffer( col_diff_, false, 0, stage_col_diff.size() * sizeof(cl_float4), stage_col_diff.data() );
        queue_.enqueueWriteBuffer( f_perm_, false, 0, stage_perm.size() * sizeof(cl_int), stage_perm.data() );
        
        f_emit_ = cl::Buffer( ctx_, CL_MEM_WRITE_ONLY, light_static.num_planes() * sizeof(cl_float4) );
        f_rad_ = cl::Buffer( ctx_, CL_MEM_READ_WRITE, light_static.num_planes() * sizeof(cl_float4) );
        f_rad2_ = cl::Buffer( ctx_, CL_MEM_READ_WRITE, light_static.num_planes() * sizeof(cl_float4) );
    
        {
            cl_float4 zero{{0.0, 0.0, 0.0, 0.0}};
            std::vector<cl_float4> stage_zero(light_static.num_planes(), zero);
            
            queue_.enqueueWriteBuffer( f_rad_, false, 0, light_static.num_planes() * sizeof(cl_float4), stage_zero.data());
            queue_.enqueueWriteBuffer( f_rad2_, false, 0, light_static.num_planes() * sizeof(cl_float4), stage_zero.data());
            
        }
        
        queue_.finish();
        
    }
    
    void load_kernel( cl::Context ctx, const std::vector<cl::Device> &used_devices ) {
        std::string source;
        {
            std::ifstream is( "rad_core.cl" );
            
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
        
        cl_int err;
        try {
            program_ = cl::Program( ctx, source, false, &err );
            
            program_.build();
        } catch( cl::Error x ) {
            std::cout << "cl error during build: " << x.what() << " " << cl_str_error(x.err()) << "\n"; 
            
            if( x.err() == CL_BUILD_PROGRAM_FAILURE ) {
                auto log = program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>( used_devices.front() );
                
                std::cout << "build error. log: \n" << log << "end of log.\n";
            }
            
            throw;
        }
        assert( err == CL_SUCCESS );
        
        kernel_ = cl::Kernel( program_, "rad_kernel", &err );
        assert( err == CL_SUCCESS );   
    }
    
    virtual void set_emit( const std::vector<vec3f> &emit ) {
        stage_.resize(emit.size());
        
        std::transform( emit.begin(), emit.end(), stage_.begin(), []( vec3f in ) {
            cl_float4 t;
            t.s[0] = in.r;
            t.s[1] = in.g;
            t.s[2] = in.b;
            return t;
        });
        queue_.enqueueWriteBuffer( f_emit_, true, 0, stage_.size() * sizeof(cl_float4), stage_.data());
    }
//     virtual bool update() = 0;
    virtual void copy( std::vector<vec3f> *out ) {
        
    }
    void run() {
        buf_dir_ = !buf_dir_;
        kernel_.setArg(0, cl_int(num_planes_));
        kernel_.setArg(1, f_fact_);
        kernel_.setArg(2, f_target_);
        kernel_.setArg(3, f_offsets_);
        kernel_.setArg(4, f_num_);
        kernel_.setArg(5, col_diff_);
        kernel_.setArg(6, f_emit_);
        
        if( !buf_dir_ ) {
            kernel_.setArg(7, f_rad_);
            kernel_.setArg(8, f_rad2_);
        } else {
            kernel_.setArg(8, f_rad_);
            kernel_.setArg(7, f_rad2_);
        }
        
        kernel_.setArg( 9, f_perm_ );
        
        queue_.enqueueNDRangeKernel( kernel_, cl::NDRange(0), cl::NDRange(num_planes_) );
        queue_.finish();
        
        pints_ += pints_iter_;
        
    }
    
    cl::Buffer f_rad() {
        cl_ubyte64 time = CL_System::get_microseconds();
        
        cl_ubyte64 dt = time - pints_last_time_;
//        std::cout << dt << "\n";
        if( dt >= 1e6 ) {
            const size_t dp = pints_ - pints_last_;
            pints_last_ = pints_;
            
            pints_last_time_ = time;
            
            std::cout << "gpu pint/s: " << dp / (dt / 1.0e6) << "\n";
        }
        
        
        if( buf_dir_ ) {
            return f_rad_;
        } else {
            return f_rad2_;
        }
    }
private:
    cl::Context ctx_;
    cl::CommandQueue queue_;
    
    cl::Program program_;
    cl::Kernel kernel_;
    
    cl::Buffer f_fact_;
    cl::Buffer f_target_;
    cl::Buffer f_offsets_;
    cl::Buffer f_num_;
    
    cl::Buffer col_diff_;
    
    cl::Buffer f_emit_;
    cl::Buffer f_rad_;
    cl::Buffer f_rad2_;
    
    cl::Buffer f_perm_;
    
    const size_t num_planes_;
    std::vector<cl_float4> stage_;
    
    bool buf_dir_;
    uint64_t pints_;
    uint64_t pints_iter_;
    uint64_t pints_last_;
    uint64_t pints_last_time_;
};

#endif


class ortho {

public:

 
    ortho() :
      pump_factor_(4) 
    {

//         glVertexPointer(4, GL_FLOAT, 0, ptr);
        //glEnableClientState( GL_VERTEX_ARRAY );
//         glColorPointer(

        //std::ifstream is( "cryistal-castle-hidden-ramp.txt" );
//         std::ifstream is( "house1.txt" );
        //std::ifstream is( "cryistal-castle-tree-wave.txt" );

//         assert( is.good() );
//         height_fields_ = crystal_bits::load_crystal(is, pump_factor_);
//         std::cout << "hf: " << height_fields_.size() << "\n";
//         
//         
//         
//         scene_static_.init_solid(height_fields_);
//         
        
//         scene_static_.init_solid_from_crystal(is, pump_factor_);
        

//         scene_static_.init_planes();

//         uint64_t scene_hash = scene_static_.hash();
//         
//         try {
//             std::ifstream is( "ff.bin" );
//             
//             
//             light_static_ = light_static( is, scene_hash );
//         } catch( std::runtime_error x ) {
//             
//             std::cerr << "load failed. recreating. error:\n" << x.what() << std::endl;
//             
//             light_static_ = setup_formfactors(scene_static_.planes(), scene_static_.solid());    
//         }
//         
//         if( !false ) {
//             std::ofstream os( "ff.bin" );
//             light_static_.write(os, scene_hash);
//         }
//         
//         
//         light_dynamic_ = light_dynamic(scene_static_.planes().size() );
        
        CL_OpenGLWindowDescription desc;
        desc.set_size( CL_Size( 1024, 768 ), true );

        desc.set_depth_size(16);
        //std::cout << "depth: " << desc.get_depth_size();
        
        wnd_ = CL_DisplayWindow(desc);

        CL_GraphicContext_GL gc = wnd_.get_gc();
//      //CL_Mat4f proj = CL_Mat4f::ortho( 0, 1024, 0, 768, 100, -100 );


        gc.set_active();
#ifdef TEST_OPENCL
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
#endif
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
        
        glShadeModel(GL_FLAT);
    }

    void setup_perspective( const player &camera ) {
        glMatrixMode(GL_PROJECTION);                        //hello

        {
            CL_Mat4f proj = CL_Mat4f::perspective( 60, 1.5, 0.2, 500 );
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

//     void render_quads() {
//         assert(false && "out of order" );
//         glBegin(GL_QUADS);
//         for ( size_t i = 0; i < scene_static_.planes().size(); ++i ) {
//             const plane &p = scene_static_.planes()[i];
// 
// //          p.energy_rgb(ls.rad_rgb(i));
// 
//             p.render();
//         }
//         glEnd();
// 
//     }

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

        //int steps = 1;

        //light_scene ls( planes_, solid_ );


        auto t_old = CL_System::get_microseconds();
        bool light_on = true;
        bool light_button_down = false;
        double delta_t = 0.01;
        player p1;

        //rad_core_ = make_unique<rad_core_threaded>( scene_static_, light_static_ );
        
//         rad_core_ = make_rad_core_threaded(scene_static_, light_static_);
        
        
//         auto rad_core2 = make_unique<rad_core_opencl>( cl_context_, cl_cqueue_, scene_static_, light_static_ );
//         rad_core2->load_kernel( cl_context_, cl_used_devices_ );
        
        
//      float min_ff = 5e-5;

        
//         vab.render(planes_.begin(), planes_.end() );
        
//         vbo_builder vbob(scene_static_.planes().size());
//         vbob.update_index_buffer(scene_static_.planes().size());
//         vbob.update_vertices( scene_static_.planes().begin(), scene_static_.planes().end());
        
         std::ifstream is( "cryistal-castle-hidden-ramp.txt" );
        //std::ifstream is( "house1.txt" );
        render_unit runit(is, vec3f( -40.0, -20.0, -40.0 ));
        
//         std::ifstream is2( "cryistal-castle-tree-wave.txt" );
//         std::ifstream is2( "cryistal-castle-hidden-ramp.txt" );
//         std::ifstream is2( "house1.txt" );
//         render_unit runit2(is2, vec3f( 60.0, -20.0, 0.0 ));
        
#if 0
        cl::BufferGL buf;
        //cl::Buffer buf;
        try {
            
            //cl_int cl_err;
            //buf = cl::Buffer( clCreateFromGLBuffer( cl_context_(), CL_MEM_WRITE_ONLY, vbob.buffers_[1], &cl_err ));
            //assert( cl_err == CL_SUCCESS );
            buf = cl::BufferGL( cl_context_, CL_MEM_WRITE_ONLY, vbob.buffers_[1] );
        
        
        
            cl_fcolor_ = cl::Buffer( cl_context_, CL_MEM_READ_ONLY, scene_static_.planes().size() * 3 * sizeof(float) );
            
            cl_kernel_.setArg(0, buf() );
            //cl_kernel_.setArg(1, cl_fcolor_ );
            cl_kernel_.setArg(1, rad_core2->f_rad() );
            cl_uint cl_color_size = scene_static_.planes().size();
            cl_kernel_.setArg(2, cl_color_size );
            
            
            
        } catch( cl::Error x ) {
            std::cerr << "cl error during gl buffer setup\ncall: " << x.what() << "\nerror code: " << cl_str_error( x.err() ) << "\n";            
            throw;
        }
#endif
        bool light_changed = true;
        
        //rad_core_ = make_unique<rad_core_lockfree>( scene_static_, light_static_ );
        
        bool do_quit = false;
        
        while ( !do_quit ) {

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

            int light_speed = 1;
            
            do_quit = keyboard.get_keycode(CL_KEY_Q);
            
            if ( keyboard.get_keycode(CL_KEY_LEFT) ) {
                light_pos.x += light_speed;
                light_changed = true;
            }
            if (  keyboard.get_keycode(CL_KEY_RIGHT) ) {
                light_pos.x -= light_speed;
                light_changed = true;
            }
            if ( keyboard.get_keycode(CL_KEY_UP) ) {
                light_pos.z += light_speed;
                light_changed = true;
            }
            if (  keyboard.get_keycode(CL_KEY_DOWN) ) {
                light_pos.z -= light_speed;
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
                //ls.reset_emit();
//                 light_dynamic_.clear_emit();
                runit.clear_emit();
//                 runit2.clear_emit();
                //ls.render_light(vec3f( 40, 50.0, light_x ), vec3f(1.0, 1.0, 1.0 ));
                if ( light_on ) {
                    //ls.render_light(light_pos, vec3f(1.0, 0.8, 0.6 ));
                    
                    //vec3f light_pos( p1.pos)
                    
                    //              vec3f light_pos = (p1.pos()* pump_factor_) - base_pos_;
                    vec3f light_weird = (light_pos * pump_factor_) - base_pos_;
//                     light_utils::render_light( light_dynamic_.emit(), scene_static_, light_weird, vec3f(1.0, 0.8, 0.6 ));
                    runit.render_light(light_weird, vec3f(1.0, 0.8, 0.6 ));
//                     runit2.render_light(light_weird, vec3f(1.0, 0.8, 0.6 ));
                }
                //ls.post();
                light_changed = false;
            }
            runit.update();
//             runit2.update();
//             rad_core2->set_emit( *light_dynamic_.emit() );
//             rad_core_->set_emit( *light_dynamic_.emit() );
            //ls.render_emit_patches();

            //steps = 1;
            //ls.do_radiosity( steps );

            //rad_core2->run();
//             rad_core_->copy( light_dynamic_.rad() );
             // stupid: transfer rgb energy fomr light scene to planes
//             for ( size_t i = 0; i < scene_static_.planes().size(); ++i ) {
//                 plane &p = const_cast<plane &>(scene_static_.planes()[i]); // FIXME: HACK!!! is the energy_rgp stored in the planes actually used? remove it!
//                 p.energy_rgb((*light_dynamic_.rad())[i]);
//             }
//             
            light_x += light_xd;


//          glPushMatrix();
            
            //CL_Mat4f proj = CL_Mat4f::look_at( 20, 20, 20, 0, 0, 0, 0.0, 1.0, 0.0 );

            //vab.update_color( planes_.begin(), planes_.end() );
//             vab.update_color( ls.rad_rgb().begin(), ls.rad_rgb().end() );
            
//             vab.setup_gl_pointers();
//             vbob.update_color(light_dynamic_.rad()->begin(), light_dynamic_.rad()->end());
#if 0
            try
            {
                cl_int err;
//                 glFinish();
                const auto &rad_colors = *light_dynamic_.rad();
//                 
                cl_cqueue_.enqueueWriteBuffer( cl_fcolor_, false, 0, rad_colors.size() * sizeof(vec3f), rad_colors.data() );
//                 
                

                err = clEnqueueAcquireGLObjects( cl_cqueue_(), 1, &(buf()), 0, 0, 0 );
                assert( err == CL_SUCCESS );


                
                //cl_kernel_.setArg(1, rad_core2->f_rad() );
                cl_kernel_.setArg(1, cl_fcolor_ );
                cl_cqueue_.enqueueNDRangeKernel( cl_kernel_, 0, cl::NDRange(scene_static_.planes().size()) );
                
                               

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
//             vbob.draw_arrays();
            runit.draw();
//             runit2.draw();
            wnd_.flip(1);



            auto t = CL_System::get_microseconds();

//             std::cout << "fps: " << 1e6 / (t - t_old) << "\n";
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
        try {
            
     
            ortho o;
            o.start();
        } catch( gl_error_exception x ) {
            std::cerr << x.what() << std::endl;
            std::cerr << "bailing out\n";
        }

        return 0;
    }

private:

    
        




//     void light_planes( const vec3i &light_pos ) {
//         for ( plane &p : planes_ ) {
// 
//             vec3f trace_pos = p.pos() + p.norm();
// 
//             const bool occ = true && util::occluded( light_pos, trace_pos, solid_ );
// 
// 
//             if ( !occ ) {
//                 vec3f d = light_pos - p.pos();
//                 d.normalize();
//                 float len = d.length();
//                 d /= len;
//                 float dot = d.dot( p.norm() );
// 
//                 if ( dot > 0 ) {
//                     p.energy( dot * (10/(2*3.1415*len*len)));
//                 } else {
//                     p.energy(0.1);
//                 }
// 
//             } else {
//                 p.energy( 0.1 );
//             }
// 
// 
//         }
//     }

#ifdef TEST_OPENCL
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
#endif

    CL_SetupCore setup_core_;




    CL_SetupDisplay display_;
    //CL_SetupSWRender swrender;

    CL_SetupGL setup_gl_;
    CL_DisplayWindow wnd_;
    GLuint texName;

    
    //std::vector<crystal_bits::matrix_ptr> height_fields_;
    
//     

    const size_t pump_factor_;
    vec3f base_pos_;
#ifdef TEST_OPENCL    
    cl::Platform cl_platform_;
    std::vector<cl::Device> cl_used_devices_;
    cl::Context cl_context_;
    
    cl::Program cl_program_;
    cl::Kernel cl_kernel_;
    
    cl::CommandQueue cl_cqueue_;
    cl::Buffer cl_fcolor_;
#endif
//     scene_static scene_static_;
//     
//     light_static light_static_;
//     light_dynamic light_dynamic_;
//     
//     std::unique_ptr<rad_core> rad_core_;
};



CL_ClanApplication app(&ortho::main);
