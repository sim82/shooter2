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

#ifndef __gl_bits_h
#define __gl_bits_h

#include <cassert>
#include <stdexcept>

#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>

#include "scene_bits.h"


#define nullptr 0
class scene_static;
class gl_program;

// the only valid reason to use the preprocessor except for include guards ;-)
#define check_gl_error {GLenum x = glGetError(); if( x != GL_NO_ERROR ) {throw gl_error_exception(x, __FILE__, __LINE__);} }
class gl_error_exception : public std::runtime_error {
public:
    gl_error_exception( GLenum err, const char *file = nullptr, int line = -1 ) 
    : std::runtime_error( err_str(err, file, line) ) {}
private:
    std::string err_str( GLenum err, const char *file, int line ) throw();
    
};

class gl_utils {
public:
    template<const int minv, const int maxv >
    static inline GLubyte clamp( int v ) {
        //return std::max( minv, std::min( maxv, v ));
        return std::min( maxv, v );
        //return v;
    }
    
};


static void android_assert( const char *file, const int line, const char *func, const char *exp ) {
    std::stringstream ss;
    
    ss << "assert: " << file << ":" << line << " " << func << ":\n" << exp;
    throw std::runtime_error( ss.str() );
}

#undef assert
#  define   assert(e)   ((e) ? (void)0 : android_assert(__FILE__, __LINE__, __func__, #e))

class vbo_builder_tristrip {
public:
#if 0
    vbo_builder_tristrip( const vbo_builder_tristrip & ) = delete;
    const vbo_builder_tristrip &operator=( const vbo_builder_tristrip & ) = delete;

    
    vbo_builder_tristrip( vbo_builder_tristrip && ) = default;
    vbo_builder_tristrip &operator=( vbo_builder_tristrip && ) = default;
#endif
    vbo_builder_tristrip() : scene_static_(0) {}
    
    vbo_builder_tristrip( const scene_static &scene );
    
    void update_color_vec3fptr( const vec3f*const first, const vec3f*const last );
    
    template<typename iiter>
    void update_color( iiter first, iiter last ) {
        assert( std::distance( first, last ) == ptrdiff_t(num_planes_) );
        // this is kind of crude, but update_color asserts that [first,last]
        // point to contiguous memory, so using an non-contig container
        // should be caught.
        update_color_vec3fptr( static_cast<vec3f*>(&(*first)), static_cast<vec3f*>(&(*last)) );
        //         const size_t num_planes_ = std::distance( first, last );


    }
    
    void draw_arrays( gl_program & program ) ;
    
    
    
private:
    const scene_static *scene_static_;
    
    size_t color_size;
    
    GLuint buffers_[2];
    GLuint index_buffer_;
    
    size_t num_planes_;
    size_t index_num_;
};
#if 0
class vbo_builder {
public:
    vbo_builder( const vbo_builder & ) = delete;
    const vbo_builder &operator=( const vbo_builder & ) = delete;
    
    vbo_builder( vbo_builder && ) = default;
    vbo_builder &operator=( vbo_builder && ) = default;
    
    vbo_builder() : color_size(0), buffers_{GLuint(-1),GLuint(-1)}, num_planes_(0) {}
    
    vbo_builder( size_t num_planes ) ;
    
    void update_index_buffer( size_t num_planes ) ;
    
    template<typename iiter>
    void update_vertices( iiter first, iiter last ) {
        //const size_t num_planes_ = std::distance( first, last );
        
        glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] );
        float *b = (float*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY );
        assert( b != nullptr );
        for( ; first != last; ++first ) {
         
            const typename iiter::reference p = *first;
            //size_t s1 = vertex_.size();
            
            p.render_vertex( b );
            b += 4 * 3;
            
        }
        
        glUnmapBuffer( GL_ARRAY_BUFFER );
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
                *(b++) = gl_utils::clamp<0,255>(255 * first->r);
                *(b++) = gl_utils::clamp<0,255>(255 * first->g);
                *(b++) = gl_utils::clamp<0,255>(255 * first->b);
                *(b++) = 255;
            }   
            
        }
        
        
        
        glUnmapBuffer( GL_ARRAY_BUFFER );
    }
    
    void draw_arrays() ;
//private:
    size_t color_size;
    
    GLuint buffers_[2];
    GLuint index_buffer_;
    
    size_t num_planes_;
};


#endif
#if 1
class gl_program {
    template<typename P>
    class compare_first_string {
    public:
        bool operator()(const P &a, const P &b ) const {
            return a.first < b.first;
        }
        
        //     bool operator()(const std::string &a, const P &b ) const {
            //         return a < b.first;
            //     }
            //     
            //     bool operator()(const P &a, const std::string &b ) const {
                //         return a.first < b;
                //     }
                
        bool operator()( const char *a, const P &b ) const {
            return a < b.first;
        }
        
        bool operator()(const P &a, const char *b ) const {
            return a.first < b;
        }
        
    };
    
    
public:
    gl_program( const gl_program & ) = delete;
    gl_program &operator=( const gl_program & ) = delete;
    
    // meeeeep: default move constructor/assignment is an error (for construction-only it probably works by luck)
    gl_program( gl_program && ) = delete;
    gl_program &operator=( gl_program && ) = delete;
    
    gl_program() : program(0) {
        
    }
    
    gl_program( const char *vertex_src, const char *fragment_source ) ;
    ~gl_program() ;
    void use() ;
    void validate() ;
    GLuint mvp_handle() {
        return gv_mvp_handle;
    }
    GLuint a_position_handle() {
        return a_position_handle_;
    }
    GLuint a_color_handle() {
        return a_color_handle_;
    }
    
    
    GLuint color_handle() {
        return color_handle_;
    }
    
    GLuint uniform_handle( const char *name ) ;
    
    GLuint get_program() const {
        return program;
    }
private:
    
    
    GLuint loadShader(GLenum shaderType, const char* pSource) ;
    
    
    GLuint program;
    
    // cached attribute handles
    GLuint a_position_handle_;
    GLuint a_color_handle_;
    GLuint gv_mvp_handle;
    GLuint color_handle_;
    
    typedef std::pair<std::string,GLuint> name_handle_pair;
    
    std::vector<name_handle_pair> uniform_handles_;
};
#endif


#endif