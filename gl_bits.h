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

#include <ClanLib/gl.h>

// #include "scene_bits.h"

class scene_static;

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

class vbo_builder_tristrip {
public:
    vbo_builder_tristrip( const vbo_builder_tristrip & ) = delete;
    const vbo_builder_tristrip &operator=( const vbo_builder_tristrip & ) = delete;
    
    vbo_builder_tristrip( vbo_builder_tristrip && ) = default;
    vbo_builder_tristrip &operator=( vbo_builder_tristrip && ) = default;
    
    vbo_builder_tristrip() : tristrip_(0) {}
    
    vbo_builder_tristrip( const scene_static::tristrip &scene );
    
    void update_color_vec3fptr( const vec3f*const first, const vec3f*const last );
    
    template<typename iiter>
    void update_color( iiter first, iiter last ) {
      //  return;
#if 0
        assert( std::distance( first, last ) == ptrdiff_t(tristrip_->idx_pairs().size()) );
        // this is kind of crude, but update_color asserts that [first,last]
        // point to contiguous memory, so using an non-contig container
        // should be caught.
        update_color_vec3fptr( static_cast<vec3f*>(&(*first)), static_cast<vec3f*>(&(*last)) );
        //         const size_t num_planes_ = std::distance( first, last );
#endif

    }
    
    void update_tex_coords( const std::vector<scene_static::texel_address> &texls ) {
      //  assert( texls.size() == num_planes_ );
    }
    
    void draw_arrays() ;
    
    
    
private:
    const scene_static::tristrip *tristrip_;
    
    size_t color_size;
    
    GLuint buffers_[2];
    GLuint index_buffer_;
    
    //size_t num_planes_;
    size_t index_num_;
};

#if 0
// deprecated
class vbo_builder {
public:
    vbo_builder( const vbo_builder & ) = delete;
    const vbo_builder &operator=( const vbo_builder & ) = delete;
    
    vbo_builder( vbo_builder && ) = default;
    vbo_builder &operator=( vbo_builder && ) = default;
    
    vbo_builder() : color_size(0), buffers_({-1,-1}), num_planes_(0) {}
    
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

class gl_texture {
public:
    gl_texture( const gl_texture & ) = delete;
    gl_texture &operator=( const gl_texture & ) = delete;
    
    gl_texture( gl_texture && ) = default;
    gl_texture &operator=( gl_texture && ) = default;
    
    gl_texture() {
        glGenTextures( 1, &tex_image_ );
        
        // select our current texture
        glBindTexture( GL_TEXTURE_2D, tex_image_ );

        
        glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
        
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
        
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

        // the texture wraps over at the edges (repeat)
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
    }
    
    void upload_rgb( const std::vector<uint8_t> &buf, GLsizei width, GLsizei height ) {
        glBindTexture( GL_TEXTURE_2D, tex_image_ );
        
        glCopyTexImage2D( GL_TEXTURE_2D, 0, GL_RGB8, 0, 0, width, height, 0 );
        
    }
    ~gl_texture() {
        glDeleteTextures( 1, &tex_image_ );
    }
private:
    GLuint tex_image_;
};

#endif