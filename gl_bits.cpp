#include <stdexcept>
#include <sstream>

#include "scene_bits.h"
#include "gl_bits.h"



vbo_builder_tristrip::vbo_builder_tristrip( const scene_static::tristrip &tristrip ) 
 : tristrip_(&tristrip)
 {
    glGenBuffers( 3, buffers_ ); check_gl_error;
    glGenBuffers( 1, &index_buffer_ ); check_gl_error;
    
    assert( sizeof(vec3f) == 3 * sizeof( GLfloat ));
    auto &ts = *tristrip_;
    const size_t vertex_size = ts.vecs().size() * 3 * sizeof(GLfloat);
    const size_t color_size = ts.vecs().size() * 4 * sizeof(GLubyte);
    const size_t tex_size = ts.vecs().size() * 2 * sizeof(GLuint);
    
    index_num_ = ts.idx().size();
    const size_t index_size = index_num_ * sizeof( GLuint );
    
    
    
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] ); check_gl_error;
    glBufferData( GL_ARRAY_BUFFER, vertex_size, ts.vecs().data(), GL_STATIC_DRAW ); check_gl_error;
//     glBufferData( GL_ARRAY_BUFFER, -1, scene.strip_vecs().data(), GL_STATIC_DRAW ); check_gl_error;
    
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] ); check_gl_error;
    glBufferData( GL_ARRAY_BUFFER, color_size, 0, GL_DYNAMIC_DRAW ); check_gl_error; check_gl_error;
    
    
    
    {
        // initialuze unused color buffer to prevent crash (unmapped memory?)
    //    glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] );
        GLubyte *b_base = (GLubyte*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY ); check_gl_error;
        std::fill( b_base, b_base + tristrip_->vecs().size() * 4 * sizeof(GLubyte), 255 );
        glUnmapBuffer( GL_ARRAY_BUFFER );
        
    }
    
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[2] ); check_gl_error;
    glBufferData( GL_ARRAY_BUFFER, tex_size, ts.tex().data(), GL_STATIC_DRAW ); check_gl_error;
    
//     glPrimitiveRestartIndex( scene_static::restart_idx ); check_gl_error;
    assert(ts.idx().size() * sizeof(GLuint) == index_size );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_buffer_ ); check_gl_error;
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, index_size, ts.idx().data(), GL_STATIC_DRAW ); check_gl_error; 
    
    
}
void vbo_builder_tristrip::update_color_vec3fptr(const vec3f * const first, const vec3f * const last) {
#if 0
    assert( tristrip_ != nullptr );
    //assert( 0 );
  //  assert( std::distance( first, last ) == ptrdiff_t(tristrip_->idx_pairs().size()) );
    
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] );
    GLubyte *b_base = (GLubyte*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY ); check_gl_error;
    std::fill( b_base, b_base + tristrip_->vecs().size() * 4 * sizeof(GLubyte), 255 );
//     glUnmapBuffer( GL_ARRAY_BUFFER );
//     return;
    
    assert( b_base != nullptr );

    //b += 4 * 3 * num_planes_;
#if 0
    const vec3f *cur = first;
    
    auto &ts = *tristrip_;
    const auto &idx_pairs = ts.idx_pairs();
    //assert( idx_pairs.size() == num_planes_ );
    
    for( ; cur != last; ++cur ) {
        auto idx = std::distance(first, cur);
        assert( idx < ptrdiff_t(idx_pairs.size()) );
        
        auto pair = idx_pairs[idx];
        
        for( uint32_t i = pair.first; i < pair.second; ++i ) {
            GLubyte *b = b_base + i * 4;
            
            *(b++) = gl_utils::clamp<0,255>(255 * cur->r);
            *(b++) = gl_utils::clamp<0,255>(255 * cur->g);
            *(b++) = gl_utils::clamp<0,255>(255 * cur->b);
            *(b++) = 255;
        }

    }
#endif


    glUnmapBuffer( GL_ARRAY_BUFFER );
#endif
}

void vbo_builder_tristrip::draw_arrays() {
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] ); check_gl_error;
    glVertexPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL)); check_gl_error;
    // glColorPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL+ 4 * 3 * num_planes_ * sizeof(GLfloat) ));
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] ); check_gl_error;
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, (GLvoid*)((char*)NULL)); check_gl_error;

    glBindBuffer( GL_ARRAY_BUFFER, buffers_[2] ); check_gl_error;
    glTexCoordPointer(2, GL_FLOAT, 0, (GLvoid*)((char*)NULL)); check_gl_error;

    
//         glBindBuffer(GL_ARRAY_BUFFER, buffer_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_); check_gl_error;
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
// This is the actual draw command
    
//      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//     glDrawElements(GL_TRIANGLE_STRIP, index_num_, GL_UNSIGNED_INT, (GLvoid*)((char*)NULL)); check_gl_error;
    glDrawElements(GL_TRIANGLE_STRIP, index_num_, GL_UNSIGNED_INT, (GLvoid*)((char*)NULL)); check_gl_error;
}


#if 0
vbo_builder::vbo_builder(size_t num_planes) : num_planes_(num_planes) {
    //GLuint b;
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
void vbo_builder::update_index_buffer(size_t num_planes) {
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_buffer_ );
    GLuint *b = (GLuint *) glMapBuffer( GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY );
    assert( b != nullptr );
    for( size_t i = 0; i < num_planes * 4; ++i ) {
        b[i] = i;
    }

    glUnmapBuffer( GL_ELEMENT_ARRAY_BUFFER );

}
void vbo_builder::draw_arrays() {
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

#endif


std::string gl_error_exception::err_str(GLenum err, const char* file, int line) throw() {
    const char *estr = nullptr;
    
    switch( err ) {
        case GL_NO_ERROR: estr = "GL_NO_ERROR"; break;
        case GL_INVALID_ENUM: estr = "GL_INVALID_ENUM"; break;
        case GL_INVALID_VALUE: estr = "GL_INVALID_VALUE"; break;
        case GL_INVALID_OPERATION: estr = "GL_INVALID_OPERATION"; break;
        case GL_STACK_OVERFLOW: estr = "GL_STACK_OVERFLOW"; break;
        case GL_STACK_UNDERFLOW: estr = "GL_STACK_UNDERFLOW"; break;
        case GL_OUT_OF_MEMORY: estr = "GL_OUT_OF_MEMORY"; break;
        default: estr = "unknown";
    };
    
    assert( estr != nullptr );
    
    try {
        std::stringstream ss;
        ss << "opengl error: " << std::hex << int(err) << " (" << estr << ")" << " at " << file << ":" << std::dec << line;

        return ss.str();
    } catch( std::bad_alloc x ) {
        return std::string();
    }
}