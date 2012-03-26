#ifndef __gl_bits_h
#define __gl_bits_h

#include <ClanLib/gl.h>


class vbo_builder {
  
    
public:
    vbo_builder( size_t num_planes ) : num_planes_(num_planes) {
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
#endif