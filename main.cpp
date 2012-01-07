/*
 * Copyright (C) 2011 Simon A. Berger
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


#include <iostream>
#include <ClanLib/application.h>
#include <ClanLib/core.h>
#include <ClanLib/sound.h>
#include <ClanLib/gui.h>
#include <ClanLib/display.h>
#include <ClanLib/swrender.h>
#include <ClanLib/gl.h>

#include <cstdio>
#include <fstream>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/array.hpp>
#include <boost/noncopyable.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>

template<class T, std::size_t N>
class farray {
    public:
    T *elems;    // fixed-size array of elements of type T

    public:
    // type definitions
    typedef T              value_type;
    typedef T*             iterator;
    typedef const T*       const_iterator;
    typedef T&             reference;
    typedef const T&       const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    farray( void *b ) : elems((T*)b) {}
    farray( T* b ) : elems(b) {}
    // iterator support
    iterator begin() { return elems; }
    const_iterator begin() const { return elems; }
    iterator end() { return elems+N; }
    const_iterator end() const { return elems+N; }

    // reverse iterator support
#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_MSVC_STD_ITERATOR) && !defined(BOOST_NO_STD_ITERATOR_TRAITS)
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#elif defined(_MSC_VER) && (_MSC_VER == 1300) && defined(BOOST_DINKUMWARE_STDLIB) && (BOOST_DINKUMWARE_STDLIB == 310)
    // workaround for broken reverse_iterator in VC7
    typedef std::reverse_iterator<std::_Ptrit<value_type, difference_type, iterator,
                                    reference, iterator, reference> > reverse_iterator;
    typedef std::reverse_iterator<std::_Ptrit<value_type, difference_type, const_iterator,
                                    const_reference, iterator, reference> > const_reverse_iterator;
#elif defined(_RWSTD_NO_CLASS_PARTIAL_SPEC) 
    typedef std::reverse_iterator<iterator, std::random_access_iterator_tag, 
            value_type, reference, iterator, difference_type> reverse_iterator; 
    typedef std::reverse_iterator<const_iterator, std::random_access_iterator_tag,
            value_type, const_reference, const_iterator, difference_type> const_reverse_iterator;
#else
    // workaround for broken reverse_iterator implementations
    typedef std::reverse_iterator<iterator,T> reverse_iterator;
    typedef std::reverse_iterator<const_iterator,T> const_reverse_iterator;
#endif

    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
    }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
    }

    // operator[]
    reference operator[](size_type i) 
    { 
        BOOST_ASSERT( i < N && "out of range" ); 
        return elems[i];
    }
    
    const_reference operator[](size_type i) const 
    {     
        BOOST_ASSERT( i < N && "out of range" ); 
        return elems[i]; 
    }

    // at() with range check
    reference at(size_type i) { rangecheck(i); return elems[i]; }
    const_reference at(size_type i) const { rangecheck(i); return elems[i]; }

    // front() and back()
    reference front() 
    { 
        return elems[0]; 
    }
    
    const_reference front() const 
    {
        return elems[0];
    }
    
    reference back() 
    { 
        return elems[N-1]; 
    }
    
    const_reference back() const 
    { 
        return elems[N-1]; 
    }

    // size is constant
    static size_type size() { return N; }
    static bool empty() { return false; }
    static size_type max_size() { return N; }
    enum { static_size = N };

   

    // direct access to data (read-only)
    const T* data() const { return elems; }
    T* data() { return elems; }

    // use array as C array (direct read/write access to data)
    T* c_array() { return elems; }

    

    // assign one value to all elements
    void assign (const T& value) { fill ( value ); }    // A synonym for fill
    void fill   (const T& value)
    {
        std::fill_n(begin(),size(),value);
    }

    // check range (may be private because it is static)
    static void rangecheck (size_type i) {
        if (i >= size()) {
            std::out_of_range e("array<>: index out of range");
            boost::throw_exception(e);
        }
    }

    T operator*() {
        return *(begin());
    }
};

template<class T>
class packed {
    T v;
public:    
    packed( char *ptr ) {
        char *tv = (char*)&v;
        std::copy( ptr, ptr + sizeof(T), tv );
    }
    
    T operator*() {
        return v;
    }
};


class CL_API_SOUND CL_SoundProvider_RawNC : public CL_SoundProvider
{
/// \name Construction
/// \{

public:
    /// \brief Constructs a sound provider based on some raw PCM data.
    ///
    /// \param sound_data Raw PCM data.
    /// \param num_samples Number of samples to be read out of sound_data.
    /// \param bytes_per_sample The size of a sample in bytes. This is 2 for 16 bit (signed), and 1 for 8 bit (unsigned).
    /// \param stereo True if sound is stereo (two channels).
    /// \param frequency Playback frequency for sample data.
    CL_SoundProvider_RawNC(
        void *sound_data,
        int num_samples,
        int bytes_per_sample,
        bool stereo,
        int frequency = 22050);

    virtual ~CL_SoundProvider_RawNC();

/// \}
/// \name Operations
/// \{

public:
    /// \brief Called by CL_SoundBuffer when a new session starts.
    /** \return The soundbuffer session to be attached to the newly started session.*/
    virtual CL_SoundProvider_Session *begin_session();

    /// \brief Called by CL_SoundBuffer when a session has finished. After this call,
    /** <p>CL_SoundBuffer will not access the session anymore. It can safely be deleted
        here (and in most cases should be delete here).</p>*/
    virtual void end_session(CL_SoundProvider_Session *session);

/// \}
/// \name Implementation
/// \{

private:
    unsigned char *m_sound_data;

    int m_num_samples;

    int m_bytes_per_sample;

    bool m_stereo;

    int m_frequency;

    friend class CL_SoundProvider_RawNC_Session;
/// \}
};

class CL_SoundProvider_RawNC_Session : public CL_SoundProvider_Session
{
/// \name Construction
/// \{

public:
    CL_SoundProvider_RawNC_Session(CL_SoundProvider_RawNC &source);
    ~CL_SoundProvider_RawNC_Session();


/// \}
/// \name Attributes
/// \{

public:
    int get_num_samples() const;
    int get_frequency() const;
    int get_num_channels() const;
    int get_position() const;


/// \}
/// \name Operations
/// \{

public:
    bool eof() const;
    void stop();
    bool play();
    bool set_position(int pos);
    bool set_end_position(int pos) { return false; }
    int get_data(float **data_ptr, int data_requested);


/// \}
/// \name Implementation
/// \{

private:
    CL_SoundProvider_RawNC source;

    int position;
    int num_samples;
    int frequency;

    bool reached_end;
/// \}
};

CL_SoundProvider_RawNC_Session::CL_SoundProvider_RawNC_Session(CL_SoundProvider_RawNC &source) :
    source(source), position(0), reached_end(false)
{
    frequency = source.m_frequency;
    num_samples = source.m_num_samples;
}

CL_SoundProvider_RawNC_Session::~CL_SoundProvider_RawNC_Session()
{
}

/////////////////////////////////////////////////////////////////////////////
// CL_SoundProvider_Raw_Session attributes:

int CL_SoundProvider_RawNC_Session::get_num_samples() const
{
    return num_samples;
}

int CL_SoundProvider_RawNC_Session::get_frequency() const
{
    return frequency;
}

int CL_SoundProvider_RawNC_Session::get_num_channels() const
{
    return (source.m_stereo) ? 2 : 1;
}

int CL_SoundProvider_RawNC_Session::get_position() const
{
    return position;
}

/////////////////////////////////////////////////////////////////////////////
// CL_SoundProvider_Raw_Session operations:

bool CL_SoundProvider_RawNC_Session::eof() const
{
    return (position >= num_samples);
}

void CL_SoundProvider_RawNC_Session::stop()
{
}

bool CL_SoundProvider_RawNC_Session::play()
{
    return true;
}
    
bool CL_SoundProvider_RawNC_Session::set_position(int pos)
{
    position = pos;
    return true;
}

int CL_SoundProvider_RawNC_Session::get_data(float **data_ptr, int data_requested)
{
    if (position + data_requested > num_samples)
    {
        data_requested = num_samples - position;
        if (data_requested < 0) return 0;
    }

    if (source.m_bytes_per_sample == 2)
    {
        if (source.m_stereo)
        {
            short *src = (short *) source.m_sound_data + position * source.m_bytes_per_sample;
            CL_SoundSSE::unpack_16bit_stereo(src, data_requested, data_ptr);
        }
        else
        {
            short *src = (short *) (source.m_sound_data + position * source.m_bytes_per_sample);
            CL_SoundSSE::unpack_16bit_mono(src, data_requested, data_ptr[0]);
//             for( int i = 0; i < data_requested; i++ ) {
//                 data_ptr[0][i] = src[i] / 32000.0;
//             }
            
            printf( "src: %p %p\n", source.m_sound_data, src  );
        }
    }
    else if (source.m_bytes_per_sample == 1)
    {
        if (source.m_stereo)
        {
            unsigned char *src = (unsigned char *) source.m_sound_data + position * source.m_bytes_per_sample;
            CL_SoundSSE::unpack_8bit_stereo(src, data_requested, data_ptr);
        }
        else
        {
            unsigned char *src = (unsigned char *) source.m_sound_data + position * source.m_bytes_per_sample;
            CL_SoundSSE::unpack_8bit_mono(src, data_requested, data_ptr[0]);
        }
    }

    position += data_requested;
    printf( "get_data: %d %d %d\n", position, data_requested, source.m_bytes_per_sample );
    return data_requested;
}






CL_SoundProvider_RawNC::CL_SoundProvider_RawNC(
    void *sound_data,
    int num_samples,
    int bytes_per_sample,
    bool stereo,
    int frequency)
{
    int data_size = num_samples * bytes_per_sample;
    if (stereo) data_size *= 2;

    m_sound_data = (unsigned char*)sound_data;
    m_num_samples = num_samples;
    m_bytes_per_sample = bytes_per_sample;
    m_stereo = stereo;
    m_frequency = frequency;
}

CL_SoundProvider_RawNC::~CL_SoundProvider_RawNC()
{
}

/////////////////////////////////////////////////////////////////////////////
// CL_SoundProvider_Raw operations:

CL_SoundProvider_Session *CL_SoundProvider_RawNC::begin_session()
{
    return new CL_SoundProvider_RawNC_Session(*this);
}

void CL_SoundProvider_RawNC::end_session(CL_SoundProvider_Session *session)
{
    delete session;
}


class timer {
    double m_start;
public:
    timer() : m_start( CL_System::get_time() / 1000.0 ) {
        
    }
    
    double elapsed() {
     
        return (CL_System::get_time() / 1000.0) - m_start;
    }
    
    
};
float compute_gaussian(float n, float theta) // theta = Blur Amount
{
        return (float)((1.0f / sqrtf(2 * (float)CL_PI * theta)) * expf(-(n * n) / (2.0f * theta * theta)));
}

void setup_gauss( const float blur, const int nsample, CL_Vec2f *offs, float *weight, float dx, float dy ) {
    weight[0] = compute_gaussian( 0, blur );
    offs[0] = CL_Vec2f( 0.0, 0.0 );
    
    float total_weight = weight[0];
    std::cout << "setup gauss\n";
    for( int i = 0; i < nsample / 2; i++ ) {
     
        float w = compute_gaussian( i + 1.0f, blur );
        
        std::cout << "w: " << w << "\n";
        weight[i * 2 + 1] = w;
        weight[i * 2 + 2] = w;
        
        total_weight += w * 2;
        
        float o = i * 2 + 1.5f;
        CL_Vec2f delta( dx * o, dy * o );
        offs[i * 2 + 1] = delta;
        offs[i * 2 + 2] = CL_Vec2d( -delta.x, -delta.y );
        
    }
    
    
    for( int i = 0; i < nsample; i++ ) {
           
        weight[i] /= total_weight;
        
    }
    
}

void setup_gauss_ass( const float blur, const int nsample, CL_Vec2f *offs, float *weight, float dx, float dy ) {
    weight[0] = compute_gaussian( 0, blur );
    offs[0] = CL_Vec2f( 0.0, 0.0 );
    
    float total_weight = weight[0];
    std::cout << "setup gauss\n";
    for( int i = 0; i < nsample; i++ ) {
     
        float w = compute_gaussian( i + 1.0f, blur );
        
        
        weight[i + 1] = w;
        total_weight += w;
        
        float o = i * 2 + 1.5f;
        std::cout << "o: " << w << "\n";
        CL_Vec2f delta( dx * o, dy * o );
        offs[i + 1] = delta;
        //offs[i * 2 + 2] = CL_Vec2d( -delta.x, -delta.y );
        
    }
    
    
    for( int i = 0; i < nsample; i++ ) {
           
        weight[i] /= total_weight;
        
    }
    
}

class Shooter2 {
public:
    static void draw_texture(CL_GraphicContext &gc, const CL_Rectf &rect, const CL_Colorf &color, const CL_Rectf &texture_unit1_coords)
    {
        CL_Vec2f positions[6] =
        {
                CL_Vec2f(rect.left, rect.top),
                CL_Vec2f(rect.right, rect.top),
                CL_Vec2f(rect.left, rect.bottom),
                CL_Vec2f(rect.right, rect.top),
                CL_Vec2f(rect.left, rect.bottom),
                CL_Vec2f(rect.right, rect.bottom)
        };

        CL_Vec2f tex1_coords[6] =
        {
                CL_Vec2f(texture_unit1_coords.left, texture_unit1_coords.top),
                CL_Vec2f(texture_unit1_coords.right, texture_unit1_coords.top),
                CL_Vec2f(texture_unit1_coords.left, texture_unit1_coords.bottom),
                CL_Vec2f(texture_unit1_coords.right, texture_unit1_coords.top),
                CL_Vec2f(texture_unit1_coords.left, texture_unit1_coords.bottom),
                CL_Vec2f(texture_unit1_coords.right, texture_unit1_coords.bottom)
        };

        CL_PrimitivesArray prim_array(gc);
        prim_array.set_attributes(0, positions);
        prim_array.set_attribute(1, color);
        prim_array.set_attributes(2, tex1_coords);
        gc.draw_primitives(cl_triangles, 6, prim_array);
}

    
    
    static int main(const std::vector<CL_String> &args) {
        
        CL_String str = cl_format( "Hello %1 %2", "world", 666 );
        
        printf( "bla: '%s'\n", str.c_str() );
        
        const char *filename = "/space/test2.wav";
        size_t file_size;
        {
            std::ifstream ist(filename);
            ist.seekg(0, std::ios::end);
            file_size = ist.tellg();
        }
        
        
        
        boost::interprocess::file_mapping fm( filename, boost::interprocess::read_only );
        boost::interprocess::mapped_region r = boost::interprocess::mapped_region( fm, boost::interprocess::read_only );
            
        char * base = (char*) r.get_address();
        farray<char, 4> chunk_id(base);
        packed<uint32_t> chunk_size( base + 4 );
        farray<char, 4> format(base + 8);
        
        bool is_riff = std::equal(chunk_id.begin(), chunk_id.end(), "RIFF" );
        printf( "is_riff: %d\n", is_riff );
        
        farray<char, 4>fmt_subchunk_id( base + 12 );
        assert( std::equal( fmt_subchunk_id.begin(), fmt_subchunk_id.end(), "fmt " ));
        
        packed<uint32_t>fmt_subchunk_size( base + 16 );
        packed<uint16_t>audio_format( base + 20 );
        packed<uint16_t>num_channels( base + 22 );
        packed<uint32_t>sample_rate( base + 24 );
        packed<uint32_t>byte_rate( base + 28 );
        packed<uint16_t>blockalign( base + 32 );
        packed<uint16_t>bits_per_sample( base + 34 );
        
        
        
        
        farray<char, 4>data_subchunk_id( base + 36 );
        assert( std::equal( data_subchunk_id.begin(), data_subchunk_id.end(), "data" ));
        
        packed<uint32_t>data_subchunk_size( base + 40 );
        char *data_begin = base + 44;
        char *data_end = data_begin + *data_subchunk_size;
        
        CL_SetupCore setup_core;
        CL_SetupSound setup_sound;
        

        CL_SoundOutput sound_output(44100, 50);

        int bytes_per_sample = ((*bits_per_sample))/ 8;
        int num_samples = *data_subchunk_size / bytes_per_sample;
        bool stereo = (*num_channels) == 2;
        printf( "num channels: %d %d %d %d\n", *num_channels, *bits_per_sample, bytes_per_sample, stereo );
        printf( "num samples: %d\n", num_samples );
        CL_SoundProvider_RawNC sp( data_begin, num_samples, bytes_per_sample, stereo, 44100 );
        CL_SoundBuffer sb( &sp );
        CL_SoundBuffer_Session session = sb.prepare();
        sb.play();
        
        CL_SetupDisplay display;
        //CL_SetupSWRender swrender;
        
        CL_SetupGL setup_gl;
        
        CL_OpenGLWindowDescription desc;
        desc.set_size( CL_Size( 1024, 768 ), true );
        
        CL_DisplayWindow wnd(desc);
        
        CL_GraphicContext gc = wnd.get_gc();
     
        CL_Texture tex_offscreen( gc, gc.get_width(), gc.get_height() );
        tex_offscreen.set_min_filter( cl_filter_linear );
        tex_offscreen.set_mag_filter( cl_filter_linear );
        
        CL_FrameBuffer fb_offscreen( gc );
        fb_offscreen.attach_color_buffer( 0, tex_offscreen );
        
        CL_ProgramObject shader = CL_ProgramObject::load( gc, "vertex_shader.glsl", "fragment_shader.glsl" );
        
        shader.bind_attribute_location(0, "Position");
        shader.bind_attribute_location(1, "Color0");
        shader.bind_attribute_location(2, "TexCoord0");
        if (!shader.link()) {
            throw CL_Exception("Unable to link shader program: Error:" + shader.get_info_log());
        }
        
        CL_Image img1(gc, "/home/sim/art/brick/brick51_1.png" );
        
        timer t1;
        while( true ) {
            gc.set_frame_buffer(fb_offscreen);
            
            
            CL_Draw::line(gc, 0, 0, 100, 100, CL_Colorf( 1.0f ,1.0f ,1.0f ));
            
            img1.draw(gc, 100, 100 );
            
            
            gc.set_texture(0, tex_offscreen);
            gc.reset_frame_buffer();
            gc.set_program_object(shader, cl_program_matrix_modelview_projection);
            shader.set_uniform1i("SourceTexture", 0);
            shader.set_uniform1f("Amount", 1.0);
            shader.set_uniform1f("Timer", 1.0);
            
            shader.set_uniform1f("abber_amount", (sin(t1.elapsed()) * 5 )/ gc.get_height());
            
            const size_t nsample = 15;
            float blur_weight[nsample];
            CL_Vec2f blur_offs[nsample];
            
            const float blur_amount = 5;//(1.0 + sin(t1.elapsed())) * 2.5;
            const float aber_amount = (sin(t1.elapsed()) * 1.0 );
            
            setup_gauss_ass( blur_amount, nsample, blur_offs, blur_weight, 0.0, aber_amount / gc.get_height());
            shader.set_uniformfv("samp_offs", 2, nsample, (float *)blur_offs);
            shader.set_uniformfv("samp_weight", 1, nsample, blur_weight);
            
            
            draw_texture(gc, CL_Rect(0, 0, gc.get_width(), gc.get_height()), CL_Colorf( 1.0f, 1.0f, 1.0f ),  CL_Rectf(0, 0, 1, 1) );//CL_Rectf(0, 0, gc.get_width(), gc.get_height()) );
            gc.reset_program_object();
            gc.reset_texture(0);
            
            wnd.flip();
            
            CL_System::sleep( 10 );
            CL_KeepAlive::process();
    //         session.play();
            //gui.exec(true);
            //usleep( 10000000 );
        }
        
        
        return 0;
    }
    
};

CL_ClanApplication app(&Shooter2::main);
 
