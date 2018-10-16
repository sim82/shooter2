#include <sstream>
#include <stdexcept>

#include "gl_bits.h"
#include "scene_bits.h"

vbo_builder_triangles::vbo_builder_triangles(const scene_static &scene)
    : scene_static_(&scene)
    , num_planes_(scene.planes().size())
{
    check_gl_error;

    glGenBuffers(2, buffers_);
    check_gl_error;
    glGenBuffers(1, &index_buffer_);
    check_gl_error;

    assert(sizeof(vec3f) == 3 * sizeof(GLfloat));

    const size_t vertex_size = scene.tri_vecs().size() * 3 * sizeof(GLfloat);
    const size_t color_size  = scene.tri_vecs().size() * 4 * sizeof(GLubyte);

    index_num_              = scene.tri_idx().size();
    const size_t index_size = index_num_ * sizeof(GLuint);

    glGenBuffers(2, buffers_);
    check_gl_error;
    glGenBuffers(1, &index_buffer_);
    check_gl_error;

    glBindBuffer(GL_ARRAY_BUFFER, buffers_[0]);
    check_gl_error;
    glBufferData(GL_ARRAY_BUFFER, vertex_size, scene.tri_vecs().data(), GL_STATIC_DRAW);
    check_gl_error;
    //     glBufferData( GL_ARRAY_BUFFER, -1, scene.strip_vecs().data(), GL_STATIC_DRAW ); check_gl_error;

    glBindBuffer(GL_ARRAY_BUFFER, buffers_[1]);
    check_gl_error;
    glBufferData(GL_ARRAY_BUFFER, color_size, 0, GL_DYNAMIC_DRAW);
    check_gl_error;
    check_gl_error;

    //     glPrimitiveRestartIndex( scene_static::restart_idx ); check_gl_error;
    assert(scene.tri_idx().size() * sizeof(GLuint) == index_size);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
    check_gl_error;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, index_size, scene.tri_idx().data(), GL_STATIC_DRAW);
    check_gl_error;
}
void vbo_builder_triangles::update_color_vec3fptr(const vec3f *const first, const vec3f *const last)
{
    assert(scene_static_ != nullptr);
    // assert( 0 );
    assert(std::distance(first, last) == ptrdiff_t(num_planes_));

    glBindBuffer(GL_ARRAY_BUFFER, buffers_[1]);
    GLubyte *b_base = (GLubyte *)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    check_gl_error;
    //     std::fill( b_base, b_base + 1000, 255 );
    //     glUnmapBuffer( GL_ARRAY_BUFFER );
    //     return;

    assert(b_base != nullptr);

    // b += 4 * 3 * num_planes_;

    const vec3f *cur = first;

    const auto &idx_pairs = scene_static_->tri_idx_pairs();
    assert(idx_pairs.size() == num_planes_);

    for (; cur != last; ++cur)
    {
        auto idx = std::distance(first, cur);
        assert(idx < ptrdiff_t(scene_static_->planes().size()));

        auto pair = idx_pairs[idx];

        for (uint32_t i = pair.first; i < pair.second; ++i)
        {
            GLubyte *b = b_base + i * 4;

            *(b++) = gl_utils::clamp<0, 255>(255 * cur->r);
            *(b++) = gl_utils::clamp<0, 255>(255 * cur->g);
            *(b++) = gl_utils::clamp<0, 255>(255 * cur->b);
            *(b++) = 255;
        }
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
}

void vbo_builder_triangles::draw_arrays(gl_program &prog)
{
#if 0
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[0] ); check_gl_error;
    glVertexPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL));
    // glColorPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL+ 4 * 3 * num_planes_ * sizeof(GLfloat) ));
    glBindBuffer( GL_ARRAY_BUFFER, buffers_[1] ); check_gl_error;
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, (GLvoid*)((char*)NULL));

//         glBindBuffer(GL_ARRAY_BUFFER, buffer_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_); check_gl_error;
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
// This is the actual draw command
    
//     glDrawElements(GL_TRIANGLE_STRIP, index_num_, GL_UNSIGNED_INT, (GLvoid*)((char*)NULL)); check_gl_error;
    glDrawElements(GL_TRIANGLE_STRIP, index_num_, GL_UNSIGNED_INT, (GLvoid*)((char*)NULL)); check_gl_error;
#else

    const GLuint VERTEX_POS_INDEX   = prog.a_position_handle();
    const GLuint VERTEX_COLOR_INDEX = prog.a_color_handle();
    ;

    glEnableVertexAttribArray(VERTEX_POS_INDEX);
    check_gl_error;
    glEnableVertexAttribArray(VERTEX_COLOR_INDEX);
    check_gl_error;

    glBindBuffer(GL_ARRAY_BUFFER, buffers_[0]);
    check_gl_error;
    glVertexAttribPointer(VERTEX_POS_INDEX, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid *)((char *)NULL));
    check_gl_error;

    glBindBuffer(GL_ARRAY_BUFFER, buffers_[1]);
    check_gl_error;
    glVertexAttribPointer(VERTEX_COLOR_INDEX, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, (GLvoid *)((char *)NULL));
    check_gl_error;

    //     glBindAttribLocation( prog.get_program(), VERTEX_POS_INDEX, "a_position" ); check_gl_error;
    //     glBindAttribLocation( prog.get_program(), VERTEX_COLOR_INDEX, "a_color" ); check_gl_error;
    //     std::cout << glGetAttribLocation( prog.get_program(), "a_position" ) << " " << glGetAttribLocation(
    //     prog.get_program(), "a_color" ) << "\n";

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
    check_gl_error;
    glDrawElements(GL_TRIANGLES, index_num_, GL_UNSIGNED_INT, (GLvoid *)((char *)NULL));
    check_gl_error;
#endif
}

vbo_builder::vbo_builder(size_t num_planes)
    : num_planes_(num_planes)
{
    // GLuint b;
    glGenBuffers(2, buffers_);
    glGenBuffers(1, &index_buffer_);

    const size_t vertex_size = num_planes_ * (4 * 3) * sizeof(float);
    glBindBuffer(GL_ARRAY_BUFFER, buffers_[0]);
    glBufferData(GL_ARRAY_BUFFER, vertex_size, 0, GL_STATIC_DRAW);

    color_size = num_planes_ * (4 * 4) * sizeof(GLubyte);
    glBindBuffer(GL_ARRAY_BUFFER, buffers_[1]);
    glBufferData(GL_ARRAY_BUFFER, color_size, 0, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_planes_ * 4 * sizeof(int), 0, GL_STATIC_DRAW);
}
void vbo_builder::update_index_buffer(size_t num_planes)
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
    GLuint *b = (GLuint *)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
    assert(b != nullptr);
    for (size_t i = 0; i < num_planes * 4; ++i)
    {
        b[i] = i;
    }

    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
}
void vbo_builder::draw_arrays()
{
    glBindBuffer(GL_ARRAY_BUFFER, buffers_[0]);
    glVertexPointer(3, GL_FLOAT, 0, (GLvoid *)((char *)NULL));
    // glColorPointer(3, GL_FLOAT, 0, (GLvoid*)((char*)NULL+ 4 * 3 * num_planes_ * sizeof(GLfloat) ));
    glBindBuffer(GL_ARRAY_BUFFER, buffers_[1]);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, (GLvoid *)((char *)NULL));

    //         glBindBuffer(GL_ARRAY_BUFFER, buffer_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    // This is the actual draw command
    glDrawElements(GL_QUADS, num_planes_ * 4, GL_UNSIGNED_INT, (GLvoid *)((char *)NULL));
}

static const char *gl_error_to_string(GLenum err)
{
    switch (err)
    {
    case GL_NO_ERROR:
        return "GL_NO_ERROR";
    case GL_INVALID_ENUM:
        return "GL_INVALID_ENUM";
    case GL_INVALID_VALUE:
        return "GL_INVALID_VALUE";
    case GL_INVALID_OPERATION:
        return "GL_INVALID_OPERATION";
    case GL_OUT_OF_MEMORY:
        return "GL_OUT_OF_MEMORY";
    default:
        return "(unknown)";
    };
}

std::string gl_error_exception::err_str(GLenum err, const char *file, int line) throw()
{
    try
    {
        std::stringstream ss;
        ss << "opengl error: " << gl_error_to_string(err) << " at " << file << ":" << line;

        return ss.str();
    }
    catch (std::bad_alloc x)
    {
        return std::string();
    }
}
gl_program::gl_program(const char *vertex_src, const char *fragment_source)
    : a_position_handle_(0)
    , a_color_handle_(1)
{
    GLuint vertexShader = loadShader(GL_VERTEX_SHADER, vertex_src);
    if (!vertexShader)
    {
        throw std::runtime_error("load vertex shader failed.\n");
    }

    GLuint pixelShader = loadShader(GL_FRAGMENT_SHADER, fragment_source);
    if (!pixelShader)
    {
        throw std::runtime_error("load fragment shader failed.\n");
    }

    program = glCreateProgram();
    if (program)
    {

        // use explicit attribute binding (more or less just for testing. the glGetAttribLocation method below should do
        // exactly the same
        glBindAttribLocation(program, a_position_handle_, "a_position");
        check_gl_error;
        glBindAttribLocation(program, a_color_handle_, "a_color");
        check_gl_error;

        glAttachShader(program, vertexShader);
        check_gl_error;

        glAttachShader(program, pixelShader);
        check_gl_error;

        glLinkProgram(program);
        check_gl_error;
        GLint linkStatus = GL_FALSE;
        glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
        if (linkStatus != GL_TRUE)
        {
            GLint bufLength = 0;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &bufLength);
            if (bufLength)
            {

                std::vector<char> buf(bufLength);
                //                     char* buf = (char*) malloc(bufLength);

                glGetProgramInfoLog(program, bufLength, NULL, buf.data());
                std::cerr << "Could not link program:\n" << buf.data() << "\n";
            }
        }

        //         a_position_handle_ = glGetAttribLocation(program, "a_position");
        //         check_gl_error;

        //         a_color_handle_ = glGetAttribLocation(program, "a_color");
        //         check_gl_error;

        //             LOGI("glGetAttribLocation(\"vPosition\") = %d\n",
        //                  gvPositionHandle);
        //
        gv_mvp_handle = glGetUniformLocation(program, "mvp_matrix");
        check_gl_error;

        std::cerr << "mvp_matrix: " << gv_mvp_handle << "\n";

        // LOGI("glGetAttribLocation(\"mvp_matrix\") = %d\n", gv_mvp_handle);

        color_handle_ = glGetUniformLocation(program, "color");
        check_gl_error;
        //             LOGI("glGetAttribLocation(\"color\") = %d\n", color_handle_);
    }
    else
    {
        throw std::runtime_error("glCreateProgram failed");
    }

    std::cerr << "program created: " << program << "\n";
    glUseProgram(program);
    check_gl_error;
}
gl_program::~gl_program()
{
    // TODO: teardown gl resources
    if (program != 0)
    {
        glDeleteProgram(program);
    }
    program = 0;
}
void gl_program::use()
{
    assert(program != 0);
    check_gl_error;
    // std::cerr << "program: " << program << "\n";
    glUseProgram(program);
    check_gl_error;
}

GLuint gl_program::uniform_handle(const char *name)
{
    // maybe this is not faster than calling glGetUniformLocation, but at least it is guaranteed not to be slow...

    std::vector<name_handle_pair>::iterator it = std::lower_bound(uniform_handles_.begin(), uniform_handles_.end(),
                                                                  name, compare_first_string<name_handle_pair>());

    if (it != uniform_handles_.end() && it->first == name)
    {
        return it->second;
    }

    GLuint h = glGetUniformLocation(program, name);

    if (h == GLuint(-1))
    {
        std::stringstream ss;
        ss << "glGetUniformLocation failed: " << name;
        throw std::runtime_error(ss.str());
    }

    bool need_sort = (it != uniform_handles_.end());
    uniform_handles_.push_back(std::make_pair(name, h));
    if (need_sort)
    {
        std::sort(uniform_handles_.begin(), uniform_handles_.end(), compare_first_string<name_handle_pair>());
    }

    return h;
}
GLuint gl_program::loadShader(GLenum shaderType, const char *pSource)
{
    GLuint shader = glCreateShader(shaderType);
    if (shader)
    {
        glShaderSource(shader, 1, &pSource, NULL);
        glCompileShader(shader);
        GLint compiled = 0;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
        if (!compiled)
        {
            GLint infoLen = 0;
            glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
            if (infoLen)
            {
                char *buf = (char *)malloc(infoLen);
                if (buf)
                {
                    glGetShaderInfoLog(shader, infoLen, NULL, buf);
                    check_gl_error;
                    std::cout << "shader info log: \n" << buf;
                    free(buf);
                }
                glDeleteShader(shader);
                shader = 0;
            }
        }
    }
    return shader;
}
