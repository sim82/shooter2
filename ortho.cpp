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

#include <ClanLib/application.h>
#include <ClanLib/core.h>
#include <ClanLib/display.h>
#include <ClanLib/gl.h>
#include <GL/gl.h>
#include <GL/glx.h>

#ifdef TEST_OPENCL
#include <CL/cl_gl.h>
#define __CL_ENABLE_EXCEPTIONS
#include "cl.hpp"
#include "cl_error_codes.h"
#include <CL/cl_gl.h>
#endif

#include "cycle.h"
#include <algorithm>
#include <array>
#include <execinfo.h>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "aligned_buffer.h"
#include "crystal_bits.h"
#include "gl_bits.h"
#include "misc_utils.h"
#include "player_bits.h"
#include "rad_core.h"
#include "scene_bits.h"
#include "vec_unit.h"

namespace ublas = boost::numeric::ublas;

static const char gVertexShader[] = "uniform mat4 mvp_matrix;\n"
                                    //     "attribute vec4 xxx;\n"
                                    "attribute vec4 a_position;\n"
                                    "attribute vec4 a_color;\n"

                                    //     "void main() {\n"
                                    //     "  gl_Position = mvp_matrix * gl_Vertex;\n"
                                    //     "  gl_FrontColor = gl_Color;\n"
                                    //     "}\n";

                                    "void main() {\n"
                                    "  gl_Position = mvp_matrix * a_position;\n"
                                    "  gl_FrontColor = a_color;\n"
                                    "}\n";
static const char gFragmentShader[] =
    //"precision mediump float;\n"
    //"uniform vec4 color;\n"
    //"attribute vec4 color;\n"
    "void main() {\n"
    //  "  gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);\n"
    "  gl_FragColor = gl_Color;\n"
    "}\n";

template <typename IDX> std::vector<IDX> ff_size_sort_permutation(const std::vector<std::vector<float>> &ffs)
{
    typedef std::pair<size_t, size_t> ss_pair;
    std::vector<ss_pair> fs;
    for (size_t i = 0; i < ffs.size(); ++i)
    {
        // std::cout << "ffs: " << i << " " << ff2s_[i].size() << "\n";
        fs.emplace_back(ffs[i].size(), i);
    }

    std::sort(fs.begin(), fs.end(), [](const ss_pair &p1, const ss_pair &p2) {
        return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first;
    });

    std::vector<IDX> perm_size(fs.size());
    perm_size.resize(fs.size());

    std::transform(fs.begin(), fs.end(), perm_size.begin(), [](const ss_pair &p) { return p.second; });

    return perm_size;
}

class light_dynamic
{
public:
    light_dynamic() {}
    light_dynamic(size_t num)
        : emit_(num)
        , rad_(num)
    {
    }

    void clear_emit() { std::fill(emit_.begin(), emit_.end(), vec3f(0.0, 0.0, 0.0)); }

    std::vector<vec3f> *emit() { return &emit_; }

    std::vector<vec3f> *rad() { return &rad_; }

private:
    std::vector<vec3f> emit_;
    std::vector<vec3f> rad_;
};

std::string hash_to_filename(uint64_t hash)
{
    std::stringstream ss;
    ss << "baked";
    for (size_t i = 0; i < 8; ++i)
    {
        size_t c = hash & 0xff;

        ss << std::hex << c;

        hash >>= 1;
    }
    ss << ".bin";
    return ss.str();
}


class render_unit
{
public:
    render_unit(std::istream &is, const vec3i &base_pos)
        : base_pos_(base_pos)
        , scene_static_(base_pos)
    {
        check_gl_error;

        assert(is.good());
        const size_t pump_factor = 4;

        base_pos_.x *= pump_factor;
        base_pos_.z *= pump_factor;
        scene_static_.init_solid_from_crystal(is, pump_factor);

        //        scene_static_.init_planes();
        scene_static_.init_tris();
        uint64_t scene_hash = scene_static_.hash();
        auto bin_name       = hash_to_filename(scene_hash);

        std::cout << "baked name: " << bin_name << "\n";
        try
        {
            std::ifstream is(bin_name.c_str());

            light_static_ = light_static(is, scene_hash);
        }
        catch (std::runtime_error x)
        {

            std::cerr << "load failed. recreating. error:\n" << x.what() << std::endl;

            light_static_ = setup_formfactors(scene_static_.planes(), scene_static_.solid());
        }

        if (!false)
        {
            std::ofstream os(bin_name.c_str());
            light_static_.write(os, scene_hash);
        }

        light_static_.do_postprocessing();
        check_gl_error;

        light_dynamic_ = light_dynamic(scene_static_.planes().size());
        rad_core_      = makeRadCoreThreaded(scene_static_, light_static_);

        //         vbob_ = vbo_builder(scene_static_.planes().size() );
        //         vbob_.update_index_buffer( scene_static_.planes().size());
        //         vbob_.update_vertices( scene_static_.planes().begin(), scene_static_.planes().end());
        //
        vbob_ts_ = vbo_builder_triangles(scene_static_);
    }

    void clear_emit() { light_dynamic_.clear_emit(); }
    void render_light(const vec3f &pos, const vec3f &color)
    {
        light_utils::render_light(light_dynamic_.emit(), scene_static_, pos - vec3f(base_pos_), color);
    }

    void update()
    {

        rad_core_->set_emit(*light_dynamic_.emit());
        rad_core_->copy(light_dynamic_.rad());

        //         vbob_.update_color( light_dynamic_.rad()->begin(), light_dynamic_.rad()->end());
        vbob_ts_.update_color(light_dynamic_.rad()->begin(), light_dynamic_.rad()->end());
    }
    void draw(gl_program &prog)
    {
        //         vbob_.draw_arrays();
        vbob_ts_.draw_arrays(prog);
    }

private:
    vec3i base_pos_;

    scene_static scene_static_;

    light_static light_static_;
    light_dynamic light_dynamic_;

    std::unique_ptr<IRadCore> rad_core_;

    //     vbo_builder vbob_;
    vbo_builder_triangles vbob_ts_;
};


class ortho
{

public:
    ortho()
        : pump_factor_(4)
    {
        CL_OpenGLWindowDescription desc;
        desc.set_size(CL_Size(1024, 768), true);

        desc.set_depth_size(16);

        wnd_ = CL_DisplayWindow(desc);

        CL_GraphicContext_GL gc = wnd_.get_gc();

        gc.set_active();

        prog1_ = gl_program(gVertexShader, gFragmentShader);

        glMatrixMode(GL_PROJECTION); // hello

        CL_Mat4f proj = CL_Mat4f::perspective(60, 1.5, 2, 200);
        glLoadMatrixf(proj.matrix);

        struct texel
        {
            GLubyte col[3];
            GLubyte alpha;

            texel()
            {
                col[0] = 128;
                col[1] = 128;
                col[2] = 128;
                alpha  = 255;
            }
        };

        std::array<texel, 64 * 64> tex_data;

        glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 0.0);
        glMatrixMode(GL_MODELVIEW);

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDepthFunc(GL_LESS);

        glShadeModel(GL_FLAT);
    }

    CL_Mat4f setup_perspective(const player &camera)
    {
        CL_Mat4f proj_p = CL_Mat4f::perspective(60, 1.5, 0.2, 500);
        const vec3f &player_pos = camera.pos();
        CL_Mat4f proj_mv        = CL_Mat4f::translate(-player_pos.x, -player_pos.y, -player_pos.z) *
                           CL_Mat4f::rotate(CL_Angle(-camera.rot_x(), cl_degrees),
                                            CL_Angle(-camera.rot_y(), cl_degrees), CL_Angle(), cl_XYZ);

        return proj_mv * proj_p;
    }

    CL_Mat4f setup_ortho()
    {

        CL_Mat4f proj_p = CL_Mat4f::ortho(-20.0, 20.0, -15.0, 15.0, 0, 200);

        CL_Mat4f proj_mv = CL_Mat4f::look_at(10, 10, 10, 0, 0, 0, 0.0, 1.0, 0.0);

        return proj_mv * proj_p;
    }

    void start()
    {

        float x1 = -10;
        float y1 = -10;

        int light_x  = 0;
        int light_xd = 1;

        vec3i light_pos(0, 40, 0);

        auto t_old             = CL_System::get_microseconds();
        bool light_on          = true;
        bool light_button_down = false;
        double delta_t         = 0.01;
        player p1;

        std::ifstream is("cryistal-castle-hidden-ramp.txt");
        //         std::ifstream is( "house1.txt" );
        render_unit runit(is, vec3f(-40.0, -20.0, -40.0));

        bool light_changed = true;
        bool do_quit = false;

        while (!do_quit)
        {

            // cube c(x1, 0, y1);

            CL_GraphicContext gc     = wnd_.get_gc();
            CL_InputDevice &keyboard = wnd_.get_ic().get_keyboard();
            CL_InputDevice &mouse    = wnd_.get_ic().get_mouse();

            p1.input(keyboard, mouse);
            p1.frame(t_old * 1.0e-3, delta_t);

            int light_speed = 1;

            do_quit = keyboard.get_keycode(CL_KEY_Q);

            if (keyboard.get_keycode(CL_KEY_LEFT))
            {
                light_pos.x += light_speed;
                light_changed = true;
            }
            if (keyboard.get_keycode(CL_KEY_RIGHT))
            {
                light_pos.x -= light_speed;
                light_changed = true;
            }
            if (keyboard.get_keycode(CL_KEY_UP))
            {
                light_pos.z += light_speed;
                light_changed = true;
            }
            if (keyboard.get_keycode(CL_KEY_DOWN))
            {
                light_pos.z -= light_speed;
                light_changed = true;
            }
            if (keyboard.get_keycode(CL_KEY_L))
            {
                if (!light_button_down)
                {
                    light_on      = !light_on;
                    light_changed = true;
                }
                light_button_down = true;
            }
            else
            {
                light_button_down = false;
            }

            glEnable(GL_CULL_FACE);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (light_x > 60 || light_x < -60)
            {
                light_xd = -light_xd;
            }

            if (light_changed)
            {
                runit.clear_emit();
                if (light_on)
                {
                    vec3f light_weird = (vec3f(light_pos) * float(pump_factor_)) - base_pos_;
                    runit.render_light(light_weird, vec3f(1.0, 0.8, 0.6));
                }
                light_changed = false;
            }
            runit.update();
            light_x += light_xd;

            int ortho_width  = 320;
            int ortho_heigth = 200;
            glViewport(gc.get_width() - ortho_width, gc.get_height() - ortho_heigth, ortho_width, ortho_heigth);
            auto mat_mvp = setup_perspective(p1);
            prog1_.use();
            glUniformMatrix4fv(prog1_.mvp_handle(), 1, GL_FALSE, mat_mvp.matrix);
            check_gl_error;

            glViewport(0, 0, gc.get_width(), gc.get_height());

            runit.draw(prog1_);
            wnd_.flip(1);

            auto t = CL_System::get_microseconds();

            delta_t = (t - t_old) * 1.0e-6;
            t_old   = t;
            CL_KeepAlive::process();

            x1 += 1;
            if (x1 > 10)
            {
                x1 = -10;
                y1 += 1;
            }
            if (y1 > 10)
            {
                y1 = -10;
            }
        }
    }

    static int main(const std::vector<CL_String> &args)
    {
        _mm_setcsr(_mm_getcsr() | _MM_FLUSH_ZERO_ON);
        //      plane p( plane::dir_zx_p, vec3f( 0.5, 0.5, 0.5 ));

        //      return 0;
        try
        {

            ortho o;
            o.start();
        }
        catch (gl_error_exception x)
        {
            std::cerr << x.what() << std::endl;
            std::cerr << "bailing out\n";
        }

        return 0;
    }

private:
    CL_SetupCore setup_core_;

    CL_SetupDisplay display_;
    // CL_SetupSWRender swrender;

    CL_SetupGL setup_gl_;
    CL_DisplayWindow wnd_;
    GLuint texName;

    // std::vector<crystal_bits::matrix_ptr> height_fields_;

    //

    const size_t pump_factor_;
    vec3f base_pos_;

    gl_program prog1_;

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
