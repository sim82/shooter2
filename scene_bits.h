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

#ifndef __scene_bits_h
#define __scene_bits_h

#include <array>
#include <iostream>

#include "crystal_bits.h"
#include "misc_utils.h"
#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/ublas/fwd.hpp>

namespace ublas = boost::numeric::ublas;

class bitmap3d : private boost::dynamic_bitset<>
{

public:
    typedef boost::dynamic_bitset<>::reference reference;

    bitmap3d(size_t x, size_t y, size_t z)
        : boost::dynamic_bitset<>(x * y * z, false)
        , x_(x)
        , y_(y)
        , z_(z)
        , slice_(x * z)
        , stride_(x)
    {
    }

    bitmap3d()
        : x_(0)
        , y_(0)
        , z_(0)
        , slice_(0)
        , stride_(0)
    {
    }

    //  bitmap3d() : bitmap3d(0, 0, 0) // aaahhhhrgggg waiting for gcc 4.7....
    //  {
    //
    //  }

    inline reference operator()(int x, int y, int z)
    {
        if (!inside(x, y, z))
        {
            throw std::runtime_error("index out of bounds");
        }

        return operator[](addr(x, y, z));
    }

    inline bool operator()(int x, int y, int z) const
    {
        if (!inside(x, y, z))
        {
            return false;
        }

        return operator[](addr(x, y, z));
    }

    inline size_t x() const { return x_; }
    inline size_t y() const { return y_; }
    inline size_t z() const { return z_; }

    inline uint64_t hash() const
    {
        bitset_hash bh;

        return bh(*this);
    }

private:
    inline bool inside(int x, int y, int z) const
    {
        return !(x < 0 || y < 0 || z < 0 || x >= int(x_) || y >= int(y_) || z >= int(z_));
    }

    inline size_t addr(size_t x, size_t y, size_t z) const { return x + stride_ * z + slice_ * y; }

    size_t x_;
    size_t y_;
    size_t z_;

    size_t slice_;
    size_t stride_;
};

class util
{
public:
    static bool occluded2(vec3i p0, vec3i p1, const bitmap3d &solid);
    static bool occluded(vec3i p0, vec3i p1, const bitmap3d &solid);
};

struct plane
{
public:
    enum dir_type
    {
        dir_xy_p,
        dir_xy_n,
        dir_yz_p,
        dir_yz_n,
        dir_zx_p,
        dir_zx_n,
    };

    bool normal_cull(const plane &other) const
    {
        //      return (dir_ == dir_xy_p && other.dir_ == dir_xy_n) ||
        //              (dir_ == dir_xy_n && other.dir_ == dir_xy_p) ||
        //              (dir_ == dir_yz_p && other.dir_ == dir_yz_n) ||
        //              (dir_ == dir_yz_n && other.dir_ == dir_yz_p) ||
        //              (dir_ == dir_zx_p && other.dir_ == dir_zx_n) ||
        //              (dir_ == dir_zx_n && other.dir_ == dir_zx_p);

        return dir_ == other.dir_;
    }

    static vec3i normali(dir_type dt);
    static vec3f normal(dir_type dt);

    static vec3f primary(dir_type dt);

    static vec3f primary0(dir_type dt);

    static vec3f primary1(dir_type dt);

    static std::array<float, 4> vgen0(dir_type dt);

    static std::array<float, 4> vgen1(dir_type dt);

    //  static vec3f col_diff( dir_type dt ) {
    //      switch( dt ) {
    //      case dir_xy_p:
    //          return vec3f(1.0, 0.5, 0.0 );
    //
    //      case dir_yz_p:
    //          return vec3f(0.0, 1.0, 0.0 );
    //
    //      case dir_xy_n:
    //      case dir_yz_n:
    //      case dir_zx_p:
    //      case dir_zx_n:
    //      default:
    //          return vec3f(0.8,0.8,0.8);
    //      }
    //  }

    static vec3f col_diff(dir_type dt);

    plane(dir_type d, const vec3f &base_pos, const vec3i &pos, float scale, float energy);

    //     void render() const {
    //
    //         for ( size_t i = 0; i < 4; ++i ) {
    //             //glColor3f( energy_, energy_, energy_ );
    //             glColor3f( energy_rgb_.r, energy_rgb_.g, energy_rgb_.b );
    //             glVertex3fv( verts_[i] );
    //
    //         }
    //     }

    template <typename oiter> void render_vertex(oiter first) const
    {
        for (size_t i = 0; i < 4; ++i)
        {
            *(first++) = verts_[i].x;
            *(first++) = verts_[i].y;
            *(first++) = verts_[i].z;
        }
    }
    template <typename oiter> oiter render_color(oiter first) const
    {
        for (size_t i = 0; i < 4; ++i)
        {
            *(first++) = energy_rgb_.r;
            *(first++) = energy_rgb_.g;
            *(first++) = energy_rgb_.b;
        }

        return first;
    }

    const vec3i &pos() const { return pos_; }

    vec3f norm() const { return normal(dir_); }

    dir_type dir() const { return dir_; }

    void energy(float e) { energy_ = e; }
    void energy_rgb(const vec3f &e) { energy_rgb_ = e; }

    const vec3f &col_diff() const { return col_diff_; }

    const std::array<vec3f, 4> &verts() const { return verts_; }

private:
    dir_type dir_;
    vec3i pos_;

    std::array<vec3f, 4> verts_;

    float energy_;
    vec3f energy_rgb_;
    vec3f col_diff_;
};

class scene_static
{
public:
    typedef std::pair<uint32_t, uint32_t> idx_pair;
    const static uint32_t restart_idx; // = 0xFFFFFFFF;

    scene_static(const vec3f &base_pos)
        : base_pos_(base_pos)
    {
    }
    scene_static() {}

    void init_solid_from_crystal(std::istream &is, size_t pump);
    void init_solid(const std::vector<crystal_bits::matrix_ptr> &slices);

    void init_planes();
    void init_strips();

    const std::vector<plane> &planes() const { return planes_; }

    const bitmap3d &solid() const { return solid_; }

    uint64_t hash() const { return solid_.hash(); }

    const std::vector<vec3f> &strip_vecs() const { return strip_vecs_; }
    const std::vector<uint32_t> &strip_idx() const { return strip_idx_; }

    const std::vector<idx_pair> &strip_idx_pairs() const { return strip_idx_pairs_; }

private:
    std::vector<plane> planes_;
    bitmap3d solid_;
    vec3f base_pos_;

    std::vector<vec3f> strip_vecs_;
    std::vector<uint32_t> strip_idx_;

    std::vector<idx_pair> strip_idx_pairs_;
};

class light_utils
{
public:
    inline static bool normal_cull(const plane &pl1, const plane &pl2)
    {
        const plane::dir_type d1 = pl1.dir();
        const plane::dir_type d2 = pl2.dir();

        const vec3i &p1 = pl1.pos();
        const vec3i &p2 = pl2.pos();

        return p1 == p2 || d1 == d2 || (d1 == plane::dir_xy_n && d2 == plane::dir_xy_p && p1.z < p2.z) ||
               (d1 == plane::dir_xy_p && d2 == plane::dir_xy_n && p1.z > p2.z) ||
               (d1 == plane::dir_yz_n && d2 == plane::dir_yz_p && p1.x < p2.x) ||
               (d1 == plane::dir_yz_p && d2 == plane::dir_yz_n && p1.x > p2.x) ||
               (d1 == plane::dir_zx_n && d2 == plane::dir_zx_p && p1.y < p2.y) ||
               (d1 == plane::dir_zx_p && d2 == plane::dir_zx_n && p1.y > p2.y);
    }

    static void render_light(std::vector<vec3f> *emitptr, const scene_static &scene, const vec3f &light_pos,
                             const vec3f &light_color);
};

class light_static
{
public:
    typedef std::pair<float, int> ff_pair;

    light_static() {}

    light_static(std::vector<std::vector<float>> &&f_fact, std::vector<std::vector<int>> &&f_target)
        : f_fact_(f_fact)
        , f_target_(f_target) /*, f_pairs_(init_pairs(f_fact_, f_target_))*/
    {
    }

    light_static(std::istream &is, uint64_t hash);

    void write(std::ostream &os, uint64_t hash);

    void do_postprocessing()
    {
        // this is the place for strange modifications like pre-multilplying
        // the target indices to turn them into direct addresses into the buffer
        // of 4-wide sse vectors.

        std::for_each(f_target_.begin(), f_target_.end(), [](std::vector<int> &target) {
            std::for_each(target.begin(), target.end(), [](int &t) { t *= 4; });

        });
    }

    size_t num_planes() const { return f_fact_.size(); }

    const std::vector<std::vector<float>> &f_fact() const { return f_fact_; }
    const std::vector<std::vector<int>> &f_target_off4() const { return f_target_; }

private:

    std::vector<std::vector<float>> f_fact_;
    std::vector<std::vector<int>> f_target_;

    //     std::vector<std::vector<ff_pair> > f_pairs_;
};

light_static setup_formfactors(const std::vector<plane> &planes_, const bitmap3d &solid_);

#endif
