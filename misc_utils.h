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


#ifndef __misc_utils_h
#define __misc_utils_h

#include <ClanLib/Core/Math/vec4.h>
#include <atomic>
#include <boost/dynamic_bitset.hpp>
#include "cycle.h"

class spinlock_mutex {
public:
    spinlock_mutex() : flag_(ATOMIC_FLAG_INIT) {}
    
    void lock() {
        while( flag_.test_and_set( std::memory_order_acquire ));
    }
    
    void unlock() {
        flag_.clear();
    }
    
private:
    std::atomic_flag flag_;
};

typedef CL_Vec3i vec3i;
typedef CL_Vec3f vec3f;

typedef CL_Vec2i vec2i;
typedef CL_Vec2f vec2f;

struct col3f_sse {
    float r;
    float g;
    float b;
    float x;

    col3f_sse( float r_, float g_, float b_ ) : r(r_), g(g_), b(b_) {}
    col3f_sse( const vec3f & v ) : r(v.r), g(v.g), b(v.b), x(0) {} // { std::cout << "assign: " << r << " " << g << " " << b << "\n";}
    const col3f_sse &operator=( const col3f_sse &other ) {
        r = other.r;
        g = other.g;
        b = other.b;
        x = 0;
        return *this;
    }

    col3f_sse() {}

//  inline operator vec3f() {
//      return vec3f(r, g, b );
//  }

};

template<typename vec_t>
typename vec_t::datatype dist_sqr( const vec_t &v1, const vec_t &v2 ) {
    vec_t vd = v1 - v2;

    return vd.dot(vd);
}



template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

template<typename Block>
class bitset_hash_iterator : public std::iterator<std::output_iterator_tag,void,void,void,void> {
    Block &hash;
    size_t i;
public:
    bitset_hash_iterator ( Block &out_hash ) : hash(out_hash), i(1) { hash = 1234; }
    
    inline bitset_hash_iterator<Block>& operator= (const bitset_hash_iterator<Block> &other ) {
        hash = other.hash;
        return *this;
    }

    inline bitset_hash_iterator<Block>& operator= (const Block &v ) {
        hash ^= v * i++;
        return *this;
    }
    
    inline bitset_hash_iterator<Block>& operator* ()
    { return *this; }
    inline bitset_hash_iterator<Block>& operator++ ()
    { return *this; }
    inline bitset_hash_iterator<Block>& operator++ (int)
    { return *this; }
};

class bitset_hash {
public:
    size_t operator()( const boost::dynamic_bitset<> &bs ) const {
    #ifndef WIN32
        // TODO: find out why the asser fails under 64bit win32
        BOOST_STATIC_ASSERT( sizeof( size_t ) == sizeof( boost::dynamic_bitset<>::block_type ) );
    #endif
        boost::dynamic_bitset<>::block_type hash = 0;
        
        to_block_range( bs, bitset_hash_iterator<boost::dynamic_bitset<>::block_type>(hash));
        
        
        
        // FIXME: what to do when size_t and Block have different width?
        // would a 'static if' with no overhead work?
        return size_t(hash);
    }
};

class tick_timer {
public:
    tick_timer() :t1_(getticks()  ) {}

    double elapsed() {
        return ::elapsed( getticks(), t1_ );
    }

private:
    ticks t1_;
};

template<typename T>
std::vector<T> apply_permutation( std::vector<T> *v1, const std::vector<size_t> &perm ) {
    //assert( v1->size() == v2->size() );
    std::vector<T> v2;
    v2.reserve(perm.size());
    assert( perm.size() == v1->size() );
    for( size_t p : perm ) {
        v2.push_back( std::move( v1->at(p )));
    }
    
    return v2;
}

#endif