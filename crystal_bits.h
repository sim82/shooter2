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

#ifndef __crystal_bits_
#define __crystal_bits_
#include <boost/numeric/ublas/fwd.hpp>
#include <memory>

namespace ublas = boost::numeric::ublas;

class crystal_bits
{
public:
    typedef std::unique_ptr<ublas::matrix<int>> matrix_ptr;

    static matrix_ptr pump(const ublas::matrix<int> &in, const size_t factor);

    static std::vector<std::vector<int>> matrix_to_intvec2d(const ublas::matrix<int> &in);
    static matrix_ptr load_crystal_slice(std::istream &is, size_t width, size_t height);
    static std::vector<matrix_ptr> load_crystal(std::istream &is, size_t pump_factor_);
};
#endif
