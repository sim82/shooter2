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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "misc_utils.h"

#include "crystal_bits.h"
crystal_bits::matrix_ptr crystal_bits::pump(const ublas::matrix< int >& in, const size_t factor) {
    ublas::matrix<int> out( in.size1() * factor, in.size2() * factor );

    for( size_t row = 0; row != in.size1(); ++row ) {
        ublas::matrix_row<ublas::matrix<int>> r( out, row * factor );



        for( size_t col = 0; col != in.size2(); ++col ) {
            for( size_t i = 0; i < factor; ++i ) {
                r[col*factor+i] = in(row, col);
            }
        }
        for( size_t i = 0; i < factor; ++i ) {
            ublas::matrix_row<ublas::matrix<int>> r2( out, row * factor + i);
            std::copy( r.begin(), r.end(), r2.begin() );
        }

    }

    return make_unique<ublas::matrix<int>>( std::move(out) );
}
std::vector<std::vector<int>> crystal_bits::matrix_to_intvec2d(const ublas::matrix< int >& in) {
    std::vector<std::vector<int>> out;
    out.reserve(in.size1());

    for( auto it1 = in.begin1(); it1 != in.end1(); ++it1 ) {
        out.emplace_back( it1.begin(), it1.end() );

        //             std::cout << "len2: " << out.back().size() << "\n";
    }
    //         std::cout << "len1: " << out.size() << "\n";
    return out;
}

crystal_bits::matrix_ptr crystal_bits::load_crystal_slice(std::istream& is, size_t width, size_t height) {
    ublas::matrix<int> slice( height, width );

    auto mapc = [](char c) {
        c = std::tolower(c);

        if ( c == ' ' ) {
            return int(0);
        } else if ( c >= 'a' && c <= 'z' ) {
            return int(1 + c - 'a');
        } else if ( c >= '0' && c <= '9' ) {
            return int(2 + 'z' - 'a' + c - '0');
        } else {
            std::cerr << "bad: " << int(c) << "\n";
            throw std::runtime_error( "bad character in map");
        }
    };


    auto it1 = slice.begin1();
    for( size_t i = 0; i < height; ++i, ++it1 ) {
        auto it2 = it1.begin();
        for( size_t j = 0; j < width; ++j, ++it2 ) {
            *it2 = mapc( is.get() );
        }

        int nl = is.get();
        assert( nl == '\n' );
    }

    return make_unique<ublas::matrix<int>>( std::move(slice));
}
std::vector< crystal_bits::matrix_ptr > crystal_bits::load_crystal(std::istream& is, size_t pump_factor_) {
    size_t width;
    size_t height;
    size_t num;

    is >> width;
    is >> height;
    is >> num;

    while( is.get() != '\n' ) {}

    std::cout << "size: " << width << " " << height << " " << num << "\n";

    std::vector<crystal_bits::matrix_ptr > out;

    for( size_t i = 0; i < num; ++i ) {
        crystal_bits::matrix_ptr mp = load_crystal_slice( is, width, height );
        out.emplace_back(pump(*mp, pump_factor_));

//         std::cout << "map: " << i << "\n";

//         for( auto it1 = out.back()->begin1(); it1 != out.back()->end1(); ++it1 ) {
//             std::copy( it1.begin(), it1.end(), std::ostream_iterator<int>(std::cout, " " ));
//             std::cout << "\n";
//         }


    }


    std::cout << "size: " << out.size() << "\n";

    return out;

    //throw "exit";

    //         size_t len = size_t(-1);
    //         while ( !is.eof() ) {
    //             std::string line;
    //
    //             std::getline(is, line);
    //
    // //          while( !is.eof() ) {
    // //              char c = is.get();
    // //              if( c == '\n' ) {
    // //                  break;
    // //              }
    // //              line.push_back(c);
    // //          }
    //
    //             if ( line.empty()) {
    //                 break;
    //             }
    //
    //             std::cout << "len: " << line.size() << "'" << std::string(line.begin(), line.end()) << "'\n";
    //
    //
    //             if ( len == size_t(-1)) {
    //                 len = line.size();
    //             } else {
    //                 assert( len == line.size() );
    //             }
    //
    //             ret.push_back(std::vector<int>(len));
    //
    //             std::transform( line.begin(), line.end(), ret.back().begin(), [](char c) {
    //                 c = std::tolower(c);
    //
    //                 if ( c == ' ' ) {
    //                     return int(0);
    //                 } else if ( c >= 'a' && c <= 'z' ) {
    //                     return int(1 + c - 'a');
    //                 } else if ( c >= '0' && c <= '9' ) {
    //                     return int(2 + 'z' - 'a' + c - '0');
    //                 } else {
    //                     std::cerr << "bad: " << int(c) << "\n";
    //                     throw std::runtime_error( "bad character in map");
    //                 }
    //             });
    //
    //         }




}
