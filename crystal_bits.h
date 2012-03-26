#ifndef __crystal_bits_
#define __crystal_bits_
#include <memory>
#include <boost/numeric/ublas/fwd.hpp>

namespace ublas = boost::numeric::ublas;

class crystal_bits {
public:
    typedef std::unique_ptr<ublas::matrix<int> > matrix_ptr;
    
    static matrix_ptr pump( const ublas::matrix<int> &in, const size_t factor ) ;
    
    static std::vector<std::vector<int>> matrix_to_intvec2d( const ublas::matrix<int> &in ) ;
    static matrix_ptr load_crystal_slice( std::istream &is, size_t width, size_t height ) ;
    static std::vector<matrix_ptr> load_crystal( std::istream &is, size_t pump_factor_ ) ;
        
};
#endif