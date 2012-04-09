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


#include <iostream>

#include <thread>
#include <mutex>
#include <atomic>
#include <ClanLib/core.h>

#include "scene_bits.h"
#include "rad_core.h"
#include "vec_unit.h"
#include "aligned_buffer.h"
#include "misc_utils.h"

class join_threads {
  
public:
    join_threads( std::vector<std::thread> *threads ) : threads_( threads ) {}
    
    ~join_threads() {
        for( size_t i = 0; i < threads_->size(); ++i ) {
            if( (*threads_)[i].joinable() ) {
                (*threads_)[i].join();
            }
        }
    }
    
private:
    std::vector<std::thread> * const threads_;
    
};

class rad_core_threaded: public rad_core {
    typedef std::mutex lock_type;
    //typedef spinlock_mutex lock_type;
public:
    class worker {

    };

    
    std::vector<std::pair<size_t,size_t>> calc_plane_distribution( const size_t num_partition ) {
        size_t num_ints = 0;
        std::vector<std::pair<size_t,size_t>> parts;
        for( auto &ff : light_static_.f_fact() ) {
            num_ints += ff.size();
        }
        
        std::cout << "num_ints: " << num_ints << "\n";
        
        size_t ints_per_part = num_ints / num_partition;
        
        size_t first = 0;
        size_t last = 0;
        size_t acc = 0;
        for( size_t i = 0; i < num_partition; ++i ) {
            while( last < light_static_.f_fact().size() && acc < ints_per_part ) {
                acc += light_static_.f_fact()[last].size();
                ++last;
            }
            acc = 0;
            parts.emplace_back( first, last );
            first = last;
        }
        
        parts.back().second = light_static_.f_fact().size();
        
        
        for( auto &p : parts ) { 
            std::cout << "part: " << p.first << " " << p.second << "\n";
        }
        return parts;
    }

    rad_core_threaded( const scene_static &scene_static, const light_static &light_static /*, const std::vector<std::pair<size_t,size_t>> &part_bounds */)
            : 
            scene_static_(scene_static),
            light_static_(light_static),
            rad_is_new_(false),
            emit_is_new_(false),
            emit_new_(light_static.num_planes()), emit_( light_static.num_planes() ), rad_( light_static.num_planes() ), rad2_( light_static.num_planes() ),
            //ffs_(ffs), ff_target_(ff_target),
            //planes_(planes),
            pints_(0),
            pints_last_(0),
            pints_last_time_(0),
            start_flag_(false),
            abort_flag_(false),
            joiner_( &threads_ )
    {

        
        
        if ( !true ) {
            threads_.push_back( std::thread( [&]() {
                work(0, scene_static_.planes().size(), 0);
            }) );

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
        } else {
#if 1
            const size_t num_planes = scene_static_.planes().size();
            const size_t num_threads = 2;
            
            auto part = calc_plane_distribution(num_threads);
            for ( size_t i = 0; i < num_threads; ++i ) {
                
                
                
                //const size_t first = num_planes / num_threads * i;
                //const size_t last = num_planes / num_threads * (i+1);
                const size_t first = part.at(i).first;
                const size_t last = part.at(i).second;
            
                std::cout << "thread: " << i << " " << first << " " << last << " " << num_planes << "\n";
                
                threads_.push_back( std::thread( [=]() {
                    work( first, last, i );
                }));
            }
#else
            
            const size_t num_threads = part_bounds.size();
            
            
            for ( size_t i = 0; i < num_threads; ++i ) {
                
                
                
                //const size_t first = num_planes / num_threads * i;
                //const size_t last = num_planes / num_threads * (i+1);
                const size_t first = part_bounds.at(i).first;
                const size_t last = part_bounds.at(i).second;
            
                std::cout << "thread: " << i << " " << first << " " << last <<  "\n";
                
                threads_.push_back( std::thread( [=]() {
                    work( first, last, i );
                }));
            }
#endif

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
//          thread1_ = std::thread( [&]() { work(planes_.size()/4, planes_.size()/4 * 2);
//          });
//          thread2_ = std::thread( [&]() { work(planes_.size()/4 * 2, planes_.size()/4 * 3);
//          });
//          thread3_ = std::thread( [&]() { work(planes_.size()/4 * 3, planes_.size()/4 * 4);
//          });
        }
        
        start_flag_.store(true);
    }

    virtual ~rad_core_threaded() {
        std::cerr << "~rad_core_threaded: request thread abort" << std::endl;
        std::cerr << "if the program hangs after this point, there is an error in the thread lifecycle management..." << std::endl;
        
        abort_flag_.store(true);
    }
    
    virtual void set_emit( const std::vector<vec3f> &emit ) {
        std::lock_guard<lock_type>lock(mtx_);

//         std::cout << "emit: " << emit_new_.size() << " " << emit.size() << "\n";
        
        assert( emit_new_.size() == emit.size() );

        emit_new_.assign( std::begin(emit), std::end(emit));
        emit_is_new_ = true;
    }

    virtual bool update() {
        return rad_is_new_;
//         std::lock_guard<std::mutex> lock(mtx_);
// 
//         if ( rad_is_new_ ) {
//             rad_is_new_ = false;
//             return true;
//         } else {
//             return false;
//         }


    }

    virtual void copy( std::vector<vec3f> *out ) {
        std::lock_guard<lock_type> lock(mtx_);
        for ( size_t i = 0; i < rad_.size(); ++i ) {
            (*out)[i].r = rad_[i].r;
            (*out)[i].g = rad_[i].g;
            (*out)[i].b = rad_[i].b;
        }

        rad_is_new_ = false;
        cl_ubyte64 time = CL_System::get_microseconds();
        
        cl_ubyte64 dt = time - pints_last_time_;
        
        if( dt >= 1e6 ) {
            const size_t dp = pints_ - pints_last_;
            pints_last_ = pints_;
            
            pints_last_time_ = time;
            
            std::cout << "pint/s: " << dp / (dt / 1.0e6) << "\n";
        }
    }




private:

    void work( size_t first, size_t last, size_t rank ) {
        while( !start_flag_.load() && !abort_flag_.load() ) {}
        
        
        hrtimer t1;
        uint64_t last_pints = 0;
        
        while (!abort_flag_.load()) {
            if( rank == 0 )
            {
                std::lock_guard<lock_type> lock(mtx_);

                if ( emit_is_new_ ) {
                    //emit_.swap(emit_new_);
                    //emit_ = std::move(emit_new_);
                    assert( emit_new_.size() == emit_.size());
                    std::copy( emit_new_.begin(), emit_new_.end(), emit_.begin() );
                    emit_is_new_ = false;
                }

            }


//             tick_timer tt;
            do_radiosity_sse(first, last);

//             std::cout << "elapsed: " << tt.elapsed() << "\n";

            if( rank == 0 )
            {
                std::lock_guard<lock_type> lock(mtx_);
                rad_is_new_ = true;
                
                double td = t1.elapsed();
                
                if( td > 1.6e9 ) {
                    auto dpints = pints_ - last_pints;
                    
                    last_pints = pints_;
                    
                    std::cout << "ticks: " << td << " " << dpints << " " << td / dpints << "\n";
                    
                    t1.reset();
                    
                }
                
            }

        }
        std::cerr << "rad_core_threaded: worker thread exit: " << rank << std::endl;
    }

    void do_radiosity_sse( size_t first, size_t last ) {
        const int steps = 1;
        const float min_ff = 0;
//      std::fill(e_rad_sse_.begin(), e_rad_sse_.end(), vec3f(0.0, 0.0, 0.0));
//      std::fill(e_rad2_sse_.begin(), e_rad2_sse_.end(), vec3f(0.0, 0.0, 0.0));

        //std::copy( emit_sse_.begin(), emit_sse_.end(), e_rad_sse_.begin() );

        const std::vector<std::vector<float> > &ffs_ = light_static_.f_fact();
        const std::vector<std::vector<int> > &ff_target_ = light_static_.f_target_off4();
//         auto &ff_pairs = light_static_.f_pairs();
        
        const std::vector<plane> &planes_ = scene_static_.planes();
        
        typedef vector_unit<float,4> vu;
        typedef vu::vec_t vec_t;
        //steps = 0;
//         std::cout << "do: " << first << " " << last << "\n";
        vec_t reflex_factor = vu::set1(1.0);
        
        
        const size_t i_unroll = 2;
        size_t pints = 0;
        for ( int i = 0; i < steps; ++i ) {
            //for( auto it = pairs_.begin(); it != pairs_.end(); ++it, ++ff_it ) {
                
                
            for ( size_t j = first; j < last; ++j ) {
                    
                //                 const size_t s = (ffs_[j].size() / 4) * 4;
                
                const bool unroll = !false;
                const size_t send_unroll = (ffs_[j].size() / i_unroll) * i_unroll;
                const size_t send = ffs_[j].size();
                size_t sstart = 0;
                vec_t rad = vu::set1(0);
                vec_t rad2 = vu::set1(0);
                vec_t rad3 = vu::set1(0);
                vec_t rad4 = vu::set1(0);
                
                const vec3f cd = planes_[j].col_diff();
                const vec_t col_diff = vu::set( 0, cd.b, cd.g, cd.r );

                auto &targets = ff_target_[j];
                auto &ffs = ffs_[j];
//                 auto &pairs = ff_pairs[j];
                
                _mm_prefetch( (char*)targets.data(), _MM_HINT_T0 );
                _mm_prefetch( (char*)ffs.data(), _MM_HINT_T0 );
#if 0        
                if( unroll ) {
                    for ( size_t k = 0; k < send_unroll; k+=i_unroll ) {
                        const vec_t ff = vu::set1( ffs[k] );
			size_t target = targets[k];
                        rad = vu::add( rad, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target))), ff ));
                        
                        const vec_t ff2 = vu::set1( ffs[k+1] );
			
			size_t target2 = targets[k+1];
                        rad2 = vu::add( rad2, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target2))), ff2 ));
			
			
#if 0
			const vec_t ff3 = vu::set1( ffs[k+2] );
                        
			size_t target3 = targets[k+2];
			rad3 = vu::add( rad3, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target3))), ff3 ));
                        const vec_t ff4 = vu::set1( ffs[k+3] );
                        
			size_t target4 = targets[k+3];
                        //auto p = pairs[k];
                        //size_t target = p.second;
                        //rad_rgb += (col_diff * e_rad_rgb_[target]) * ff2s_[j][k];
                        
                        //                     if ( false && ffs[k] < min_ff ) {
                            //                         continue;
                        //                     }
                        
                        //const vec_t ff = vu::set1( p.first );
                        rad4 = vu::add( rad4, vu::mul( vu::mul( col_diff, vu::load( (float*) rad_(target4))), ff4 ));
#endif                   
                    }
                    sstart = send_unroll;
                }
#endif          

                const float * const rad_base = (float*)rad_(0);
                for ( size_t k = sstart; k < send; ++k ) {
                    size_t target = targets[k]; // target hast to be pre-multiplied by 4 in light_static::do_preporcessing!
                    
                    const vec_t ff = vu::set1( ffs[k] );
//                     const vec_t ff_diff = ;
                    //const vec_t cd_ff = vu::mul( col_diff, ff);
//                     rad = vu::add( rad, vu::mul( vu::load( (float*) rad_(target)), vu::mul( ff, col_diff ) ));
                    const vec_t tcol = vu::load( rad_base + target );
#if 1                    
                     rad = vu::add( rad, vu::mul( tcol, vu::mul( ff, col_diff )));
#else
                    
                    // this version seems to be a bit faster on core2 class cpus but is slower in sandy bridge...
                    rad = vu::add( rad, vu::mul (vu::mul( tcol, ff ), col_diff ) );
#endif
                }
                
                
                pints += send;
#if 0         
                if( unroll ) {
                    rad = vu::add( rad, rad2 );
#if 0
		    rad3 = vu::add( rad3, rad4 );
                    rad = vu::add( rad, rad3 );
#endif
                }
#endif
                vu::store( vu::add( vu::load((float*)emit_(j)), vu::mul(rad, reflex_factor)), (float*)rad2_(j));
//              std::cout << "col: " << rad2_[j].r << " " << cd.r <<  "\n";
                //e_rad2_rgb_[j] = emit_rgb_[j] + rad_rgb;// * reflex_factor;

            }


            {
                std::lock_guard<lock_type>lock(mtx_);
                //rad_.swap(rad2_);
                std::copy( &rad2_[first], &rad2_[last], &rad_[first]);
                pints_ += pints;
            }
            //rad_ = rad2__;


        }
        
        
        //e_rad_rgb_.assign( e_rad_sse_.begin(), e_rad_sse_.end() );
        

    }

    lock_type mtx_;
//  std::thread thread0_;
//  std::thread thread1_;
//  std::thread thread2_;
//  std::thread thread3_;

    std::vector<std::thread> threads_;

    const scene_static &scene_static_;
    const light_static &light_static_;
    
    bool rad_is_new_;
    bool emit_is_new_;

    aligned_buffer<col3f_sse> emit_new_;
    aligned_buffer<col3f_sse> emit_;
    aligned_buffer<col3f_sse> rad_;
    aligned_buffer<col3f_sse> rad2_;
    
    
    
    uint64_t pints_;
    uint64_t pints_last_;
    cl_ubyte64 pints_last_time_;
    
    std::atomic<bool> start_flag_;
    std::atomic<bool> abort_flag_;
    
    join_threads joiner_;
};

// deactivated after change of light_static::f_target addressing scheme
#if 0
class rad_core_lockfree: public rad_core {
    //typedef std::mutex lock_type;
    typedef spinlock_mutex lock_type;
public:
    class worker {

    };

    
    std::vector<std::pair<size_t,size_t>> calc_plane_distribution( const size_t num_partition ) {
        size_t num_ints = 0;
        std::vector<std::pair<size_t,size_t>> parts;
        for( auto &ff : light_static_.f_fact() ) {
            num_ints += ff.size();
        }
        
        std::cout << "num_ints: " << num_ints << "\n";
        
        size_t ints_per_part = num_ints / num_partition;
        
        size_t first = 0;
        size_t last = 0;
        size_t acc = 0;
        for( size_t i = 0; i < num_partition; ++i ) {
            while( last < light_static_.f_fact().size() && acc < ints_per_part ) {
                acc += light_static_.f_fact()[last].size();
                ++last;
            }
            acc = 0;
            parts.emplace_back( first, last );
            first = last;
        }
        
        parts.back().second = light_static_.f_fact().size();
        
        
        for( auto &p : parts ) { 
            std::cout << "part: " << p.first << " " << p.second << "\n";
        }
        return parts;
    }

    rad_core_lockfree( const scene_static &scene_static, const light_static &light_static /*, const std::vector<std::pair<size_t,size_t>> &part_bounds */)
            : 
            scene_static_(scene_static),
            light_static_(light_static),
            
            emit_( light_static.num_planes() ), rad_( light_static.num_planes() ), rad2_( light_static.num_planes() ),
            //ffs_(ffs), ff_target_(ff_target),
            //planes_(planes),
            pints_(0),
            pints_last_(0),
            pints_last_time_(0)
    {

        
        
        if ( !true ) {
            threads_.push_back( std::thread( [&]() {
                work(0, scene_static_.planes().size(), 0);
            }) );

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
        } else {
#if 1
            const size_t num_planes = scene_static_.planes().size();
            const size_t num_threads = 1;
            
            auto part = calc_plane_distribution(num_threads);
            for ( size_t i = 0; i < num_threads; ++i ) {
                
                
                
                //const size_t first = num_planes / num_threads * i;
                //const size_t last = num_planes / num_threads * (i+1);
                const size_t first = part.at(i).first;
                const size_t last = part.at(i).second;
            
                std::cout << "thread: " << i << " " << first << " " << last << " " << num_planes << "\n";
                
                threads_.push_back( std::thread( [=]() {
                    work( first, last, i );
                }));
            }
#else
            
            const size_t num_threads = part_bounds.size();
            
            
            for ( size_t i = 0; i < num_threads; ++i ) {
                
                
                
                //const size_t first = num_planes / num_threads * i;
                //const size_t last = num_planes / num_threads * (i+1);
                const size_t first = part_bounds.at(i).first;
                const size_t last = part_bounds.at(i).second;
            
                std::cout << "thread: " << i << " " << first << " " << last <<  "\n";
                
                threads_.push_back( std::thread( [=]() {
                    work( first, last, i );
                }));
            }
#endif

//          thread0_ = std::thread( [&]() { work(0, planes_.size()/4);
//          });
//          thread1_ = std::thread( [&]() { work(planes_.size()/4, planes_.size()/4 * 2);
//          });
//          thread2_ = std::thread( [&]() { work(planes_.size()/4 * 2, planes_.size()/4 * 3);
//          });
//          thread3_ = std::thread( [&]() { work(planes_.size()/4 * 3, planes_.size()/4 * 4);
//          });
        }
    }

    virtual void set_emit( const std::vector<vec3f> &emit ) {
        

        //assert( emit_new_.size() == emit.size() );

        emit_.assign( std::begin(emit), std::end(emit));
        //emit_is_new_ = true;
    }

    virtual bool update() {
        //return rad_is_new_;
        
        return true;
//         std::lock_guard<std::mutex> lock(mtx_);
// 
//         if ( rad_is_new_ ) {
//             rad_is_new_ = false;
//             return true;
//         } else {
//             return false;
//         }


    }

    virtual void copy( std::vector<vec3f> *out ) {
        col3f_sse * rb_src = nullptr;
        
        
        bool dir = buf_dir_.load();
        
        if( !dir ) {
            rb_src = rad_.data();
         
        } else {
            rb_src = rad2_.data();
         
        }
        for ( size_t i = 0; i < rad_.size(); ++i ) {
            (*out)[i].r = rb_src[i].r;
            (*out)[i].g = rb_src[i].g;
            (*out)[i].b = rb_src[i].b;
        }

        buf_dir_.store(!dir);
        
        
        cl_ubyte64 time = CL_System::get_microseconds();
        
        cl_ubyte64 dt = time - pints_last_time_;
        
        if( dt >= 1e6 ) {
            const size_t pt = pints_.load();
            
            const size_t dp = pt - pints_last_;
            pints_last_ = pt;
            
            pints_last_time_ = time;
            
            std::cout << "pint/s: " << dp / (dt / 1.0e6) << "\n";
        }
    }




private:

    void work( size_t first, size_t last, size_t rank ) {
        while (true) {
            //std::copy( emit_new_.begin(), emit_new_.end(), emit_.begin() );
        

//             tick_timer tt;
            do_radiosity_sse(first, last);

//             std::cout << "elapsed: " << tt.elapsed() << "\n";

//             if( rank == 0 )
//             {
//                 std::lock_guard<lock_type> lock(mtx_);
//                 rad_is_new_ = true;
//             }

        }
    }

    void do_radiosity_sse( size_t first, size_t last ) {
        const int steps = 1;
        const float min_ff = 0;
//      std::fill(e_rad_sse_.begin(), e_rad_sse_.end(), vec3f(0.0, 0.0, 0.0));
//      std::fill(e_rad2_sse_.begin(), e_rad2_sse_.end(), vec3f(0.0, 0.0, 0.0));

        //std::copy( emit_sse_.begin(), emit_sse_.end(), e_rad_sse_.begin() );

        const std::vector<std::vector<float> > &ffs_ = light_static_.f_fact();
        const std::vector<std::vector<int> > &ff_target_ = light_static_.f_target();
        const std::vector<plane> &planes_ = scene_static_.planes();
        
        typedef vector_unit<float,4> vu;
        typedef vu::vec_t vec_t;
        //steps = 0;
//         std::cout << "do: " << first << " " << last << "\n";
        vec_t reflex_factor = vu::set1(1.0);
        
        size_t pints = 0;
        for ( int i = 0; i < steps; ++i ) {
            //for( auto it = pairs_.begin(); it != pairs_.end(); ++it, ++ff_it ) {


            for ( size_t j = first; j < last; ++j ) {

                col3f_sse *rb_src = nullptr;
                col3f_sse *rb_dst = nullptr;
                
                if( !buf_dir_.load() ) {
                    rb_src = rad_.data();
                    rb_dst = rad2_.data();
                } else {
                    rb_src = rad2_.data();
                    rb_dst = rad_.data();
                }
                    
                
                
                const size_t s = ffs_[j].size();
                vec_t rad = vu::set1(0);

                const vec3f cd = planes_[j].col_diff();
                const vec_t col_diff = vu::set( 0, cd.b, cd.g, cd.r );



                for ( size_t k = 0; k < s; ++k ) {
                    size_t target = ff_target_[j][k];
                    //rad_rgb += (col_diff * e_rad_rgb_[target]) * ff2s_[j][k];

                    if ( false && ffs_[j][k] < min_ff ) {
                        continue;
                    }

                    const vec_t ff = vu::set1( ffs_[j][k] );

                    rad = vu::add( rad, vu::mul( vu::mul( col_diff, vu::load( (float*) (rb_src+target))), ff ));

                }
                pints += s;
                vu::store( vu::add( vu::load((float*)emit_(j)), vu::mul(rad, reflex_factor)), (float*)(rb_dst+j));
//              std::cout << "col: " << rad2_[j].r << " " << cd.r <<  "\n";
                //e_rad2_rgb_[j] = emit_rgb_[j] + rad_rgb;// * reflex_factor;

            }

            pints_.fetch_add(pints);
            
            //rad_ = rad2__;


        }

        //e_rad_rgb_.assign( e_rad_sse_.begin(), e_rad_sse_.end() );



    }

    lock_type mtx_;
//  std::thread thread0_;
//  std::thread thread1_;
//  std::thread thread2_;
//  std::thread thread3_;

    std::vector<std::thread> threads_;

    const scene_static &scene_static_;
    const light_static &light_static_;
    
//     bool rad_is_new_;
//     bool emit_is_new_;

//     aligned_buffer<col3f_sse> emit_new_;
    aligned_buffer<col3f_sse> emit_;
    aligned_buffer<col3f_sse> rad_;
    aligned_buffer<col3f_sse> rad2_;
    
    std::atomic<bool> buf_dir_;
    
    std::atomic<size_t> pints_;
    size_t pints_last_;
    cl_ubyte64 pints_last_time_;
};
#endif

std::unique_ptr<rad_core> make_rad_core_threaded(const scene_static &scene_static, const light_static &light_static) {
    return make_unique<rad_core_threaded>( scene_static, light_static );
  //  return make_unique<rad_core_lockfree>( scene_static, light_static );
}

