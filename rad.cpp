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


#include <ClanLib/application.h>
#include <ClanLib/core.h>
#include <ClanLib/sound.h>
#include <ClanLib/gui.h>
#include <ClanLib/display.h>
#include <ClanLib/swrender.h>
#include <ClanLib/gl.h>

#define BOOST_UBLAS_NDEBUG

#include <boost/numeric/ublas/matrix.hpp>
#include <algorithm>

using namespace boost::numeric::ublas;
class timer {
    double m_start;
public:
    timer() : m_start( CL_System::get_time() / 1000.0 ) {
        
    }
    
    double elapsed() {
     
        return (CL_System::get_time() / 1000.0) - m_start;
    }
    
    
};



template<typename T>
struct vec2_gen {
    T v[2];

    vec2_gen() {}

    vec2_gen( T x_, T y_ ) {
        v[0] = x_;
        v[1] = y_;
    }
    
    inline T &operator[]( size_t i ) {
        return v[i];
    }
    inline const T &operator[]( size_t i ) const {
        return v[i];
    }
    T &x() {
        return v[0];
    }
    T &y() {
        return v[1];
    }
    
    const T &x() const {
        return v[0];
    }
    const T &y() const {
        return v[1];
    }
    vec2_gen<T> norm() const {
        T l = sqrt(*this ^ *this);
        return vec2_gen<T>( v[0] / l, v[1] / l );
    }

    const vec2_gen<T> &operator+=( const vec2_gen<T> &other ) {
		v[0] += other.v[0];
		v[1] += other.v[1];


		return *this;
	}

};
template<typename T>
vec2_gen<T> operator+(const vec2_gen<T> &v1, const vec2_gen<T> &v2 ) {
    return vec2_gen<T>(v1.x() + v2.x(), v1.y() + v2.y());
}

template<typename T>
vec2_gen<T> operator-(const vec2_gen<T> &v1, const vec2_gen<T> &v2 ) {
    return vec2_gen<T>(v1.x() - v2.x(), v1.y() - v2.y());
}

template<typename T>
vec2_gen<T> operator*( T v1, const vec2_gen<T> &v2 ) {
    return vec2_gen<T>(v1 * v2.x(), v1 * v2.y());
}

template<typename T>
vec2_gen<T> operator*(const vec2_gen<T> &v1, T v2 ) {
    return vec2_gen<T>(v1.x() * v2, v1.y() * v2 );
}




template<typename T>
T operator^(const vec2_gen<T> &v1, const vec2_gen<T> &v2 ) {
    return v1.x() * v2.x() + v1.y() * v2.y();
}

template<typename T>
inline std::ostream &operator<<( std::ostream &os, const vec2_gen<T> &v ) {
    os << "[" << v.x() << " " << v.y() << "]";
    return os;
}


template<typename T>
struct vec3_gen {
    T v[3];

    vec3_gen() {}
    vec3_gen( T x_, T y_, T z_ ) {
        v[0] = x_;
        v[1] = y_;
        v[2] = z_;
    }
    
    vec3_gen( const vec2_gen<T> &other2 ) {
        v[0] = other2.x();
        v[1] = other2.y();
        v[2] = 0;
    }
    
    inline T &operator[]( size_t i ) {
        return v[i];
    }
    
    inline const T &operator[]( size_t i ) const {
        return v[i];
    }
    
    T &x() {
        return v[0];
    }
    T &y() {
        return v[1];
    }
    T &z() {
        return v[2];
    }
    
    const T &x() const {
        return v[0];
    }
    const T &y() const {
        return v[1];
    }
    const T &z() const {
        return v[2];
    }
    
    vec3_gen<T> norm() const {
    	T l = sqrt(*this ^ *this);
		return vec3_gen<T>( v[0] / l, v[1] / l, v[2] / l );
	}

    vec3_gen<T> operator-() {
    	return vec3_gen<T>( -v[0], -v[1], -v[2] );
    }


    const vec3_gen<T> &operator+=( const vec3_gen<T> &other ) {
    	v[0] += other.v[0];
    	v[1] += other.v[1];
    	v[2] += other.v[2];

    	return *this;
    }
};

template<typename T>
vec3_gen<T> operator*(const vec3_gen<T> &v1, const vec3_gen<T> &v2 ) {
    return vec3_gen<T>(v1.x() * v2.x(), v1.y() * v2.y(), v1.z() * v2.z());
}

template<typename T>
vec3_gen<T> operator*(const vec3_gen<T> &v1, T v ) {
    return vec3_gen<T>(v1.x() * v, v1.y() * v, v1.z() * v );
}

template<typename T>
vec3_gen<T> operator*( T v, const vec3_gen<T> &v2 ) {
    return vec3_gen<T>(v * v2.x(), v * v2.y(), v * v2.z());
}



template<typename T>
vec3_gen<T> operator+(const vec3_gen<T> &v1, const vec3_gen<T> &v2 ) {
    return vec3_gen<T>(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

template<typename T>
vec3_gen<T> operator-(const vec3_gen<T> &v1, const vec3_gen<T> &v2 ) {
    return vec3_gen<T>(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}



template<typename T>
T operator^(const vec3_gen<T> &v1, const vec3_gen<T> &v2 ) {
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}
template<typename T>
inline std::ostream &operator<<( std::ostream &os, const vec3_gen<T> &v ) {
    os << "[" << v.x() << " " << v.y() << " " << v.z() << "]";
    return os;
}

typedef vec3_gen<float> vec3f;
typedef vec2_gen<float> vec2f;

typedef vec2_gen<double> vec2d;
typedef vec3_gen<double> vec3d;

class patch {
    
public:
    enum dir_type {
        dir_flat,
        dir_up,
        dir_down,
        dir_left,
        dir_right
    };
    
private:
    vec2f pos_;
    
    dir_type dir_;
    
    vec3f norm_;
    
    void setup_norm() {
        switch( dir_ ) {
        case dir_flat:
            norm_.x() = 0;
            norm_.y() = 0;
            norm_.z() = 1;
            break;
        case dir_up:
            norm_.x() = 0;
            norm_.y() = -1;
            norm_.z() = 0;
            break;
        case dir_down:
            norm_.x() = 0;
            norm_.y() = 1;
            norm_.z() = 0;
            break;      
        case dir_left:
            norm_.x() = -1;
            norm_.y() = 0;
            norm_.z() = 0;
            break;           
        case dir_right:
            norm_.x() = 1;
            norm_.y() = 0;
            norm_.z() = 0;
            break;            
            
                
        }
    }
public:
    patch( int x, int y, dir_type dx ) : pos_(x,y), dir_(dx)
    {
        setup_norm();
    }
    inline int get_x() const { return pos_.x(); }
    inline int get_y() const { return pos_.y(); }
    
    inline const vec3f & norm() const { return norm_; }
    inline const vec2f & pos() const { return pos_; }

    inline vec2f trace_pos() const {
    	return vec2f( pos_.x() + norm_.x(), pos_.y() + norm_.y() );
    }

    inline dir_type dir() const {
    	return dir_;
    }
};

class particle {

public:
	particle() : active_(false) {}

	particle( vec2f pos, vec2f vel ) : pos_(pos), vel_(vel), active_(true) {

	}


	void timestep( float dt ) {
		if( !active_ ) {
			return;
		}

		vel_[1] += 9.81 * dt * 0.5;
		pos_ += vel_ * dt;
	}


	const vec2f &pos() const {
		return pos_;
	}

	bool active() const {
		return active_;
	}
private:
	vec2f pos_;
	vec2f vel_;

	bool active_;

};

float randf() {
	return float(std::rand()) / RAND_MAX;
}

class particle_system {
public:
	particle_system( size_t size, float rate, const vec2f &pos, const vec2f &vel )
		: parts_(size, particle(vec2f(0.0, 0.0), vec2f(0.0,0.0)) ),
		  pos_(pos), vel_(vel), ptr_(0), dcum_(0)
	{

	}

	void timestep( float dt ) {
		//size_t change = rate_ * dt;
		size_t change = 1;

		dcum_ += dt;

		if( dcum_ > 2.0 ) {
			for( size_t i = 0; i < change; ++i ) {
				if( ptr_ >= parts_.size() ) {
					ptr_ = 0;
				}

				parts_.at(ptr_) = particle(pos_, vel_ * float(2.0 + 0.5 * randf()));
				++ptr_;
			}
			dcum_ = 0;
		}
		std::for_each( parts_.begin(), parts_.end(), [&]( particle & part ) {
			std::cout << "part: " << part.pos() << "\n";
			part.timestep(dt);
		});
	}


	void draw( CL_GraphicContext gc ) {



		const float screen_width = 60;
		const float screen_height = 40;

		const float scale_x = gc.get_width() / screen_width;
		const float scale_y = gc.get_height() / screen_height;

		for( const particle & part : parts_ ) {
			const vec2f &pos = part.pos();

			CL_Draw::point( gc, CL_Pointf(pos.x() * scale_x, pos.y() * scale_y ), CL_Colorf( 1.0f, 0.5f, 0.0f ));
		}
	}



	void active_parts( std::vector<vec2f> *out ) {
		for( const particle & part : parts_ ) {
			if( part.active() ) {
				out->push_back(part.pos());
			}
		}
	}

private:
	std::vector<particle> parts_;
	float rate_;
	vec2f pos_;
	vec2f vel_;

	size_t ptr_;
	float dcum_;
};

class light {
public:


	void pos( const vec2f & p ) {
		pos_ = p;
	}
	const vec2f &pos() const {
		return pos_;
	}

	const vec3f &color() const {
		return color_;
	}
private:
	vec2f pos_;
	vec3f color_;
};


class patch_scene {


public:

	patch_scene( size_t width, size_t height )
		: width_(width), height_(height)

	{

		solid_.resize( height_, width_ );
		for( auto it1 = solid_.begin1(); it1 != solid_.end1(); ++it1 ) {
			std::fill( it1.begin(), it1.end(), false );
		}
	}

	void solid( size_t x, size_t y, bool s ) {
		solid_( y, x ) = s;
	}

	void setup_patches() {
		flat_patches_.clear();
		up_patches_.clear();
		down_patches_.clear();
		left_patches_.clear();
		right_patches_.clear();

		assert( solid_.size1() == height_ );
		assert( solid_.size2() == width_ );

		for( size_t y = 0; y < height_; ++y ) {
			for( size_t x = 0; x < width_; ++x ) {
				if( solid_(y, x) ) {

					std::cout << "solid: " << x << " " << y << "\n";
					if( x > 0 ) {
						if( !solid_(y,x-1) ) {
							left_patches_.push_back(patch(x, y, patch::dir_left) );
							ortho_patches_.push_back( left_patches_.back() );
							patches_.push_back( ortho_patches_.back() );
						}
					}
					if( x < width_ -1 ) {
						if( !solid_(y, x+1) ) {
							right_patches_.push_back(patch(x, y, patch::dir_right) );
							ortho_patches_.push_back( right_patches_.back() );
							patches_.push_back( ortho_patches_.back() );
						}
					}
					if( y > 0 ) {
						if( !solid_(y-1,x) ) {
							up_patches_.push_back( patch(x, y, patch::dir_up) );
							ortho_patches_.push_back( up_patches_.back() );
							patches_.push_back( ortho_patches_.back() );
						}
					}
					if( y < height_ - 1 ) {
						if( !solid_(y+1,x) ) {
							down_patches_.push_back( patch(x, y, patch::dir_down) );
							ortho_patches_.push_back( down_patches_.back() );
							patches_.push_back( ortho_patches_.back() );
						}
					}


				} else {
					flat_patches_.push_back(patch(x, y, patch::dir_flat));
					patches_.push_back( flat_patches_.back() );
				}
			}
		}



	}


	void setup_formfactors() {

		for( size_t i = 0; i < patches_.size(); ++i ) {
			vec2f p1 = patches_.at(i).pos();



			vec2f trace_p1 = patches_[i].trace_pos();
			for( size_t j = 0; j < patches_.size(); ++j ) {

	//            		std::cerr << "i: " << i << " " << j << "\n";
				vec2f trace_p2 = patches_[j].trace_pos();

				vec2f p2 = patches_[j].pos();

				double d2 = dist_sqr( p1.x(), p1.y(), p2.x(), p2.y() );

				bool dist_cull = false;


				vec3f p1_3d = p1;
				if( patches_[i].dir() != patch::dir_flat ) {
					p1_3d[2] = 1.0;
				}

				vec3f norm1 = patches_[i].norm();
				float ff = 0;
				if( d2 < 1 || d2 > 30 ) {
					dist_cull = true;
				} else {


					vec3f p2_3d = p2;
					if( patches_[j].dir() != patch::dir_flat ) {
						p2_3d[2] = 1.0;
					}

					const vec3f &norm2 = patches_[j].norm();

					vec3f dn = (p1_3d - p2_3d).norm();
					//std::cout << p1_3d << " " << p2_3d << " " << dn << "\n";
					//.norm();
					float ff1 = std::max( 0.0f, norm1 ^ -dn);
					float ff2 = std::max( 0.0f, norm2 ^ dn);

					ff = ff1 * ff2;
					ff = std::max( 0.2f, ff );
	//						std::cout << "ff: " << ff << "\n";
					//dist_cull = ff < 0.01;
				}



				//dist_cull = false;
				if( !dist_cull && i != j && !occluded( trace_p1.x(), trace_p1.y(), trace_p2.x(), trace_p2.y(), solid_ ) ) {

					pairs.push_back(std::make_pair(i,j));

					ff = std::max( 0.2f, ff );
					ffs.push_back(ff / (3.1415 * d2));
				}
			}
		}

		e_emit.resize(patches_.size());
		e_rad.resize(patches_.size());

		e_emit_rgb.resize(patches_.size());


		col_diff.resize(patches_.size());

	}


	void light_particles( const std::vector<vec2f> &parts )
	{


		std::fill( e_emit.begin(), e_emit.end(), 0.0 );

		for( vec2f pos : parts ) {
		if( pos.x() > width_ || pos.y() > height_ ) {
			continue;
		}

		for( size_t i = 0; i < patches_.size(); ++i ) {
				const patch &p = patches_[i];



				const int light_x = pos.x();
				const int light_y = pos.y();

				const int x = p.get_x();
				const int y = p.get_y();
				if( p.dir() == patch::dir_flat ) {

					bool occ = occluded( x, y, light_x, light_y, solid_ );



					float d2 = dist_sqr( x, y, light_x, light_y );
					if( !occ && d2 > 0.00001 ) {
						e_emit[i] += 10.0/(3.14159 * d2);
					}


					col_diff[i] = vec3f(1, 1, 1);
					e_emit_rgb[i] = col_diff[i] * e_emit[i];

				} else {


					vec2f norm( p.norm().x(), p.norm().y() );
					bool occ = occluded( p.get_x() + norm.x(), p.get_y() + norm.y(), light_x, light_y, solid_ );


					if( !occ ) {


						vec2f d_pl = (pos - p.pos());

						float d2 = d_pl ^ d_pl;

						float dot = norm ^ d_pl.norm();

						if( dot > 0 ) {
		//							std::cout << dot << "\n";
							e_emit[i] += (10 * dot) / (3.14159 * d2);
						}
					}


					col_diff[i] = vec3f(1.0,0.5,0.0);

					e_emit_rgb[i] = col_diff[i] * e_emit[i];

				}


			}
		}

	}

	void do_radiosity() {
		float_vec e_rad2( e_rad.size(), 0.0);
		std::fill(e_rad.begin(), e_rad.end(), 0.0);

		e_rad_rgb.clear();
		e_rad_rgb.resize(patches_.size(), vec3f(0.0, 0.0, 0.0));
		std::vector<vec3f> e_rad2_rgb(patches_.size(), vec3f(0.0, 0.0, 0.0));
//			std::fill(e_rad_rg.begin(), e_rad.end(), 0.0);



#if 1


		std::cout << "pairs: " << pairs.size() << "\n";
		assert( !pairs.empty());

		for( int i = 0; i < 5; ++i ) {


			double rad = 0;
			vec3f rad_rgb(0,0,0);

			size_t prev = pairs.front().first;
			auto ff_it = ffs.begin();

			assert( ffs.size() == pairs.size() );
			for( auto it = pairs.begin(); it != pairs.end(); ++it, ++ff_it ) {


				size_t i = it->first;
				size_t j = it->second;

				if( i != prev ) {
					e_rad2[prev] = e_emit[prev] + 0.5 * rad;
					e_rad2_rgb[prev] = e_emit_rgb[prev] + 0.5f * rad_rgb;

					rad = 0;
					rad_rgb = vec3f(0,0,0);

					prev = i;
				}


// 				vec2f p1 = patches_[i].pos();
//					vec3f p1_3d = p1;
//					if( patches_[i].dir() != patch::dir_flat ) {
//						p1_3d[2] = 1.0;
//					}
//
//					vec3f norm1 = patches_[i].norm();


// 				const vec2f &p2 = patches_[j].pos();
// 				double d2 = dist_sqr( p1.x(), p1.y(), p2.x(), p2.y() );
//
//					if( d2 < 1 ) {
//						continue;
//					}
//
//
//
//
//					vec3f p2_3d = p2;
//					if( patches_[j].dir() != patch::dir_flat ) {
//						p2_3d[2] = 1.0;
//					}
//
//					const vec3f &norm2 = patches_[j].norm();
//
//					vec3f dn = (p1_3d - p2_3d).norm();
//					//std::cout << p1_3d << " " << p2_3d << " " << dn << "\n";
//					//.norm();
//					float ff1 = std::max( 0.0f, norm1 ^ -dn);
//					float ff2 = std::max( 0.0f, norm2 ^ dn);
//
//					float ff = ff1 * ff2;
////							if( ff > 0.0001 ) {
////								std::cout << "ff: " << ff << "\n";
////							}
//
////							if( ff > 0.00001 ) {
////								std::cout << "ff: " << ff1 << " " << ff2 << " " << dn << " " << norm1 << " " << norm2 << "\n";
////							}
//
//					//ff += 0.1;
//
//					ff = std::max(0.15f,ff);
				//rad += 1.0/(3.14159 * d2) * e_rad[j] * *ff_it;

				//rad_rgb += float(*ff_it/(3.14159 * d2)) * (e_rad_rgb[j] * col_diff[j]);
				rad_rgb += *ff_it * (e_rad_rgb[j] * col_diff[j]);

			}
			//e_rad2[prev] = e_emit[prev] + 0.5 * rad;
			e_rad2_rgb[prev] = e_emit_rgb[prev] + 0.5f * rad_rgb;

			e_rad.swap(e_rad2);
			e_rad_rgb.swap(e_rad2_rgb);
			//m = m2;
		}

#else
		for( int i = 0; i < 10; ++i ) {
			for( auto it1 = patches_.begin(); it1 != patches_.end(); ++it1 ) {
				const size_t i1 = std::distance(patches_.begin(), it1 );

				vec2f p1 = it1->pos();

				vec2f trace_p1 = it1->trace_pos();

				vec3f p1_3d = p1;
				if( it1->dir() != patch::dir_flat ) {
					p1_3d[2] = 1.0;
				}

				vec3f norm1 = it1->norm();

				double rad = 0;
				vec3f rad_rgb(0,0,0);
				for( auto it2 = patches_.begin(); it2 != patches_.end(); ++it2 ) {
			//		double ff = form_factor( *it1, *it2 );

							  //  std::cout << "dot: " << ff << "\n";
					const size_t i2 = std::distance(patches_.begin(), it2 );



					if( it1 == it2 ) {
						continue;
					}

					const vec2f &p2 = it2->pos();
					double d2 = dist_sqr( p1.x(), p1.y(), p2.x(), p2.y() );

					if( d2 < 1 ) {
						continue;
					}


					vec2f trace_p2 = it2->trace_pos();

					vec3f p2_3d = p2;
					if( it2->dir() != patch::dir_flat ) {
						p2_3d[2] = 1.0;
					}

					const vec3f &norm2 = it2->norm();

					vec3f dn = (p1_3d - p2_3d).norm();
					//std::cout << p1_3d << " " << p2_3d << " " << dn << "\n";
					//.norm();





					bool occ = occluded( trace_p1.x(), trace_p1.y(), trace_p2.x(), trace_p2.y(), solid_ );

					//assert( m2( p2[1], p2[0] ) == 0 );
//						if( m2( p2[1], p2[0] ) != 0 ) {
//							std::cout << "m2: " << m2( p2[1], p2[0] ) << " " << p2[1] << " " << p2[0] << "\n";
//						}

					if( !occ ) {
						float ff1 = std::max( 0.0f, norm1 ^ -dn);
						float ff2 = std::max( 0.0f, norm2 ^ dn);

						float ff = ff1 * ff2;
//							if( ff > 0.0001 ) {
//								std::cout << "ff: " << ff << "\n";
//							}

//							if( ff > 0.00001 ) {
//								std::cout << "ff: " << ff1 << " " << ff2 << " " << dn << " " << norm1 << " " << norm2 << "\n";
//							}

						ff += 0.1;
						rad += 1.0/(3.14159 * d2) * e_rad[i2] * ff;

						rad_rgb += float(ff/(3.14159 * d2)) * (e_rad_rgb[i2] * col_diff[i2]);

//							if( rad != rad ) {
//								std::cout << d2 << " " << e_rad[i2] << " " << ff << "\n";
//								throw std::runtime_error("");
//							}
						//assert( rad == rad );
					}

				}

//					std::cout << "rad: " << rad << "\n";
				e_rad2[i1] = e_emit[i1] + 0.5 * rad;
				e_rad2_rgb[i1] = e_emit_rgb[i1] + 0.5f * rad_rgb;
			}
			e_rad.swap(e_rad2);
			e_rad_rgb.swap(e_rad2_rgb);
			//m = m2;
		}
#endif
		          // e_rad = e_emit;


	}


	const std::vector<vec3f> &get_e_rad_rgb() const {
		return e_rad_rgb;
	}


	void draw_tiles( CL_GraphicContext gc ) {
		float tile_width = 4;
		float tile_height = 4;

		for( auto it = patches_.begin(); it != patches_.end(); ++it ) {
			const size_t i = std::distance(patches_.begin(), it);
			if( it->dir() == patch::dir_flat ) {
				const int x = it->get_x();
				const int y = it->get_y();

				float x1 = tile_width * x;
				float x2 = x1 + tile_width;

				float y1 = tile_height * y;
				float y2 = y1 + tile_height;



// 				float v = e_rad[i];
				auto v_rgb = e_rad_rgb[i];

				CL_Draw::fill(gc, x1, y1, x2, y2, CL_Colorf( v_rgb[0], v_rgb[1], v_rgb[2] ));
			} else {
				float x1 = it->get_x() * tile_width + (tile_width/2);
				float y1 = it->get_y() * tile_height + (tile_height/2);

				float x2 = x1 + it->norm().x() * (tile_width/2);
				float y2 = y1 + it->norm().y() * (tile_height/2);



// 				float v = e_rad[i];
				auto v_rgb = e_emit_rgb[i];

				CL_Draw::line(gc, x1, y1, x2, y2, CL_Colorf( v_rgb.x(), v_rgb.y(), v_rgb.z() ) );
			}
		}
	}

private:

	static bool occluded(int x0, int y0, int x1, int y1, const matrix<bool> &solid)
	{
//         int x0 = 1;
//         int x1 = 10;
//
//         int y0 = 1;
//         int y1 = 5;
		int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
		int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
		int err = dx+dy, e2; /* error value e_xy */

		while(true){  /* loop */
			//setPixel(x0,y0);
			//m(y0,x0) = 1.0;
			if( solid(y0,x0) ) {
				return true;
			}

			if (x0==x1 && y0==y1) break;
			e2 = 2*err;
			if (e2 >= dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
			if (e2 <= dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
		}

		return false;
	}

	static float dist_sqr( int x0, int y0, int x1, int y1 ) {
		int dx = abs( x0 - x1 );
		int dy = abs( y0 - y1 );

		return float( dx * dx + dy * dy );
	}

	typedef std::vector<patch> patch_vec;

	patch_vec flat_patches_;
	patch_vec up_patches_;
	patch_vec down_patches_;
	patch_vec left_patches_;
	patch_vec right_patches_;
	patch_vec ortho_patches_;
	patch_vec patches_;

	typedef std::vector<float> float_vec;




	matrix<bool> solid_;
	size_t width_;
	size_t height_;

	std::vector<std::pair<size_t,size_t> > pairs;
	std::vector<float> ffs;

	float_vec e_emit;//(patches_.size());
	float_vec e_rad;//(patches_.size());

	std::vector<vec3f> e_emit_rgb;//(patches_.size());


	std::vector<vec3f> col_diff;//(patches_.size());

	std::vector<vec3f> e_rad_rgb;//(patches_.size());

};

class rad {
    

    CL_SetupCore setup_core_;
    CL_SetupSound setup_sound_;
    
    

    CL_SetupDisplay display_;
    //CL_SetupSWRender swrender;
    
    CL_SetupGL setup_gl_;
    CL_DisplayWindow wnd_;
    
    
    patch_scene pc_;
    
    double form_factor( const patch& p1, const patch& p2) {
        vec2f pos1 = p1.pos();
        vec2f pos2 = p2.pos();
        
        vec3f norm1 = p1.norm();
        
        vec2f d = pos2 - pos1;
        
        vec3f dnorm = d.norm();
        
        //std::cout << "d: " << d << " " << dnorm << "\n";
        
        double dot = dnorm ^ norm1;
        
        return dot;
        
    }
public:
    
    rad() : pc_(60 * 4, 40 * 4) {}

    void setup() {
        CL_OpenGLWindowDescription desc;
        desc.set_size( CL_Size( 1024, 768 ), true );    
        
        wnd_ = CL_DisplayWindow(desc);
        
        
        
        for( size_t i = 0; i < 5; i++ ) {
            pc_.solid( i + 10, 10, true );

            pc_.solid( 8, 13 + i, true );
        }
        
        
        for( size_t i = 0; i < 40; ++i ) {
        	pc_.solid(0, i, true );
        	pc_.solid(i, 0, true );
        }

        pc_.setup_patches();
        pc_.setup_formfactors();
    }
    
    void mainloop() {
        CL_GraphicContext gc = wnd_.get_gc();
        
        
        
        

//        for( auto it = m.begin1(); it != m.end1(); ++it ) {
//            std::fill( it.begin(), it.end(), 0.0 );
//        }
//
//
//         int light_x = 2;
//         int light_y = 10;
//         int ld = 1;







		particle_system ps( 20, 1, vec2f(0.0f, 0.0f), vec2f( 1.0, 0.0 ));

		float dt = 0.01;

        while ( true ) {
            timer t1;
            


            
#if 0

            for( size_t i = 0; i < patches_.size(); ++i ) {
            	const patch &p = patches_[i];

            	const int x = p.get_x();
				const int y = p.get_y();
            	if( p.dir() == patch::dir_flat ) {

    				bool occ = occluded( x, y, light_x, light_y, solid_ );



    				float d2 = dist_sqr( x, y, light_x, light_y );
    				if( !occ && d2 > 0.00001 ) {
    					e_emit[i] = 100.0/(3.14159 * d2);
    				} else {
    					e_emit[i] = 0.0;
    				}


    				col_diff[i] = vec3f(1, 1, 1);
    				e_emit_rgb[i] = col_diff[i] * e_emit[i];

            	} else {


            		vec2f norm( p.norm().x(), p.norm().y() );
            		bool occ = occluded( p.get_x() + norm.x(), p.get_y() + norm.y(), light_x, light_y, solid_ );


            		if( !occ ) {


						vec2f d_pl = (light - p.pos());

						float d2 = d_pl ^ d_pl;

						float dot = norm ^ d_pl.norm();

						if( dot > 0 ) {
//							std::cout << dot << "\n";
							e_emit[i] = (100 * dot) / (3.14159 * d2);
						} else {
							e_emit[i] = 0;
						}
            		} else {
            			e_emit[i] = 0;
            		}

            		if( x == 0 || y == 0 ) {
            			//e_emit[i] = std::max( e_emit[i], 1.0f );
            			col_diff[i] = vec3f(0,0,1);
            		} else {
            			col_diff[i] = vec3f(1,1,0);
            		}

            		e_emit_rgb[i] = col_diff[i] * e_emit[i];

            	}


            }
#endif

            std::vector<vec2f> parts;
			ps.active_parts( &parts );

			pc_.light_particles(parts);

			pc_.do_radiosity();

#if 0
            for( int y = 0; y < 20; ++y ) {
                for( int x = 0; x < 30; ++x ) {
                    float x1 = 32 * x;
                    float x2 = x1 + 32;
                    
                    float y1 = 32 * y;
                    float y2 = y1 + 32;
                    
                    //float v = (y / 20.0) * 0.5 + (x / 30.0) * 0.5;
                    float v = m(y,x);
                    CL_Draw::fill(gc, x1, y1, x2, y2, CL_Colorf(v, v, v));        
//                     CL_Draw::box(gc, x1, y1, x2, y2, CL_Colorf::red);        
                }
            }
            std::cout << "t2: " << t1.elapsed() << "\n";

#endif
            

#if 1
            pc_.draw_tiles( gc );
            
            

#endif

            ps.timestep(dt);
            ps.draw( gc );

            std::cout << "t1: " << t1.elapsed() << "\n";
            wnd_.flip();
            
            dt = t1.elapsed();
//            CL_System::sleep( 30 );
            CL_KeepAlive::process();
        }
    }
    
    static int main(const std::vector<CL_String> &args) {
        
        
        rad r;
        r.setup();
        r.mainloop();
        
        
        return 0;
    }
    

};
CL_ClanApplication app(&rad::main);
