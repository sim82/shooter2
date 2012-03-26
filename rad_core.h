#ifndef __rad_core_h
#define __rad_core_h

#include <memory>
#include <vector>
#include "misc_utils.h"
class scene_static;
class light_static;

class rad_core {
public:
    virtual void set_emit( const std::vector<vec3f> &emit ) = 0;
//     virtual bool update() = 0;
    virtual void copy( std::vector<vec3f> *out ) = 0;
};


std::unique_ptr<rad_core> make_rad_core_threaded(const scene_static &scene_static, const light_static &light_static);

#endif