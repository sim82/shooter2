#ifndef __player_bits_h
#define __player_bits_h

#include <ClanLib/core.h>
#include <ClanLib/display.h>
#include "misc_utils.h"

class input_mapper {
public:

    void add_mapping( int keycode, bool *indicator ) ;

    void input( const CL_InputDevice &dev ) ;

private:

    struct mapping {
        int keycode_;
        bool *indicator_;

        mapping( int keycode, bool *indicator ) : keycode_(keycode), indicator_(indicator) {
            *indicator_ = false;
        }

    };

    std::vector<mapping> mappings_;

};

class mouse_mapper {
public:
    mouse_mapper() : valid_(false), delta_x_(0), delta_y_(0) {}

    void input( const CL_InputDevice &mouse ) ;

    float delta_x() const {return delta_x_; }

    float delta_y() const {return delta_y_; }

private:

    bool valid_;
    int old_x_;
    int old_y_;

    float delta_x_;
    float delta_y_;
};

class player {
public:

    player() ;

    void input( const CL_InputDevice &keyboard, const CL_InputDevice &mouse ) ;

    void frame( double time, double dt ) ;
    const vec3f &pos() const { return pos_; }
    float rot_x() const { return rot_x_; }
    float rot_y() const { return rot_y_; }

private:
    bool i_forward_;
    bool i_backward_;
    bool i_left_;
    bool i_right_;

    input_mapper input_mapper_;
    mouse_mapper mouse_mapper_;

    vec3f pos_;
    float rot_y_;
    float rot_x_;

};


#endif