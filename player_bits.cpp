
#include <ClanLib/core.h>
#include "player_bits.h"

player::player() {
    input_mapper_.add_mapping( CL_KEY_W, &i_forward_);
    input_mapper_.add_mapping( CL_KEY_S, &i_backward_);
    input_mapper_.add_mapping( CL_KEY_A, &i_left_);
    input_mapper_.add_mapping( CL_KEY_D, &i_right_);

    pos_.x = 0;
    pos_.y = 0;
    pos_.z = 0;

    rot_y_ = 0.0;
    rot_x_ = 0.0;
}
void player::input(const CL_InputDevice& keyboard, const CL_InputDevice& mouse) {
    input_mapper_.input(keyboard);

    mouse_mapper_.input(mouse);
}
void player::frame(double time, double dt) {


    // set up player head to world rotation from x/y rotations
    rot_y_ -= mouse_mapper_.delta_x(); // NOTE: rotations are CCW, so clockwise head rotation (to the right) means negative mouse delta
    rot_x_ -= mouse_mapper_.delta_y();

//          CL_Mat4f rot = CL_Mat4f::rotate(CL_Angle(rot_x_, cl_degrees), CL_Angle(rot_y_, cl_degrees), CL_Angle(), cl_YXZ );

    CL_Quaternionf rot_quat(CL_Angle(rot_x_, cl_degrees), CL_Angle(rot_y_, cl_degrees), CL_Angle(), cl_YXZ );



    CL_Vec4f trans_vec(0.0, 0.0, 0.0, 1.0);

    const float move_speed = 4.0; // 4 m/s

    if ( i_forward_ ) {
        trans_vec.z -= move_speed * dt; // NOTE: translation is in 'player head coordinate system' so forward is -z
    }
    if ( i_backward_ ) {
        trans_vec.z += move_speed * dt;
    }
    if ( i_left_ ) {
        trans_vec.x -= move_speed * dt;
    }
    if ( i_right_ ) {
        trans_vec.x += move_speed * dt;
    }

    // rotate translation vector from 'player head coordinates' into world coordinates

    trans_vec = rot_quat.rotate_vector(trans_vec);

    pos_.x += trans_vec.x;
    pos_.y += trans_vec.y;
    pos_.z += trans_vec.z;



}
void mouse_mapper::input(const CL_InputDevice& mouse) {
    int new_x = mouse.get_x();
    int new_y = mouse.get_y();

    if ( valid_ ) {
        delta_x_ = new_x - old_x_;
        delta_y_ = new_y - old_y_;
    }

    old_x_ = new_x;
    old_y_ = new_y;
    valid_ = true;
}
void input_mapper::add_mapping(int keycode, bool* indicator) {
    mappings_.push_back(mapping( keycode, indicator ));
}
void input_mapper::input(const CL_InputDevice& dev) {
for ( const mapping &m : mappings_ ) {
        *m.indicator_ = dev.get_keycode(m.keycode_);
    }
}
