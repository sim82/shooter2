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

#include <thread>
#include <future>
#include <ClanLib/core.h>
#include <ClanLib/Core/XML/xml_token.h>
#include <ClanLib/Core/XML/xml_writer.h>
#include <ClanLib/Core/XML/xml_tokenizer.h>

#include "player_bits.h"

player::player() {
    input_mapper_.add_mapping( CL_KEY_W, &i_forward_);
    input_mapper_.add_mapping( CL_KEY_S, &i_backward_);
    input_mapper_.add_mapping( CL_KEY_A, &i_left_);
    input_mapper_.add_mapping( CL_KEY_D, &i_right_);
    input_mapper_.add_mapping( CL_KEY_SHIFT, &i_run_);
    
    input_mapper_.write_mappings( "test.xml" );
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

    
    
    const float base_speed = 4; // 4 m/s
    const float run_multiplier = 4;
    const float move_speed = i_run_ ? base_speed * run_multiplier : base_speed;

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
void input_mapper::write_mappings(const char* filename) {

    CL_DomDocument document;
    CL_DomElement root = document.create_element("my-structure");
    document.append_child(root);
    root.set_child_string("name", "xxx");
    root.set_child_string("scale", CL_StringHelp::float_to_text(666));
    root.set_child_bool("enabled", false);
    CL_File file(filename, CL_File::create_always, CL_File::access_read_write);
    document.save(file, false);
    
    auto x = std::async( [] { return std::string( "bla bla bla" ); } );
    
    std::cout << "bla: " << x.get() << "\n";
    
}
