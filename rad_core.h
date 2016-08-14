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

#pragma once

#include "misc_utils.h"
#include <memory>
#include <vector>
class scene_static;
class light_static;

class IRadCore
{
public:
    virtual void set_emit(const std::vector<vec3f> &emit) = 0;
    //     virtual bool update() = 0;
    virtual void copy(std::vector<vec3f> *out) = 0;
    virtual ~IRadCore() {}
};

std::unique_ptr<IRadCore> makeRadCoreThreaded(const scene_static &scene_static, const light_static &light_static);

