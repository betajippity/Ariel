// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: geom.hpp. Adapted from Takua Render.
// Defines the generic geometry pure virtual class that all geometry will inherit from

#ifndef MESHGEN_HPP
#define MESHGEN_HPP

#include <glm/glm.hpp>
#include "obj/obj.hpp"

namespace geomCore {
//====================================
// Class Declarations
//====================================

class MeshGen{
    public:
        //Initializers
        MeshGen(){};
        virtual ~MeshGen(){};

        //Interactions
        virtual void Tesselate(objCore::Obj* o) = 0;
};
}

#endif
