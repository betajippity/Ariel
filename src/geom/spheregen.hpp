// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: spheregen.hpp. Adapted from Takua Render.
// Defines the unit sphere geometry class, inherits from the generic geom class

#ifndef SPHEREGEN_HPP
#define SPHEREGEN_HPP

#include "meshgen.hpp"

namespace geomCore {
//====================================
// Class Declarations
//====================================
	
class SphereGen: public MeshGen {
	public:
		//Initializers
		SphereGen();
		SphereGen(const unsigned int& subdivCount);
		~SphereGen();

		//Getters
		void Tesselate(objCore::Obj* o);
		void Tesselate(objCore::Obj* o, const glm::vec3& center, const float& radius);

		//Data
		int m_subdivs;

	private:
		glm::vec3 GetPointOnSphereByAngles(const float& angle1, const float& angle2);
};
}

#endif
