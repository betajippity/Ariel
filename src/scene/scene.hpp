// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.hpp
// Scene class, meant for storing the scene state and checking particle positions

#ifndef SCENE_HPP
#define SCENE_HPP

#include "../utilities/utilities.h"
#include "../grid/macgrid.inl"
#include "../geom/geom.inl"
#include <vector>

using namespace std;
using namespace glm;

namespace sceneCore {
//====================================
// Class Declarations
//====================================

class scene {
	public:
		scene();
		~scene();

		void addGeom(geomCore::geom* object);
		void generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
							   const float& density);

	private:

		void addParticle(const vec3& pos, const geomtype& type, const float& thickness, 
						 vector<fluidCore::particle*>& particles);

		vector<geomCore::geom*> objects;

};
}

#endif
