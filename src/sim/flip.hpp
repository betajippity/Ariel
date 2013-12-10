// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: flip.hpp
// A flip simulator!

#ifndef FLIP_HPP
#define FLIP_HPP

#include "../grid/macgrid.inl"
#include "../grid/particlegrid.hpp"
#include "../scene/scene.hpp"
#include <vector>
#include <omp.h>

using namespace std;
using namespace glm;

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class flipsim{
	public:
		flipsim(const vec3& maxres, sceneCore::scene* scene, const float& density);
		~flipsim();

		void init();

		vector<particle*>* getParticles();
		vec3 getDimensions();

	private:
		void computeDensity();

		vec3 dimensions;
		vector<particle*> particles;
		macgrid mgrid;
		particlegrid* pgrid;

		float density;
		float max_density;

		sceneCore::scene* scene;

};
}

#endif
