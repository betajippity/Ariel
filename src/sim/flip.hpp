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
		void step();

		vector<particle*>* getParticles();
		vec3 getDimensions();
		sceneCore::scene* getScene();

	private:
		void computeDensity();
		void applyExternalForces();
		void project();

		bool isCellFluid(const int& x, const int& y, const int& z);

		vec3 dimensions;
		vector<particle*> particles;
		macgrid mgrid;
		particlegrid* pgrid;

		int subcell;
		float density;
		float max_density;

		sceneCore::scene* scene;

		int timestep;
		float stepsize;

};
}

#endif
