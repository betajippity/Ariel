// Ariel: FLIP Fluid Simulator
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
		flipsim(const vec3& maxres, sceneCore::scene* scene, const float& density, const gridtype& type,
				const bool& verbose);
		~flipsim();

		void init();
		void step(bool saveVDB, bool saveOBJ);

		vector<particle*>* getParticles();
		vec3 getDimensions();
		sceneCore::scene* getScene();

		int frame;

	private:
		void computeDensity();
		void applyExternalForces();
		void subtractPreviousGrid();
		void storePreviousGrid();
		void subtractPressureGradient();
		void extrapolateVelocity();
		void project();
		void solvePicFlip();
		void advectParticles();
		bool isCellFluid(const int& x, const int& y, const int& z);

		vec3 dimensions;
		vector<particle*> particles;
		macgrid mgrid;
		macgrid mgrid_previous;
		particlegrid* pgrid;

		int subcell;
		float density;
		float max_density;
		float densitythreshold;
		float picflipratio;

		sceneCore::scene* scene;

		float stepsize;

		gridtype type;
		bool verbose;
};
}

#endif
