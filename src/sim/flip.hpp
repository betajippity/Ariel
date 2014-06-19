// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: flip.hpp
// A flip simulator!

#ifndef FLIP_HPP
#define FLIP_HPP

#include <tbb/tbb.h>
#include "../grid/macgrid.inl"
#include "../grid/particlegrid.hpp"
#include "../scene/scene.hpp"

namespace fluidCore {
//====================================
// Class Declarations
//====================================

class flipsim{
	public:
		flipsim(const glm::vec3& maxres, sceneCore::Scene* scene, const float& density, 
				const bool& verbose);
		~flipsim();

		void init();
		void step(bool saveVDB, bool saveOBJ, bool savePARTIO);

		std::vector<particle*>* getParticles();
		glm::vec3 getDimensions();
		sceneCore::Scene* getScene();

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

		glm::vec3 dimensions;
		std::vector<particle*> particles;
		macgrid mgrid;
		macgrid mgrid_previous;
		ParticleGrid* pgrid;

		int subcell;
		float density;
		float max_density;
		float densitythreshold;
		float picflipratio;

		sceneCore::Scene* scene;

		float stepsize;

		bool verbose;
};

class fliptask: public tbb::task {
	public:
		fliptask(flipsim* sim, bool dumpVDB, bool dumpOBJ, bool dumpPARTIO);

		tbb::task* execute();
	private:
		flipsim* m_sim;
		bool m_dumpPARTIO;
		bool m_dumpOBJ;
		bool m_dumpVDB;
};
}

#endif
