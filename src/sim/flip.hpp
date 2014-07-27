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

class FlipSim{
	public:
		FlipSim(const glm::vec3& maxres, const float& density, const float& stepsize, 
				sceneCore::Scene* scene, const bool& verbose);
		~FlipSim();

		void Init();
		void Step(bool saveVDB, bool saveOBJ, bool savePARTIO);

		std::vector<Particle*>* GetParticles();
		glm::vec3 GetDimensions();
		sceneCore::Scene* GetScene();

		int							m_frame;

	private:
		void ComputeDensity();
		void ApplyExternalForces();
		void SubtractPreviousGrid();
		void StorePreviousGrid();
		void SubtractPressureGradient();
		void ExtrapolateVelocity();
		void Project();
		void SolvePicFlip();
		void AdvectParticles();
		bool IsCellFluid(const int& x, const int& y, const int& z);

		glm::vec3					m_dimensions;
		std::vector<Particle*>		m_particles;
		MacGrid						m_mgrid;
		MacGrid						m_mgrid_previous;
		ParticleGrid*				m_pgrid;

		int							m_subcell;
		float						m_density;
		float						m_max_density;
		float						m_densitythreshold;
		float						m_picflipratio;

		sceneCore::Scene*			m_scene;

		bool						m_verbose;
		float						m_stepsize;
};

class FlipTask: public tbb::task {
	public:
		FlipTask(FlipSim* sim, bool dumpVDB, bool dumpOBJ, bool dumpPARTIO);

		tbb::task* execute();
	private:
		FlipSim*					m_sim;
		bool						m_dumpPARTIO;
		bool						m_dumpOBJ;
		bool						m_dumpVDB;
};
}

#endif
