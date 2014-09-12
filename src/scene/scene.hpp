// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.hpp
// Scene class, meant for storing the scene state and checking particle positions

#ifndef SCENE_HPP
#define SCENE_HPP

#include <vector>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include "../utilities/utilities.h"
#include "../grid/macgrid.inl"
#include "../geom/mesh.hpp"
#include "../grid/particlegrid.hpp"
#include "../grid/levelset.hpp"
#include "../spatial/bvh.hpp"

namespace sceneCore {
//====================================
// Class Declarations
//====================================

class Scene {
	friend class SceneLoader;
	public:
		Scene();
		~Scene();

		// void AddSolidObject(objCore::Obj* object, const int& startFrame, const int& endFrame);
		void GenerateParticles(std::vector<fluidCore::Particle*>& particles, 
							   const glm::vec3& dimensions, const float& density, 
							   fluidCore::ParticleGrid* pgrid, const int& frame);

		void AddExternalForce(glm::vec3 force);
		std::vector<glm::vec3>& GetExternalForces();

		fluidCore::LevelSet* GetSolidLevelSet();
		fluidCore::LevelSet* GetLiquidLevelSet();

		void BuildLevelSets(const int& frame);
		void SetPaths(const std::string& imagePath, const std::string& meshPath, 
					  const std::string& vdbPath, const std::string& partioPath);

		void ExportParticles(std::vector<fluidCore::Particle*> particles, 
							 const float& maxd, const int& frame, const bool& VDB, 
							 const bool& OBJ, const bool& PARTIO);

		std::vector<geomCore::Geom*>& GetSolidGeoms();
		std::vector<geomCore::Geom*>& GetLiquidGeoms();

		rayCore::Intersection IntersectSolidGeoms(const rayCore::Ray& r);
		bool CheckPointInsideSolidGeom(const glm::vec3& p, const float& frame, 
									   unsigned int& solidGeomID);
		bool CheckPointInsideLiquidGeom(const glm::vec3& p, const float& frame, 
									    unsigned int& liquidGeomID);
		bool CheckPointInsideGeomByID(const glm::vec3& p, const float& frame, 
									  const unsigned int& geomID);

		unsigned int GetLiquidParticleCount();

		std::string						                            m_imagePath;
		std::string						                            m_meshPath;
		std::string						                            m_vdbPath;
		std::string						                            m_partioPath;

		tbb::mutex                                                  m_particleLock;

	private:
		void AddLiquidParticle(const glm::vec3& pos, const glm::vec3& vel, const float& thickness, 
							   const float& scale, const int& frame, 
							   const unsigned int& liquidGeomID);
		void AddSolidParticle(const glm::vec3& pos, const float& thickness, const float& scale, 
							  const int& frame, const unsigned int& solidGeomID);

		fluidCore::LevelSet*			                            m_solidLevelSet;
		fluidCore::LevelSet*			                            m_liquidLevelSet;
		std::vector<glm::vec3>			                            m_externalForces;

		std::vector<geomCore::GeomTransform>						m_geomTransforms;
		std::vector<spaceCore::Bvh<objCore::Obj> >					m_meshFiles;
		std::vector<spaceCore::Bvh<objCore::InterpolatedObj> >		m_animMeshes;
		std::vector<geomCore::Geom>									m_geoms;
		std::vector<geomCore::MeshContainer>						m_meshContainers;
		std::vector<geomCore::AnimatedMeshContainer>				m_animmeshContainers;
		std::vector<geomCore::Geom*>								m_solids;
		std::vector<geomCore::Geom*>								m_liquids;	
		std::vector<glm::vec3>										m_liquidStartingVelocities;

		tbb::concurrent_vector<fluidCore::Particle*>				m_liquidParticles;
		tbb::concurrent_vector<fluidCore::Particle*>				m_permaSolidParticles;
		tbb::concurrent_vector<fluidCore::Particle*>				m_solidParticles;
	
		unsigned int 												m_liquidParticleCount;

};
}

#endif

