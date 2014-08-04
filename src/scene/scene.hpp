// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.hpp
// Scene class, meant for storing the scene state and checking particle positions

#ifndef SCENE_HPP
#define SCENE_HPP

#include <vector>
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

		void AddSolidObject(objCore::Obj* object, const int& startFrame, const int& endFrame);
		void AddLiquidObject(objCore::Obj* object, const int& startFrame, const int& endFrame);
		void GenerateParticles(std::vector<fluidCore::Particle*>& particles, 
							   const glm::vec3& dimensions, const float& density, 
							   fluidCore::ParticleGrid* pgrid, const int& frame);

		void AddExternalForce(glm::vec3 force);
		std::vector<glm::vec3>& GetExternalForces();

		std::vector<objCore::Obj*>& GetSolidObjects();
		std::vector<objCore::Obj*>& GetLiquidObjects();

		fluidCore::LevelSet* GetSolidLevelSet();
		fluidCore::LevelSet* GetLiquidLevelSet();

		void BuildLevelSets(const int& frame);
		void SetPaths(const std::string& imagePath, const std::string& meshPath, 
					  const std::string& vdbPath, const std::string& partioPath);

		glm::vec2 GetSolidFrameRange(const int& index);
		glm::vec2 GetLiquidFrameRange(const int& index);

		void ProjectPointsToSolidSurface(std::vector<glm::vec3>& points);

		void ExportParticles(std::vector<fluidCore::Particle*> particles, const float& maxd, 
							 const int& frame, const bool& VDB, const bool& OBJ, 
							 const bool& PARTIO);

		//new stuff
		std::vector<geomCore::Geom*>& GetSolidGeoms();
		std::vector<geomCore::Geom*>& GetLiquidGeoms();

		std::string						m_imagePath;
		std::string						m_meshPath;
		std::string						m_vdbPath;
		std::string						m_partioPath;

	private:
		void AddParticle(const glm::vec3& pos, const geomtype& type, const float& thickness, 
						 const float& scale, std::vector<fluidCore::Particle*>& particles, 
						 const int& frame);

		std::vector< objCore::Obj* >		m_solidObjects;
		std::vector< objCore::Obj* >		m_liquidObjects;

		fluidCore::LevelSet*			m_solidLevelSet;
		fluidCore::LevelSet*			m_liquidLevelSet;

		fluidCore::LevelSet*			m_permaSolidLevelSet;
		fluidCore::LevelSet*			m_permaLiquidLevelSet;

		bool							m_permaLiquidSDFActive;
		bool							m_permaSolidSDFActive;

		std::vector<glm::vec2>			m_solidObjectFrameRanges;
		std::vector<glm::vec2>			m_liquidObjectFrameRanges;

		std::vector<glm::vec3>			m_externalForces;

		//new stuff
		std::vector<geomCore::GeomTransform>						m_geomTransforms;
		std::vector<spaceCore::Bvh<objCore::Obj> >					m_meshFiles;
		std::vector<spaceCore::Bvh<objCore::InterpolatedObj> >		m_animMeshes;
		std::vector<geomCore::Geom>									m_geoms;
		std::vector<geomCore::MeshContainer>						m_meshContainers;
		std::vector<geomCore::AnimatedMeshContainer>				m_animmeshContainers;
		std::vector<geomCore::Geom*>								m_solids;
		std::vector<geomCore::Geom*>								m_liquids;
};
}

#endif
