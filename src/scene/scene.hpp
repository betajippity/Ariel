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
#include "../geom/geom.inl"
#include "../grid/particlegrid.hpp"
#include "../grid/levelset.hpp"

namespace sceneCore {
//====================================
// Class Declarations
//====================================

class scene {
	public:
		scene();
		~scene();

		void addSolidObject(objCore::objContainer* object, int startFrame, int endFrame);
		void addLiquidObject(objCore::objContainer* object, int startFrame, int endFrame);
		void generateParticles(std::vector<fluidCore::particle*>& particles, 
							   const glm::vec3& dimensions, const float& density, 
							   fluidCore::particlegrid* pgrid, const int& frame);

		std::vector<objCore::objContainer*>& getSolidObjects();
		std::vector<objCore::objContainer*>& getLiquidObjects();

		fluidCore::levelset* getSolidLevelSet();
		fluidCore::levelset* getLiquidLevelSet();

		void buildLevelSets(const int& frame);
		void setPaths(const std::string& imagePath, const std::string& meshPath, 
					  const std::string& vdbPath, const std::string& partioPath);

		glm::vec2 getSolidFrameRange(const int& index);
		glm::vec2 getLiquidFrameRange(const int& index);

		void projectPointsToSolidSurface(std::vector<glm::vec3>& points);

		void exportParticles(std::vector<fluidCore::particle*> particles, const float& maxd, 
							 const int& frame, const bool& VDB, const bool& OBJ, 
							 const bool& PARTIO);

		std::string imagePath;
		std::string meshPath;
		std::string vdbPath;
		std::string partioPath;

	private:
		void addParticle(const glm::vec3& pos, const geomtype& type, const float& thickness, 
						 const float& scale, std::vector<fluidCore::particle*>& particles, 
						 const int& frame);

		std::vector<objCore::objContainer*> solidObjects;
		std::vector<objCore::objContainer*> liquidObjects;

		fluidCore::levelset* solidLevelSet;
		fluidCore::levelset* liquidLevelSet;

		fluidCore::levelset* permaSolidLevelSet;
		fluidCore::levelset* permaLiquidLevelSet;

		bool permaLiquidSDFActive;
		bool permaSolidSDFActive;

		std::vector<glm::vec2> solidObjectFrameRanges;
		std::vector<glm::vec2> liquidObjectFrameRanges;
};
}

#endif
