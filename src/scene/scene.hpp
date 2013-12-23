// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: scene.hpp
// Scene class, meant for storing the scene state and checking particle positions

#ifndef SCENE_HPP
#define SCENE_HPP

#include "../utilities/utilities.h"
#include "../grid/macgrid.inl"
#include "../geom/geom.inl"
#include "../grid/particlegrid.hpp"
#include "../grid/levelset.hpp"
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

		void addSolidObject(objCore::objContainer* object);
		void addLiquidObject(objCore::objContainer* object);
		void generateParticles(vector<fluidCore::particle*>& particles, const vec3& dimensions, 
							   const float& density, fluidCore::particlegrid* pgrid);

		vector<objCore::objContainer*>& getSolidObjects();
		vector<objCore::objContainer*>& getLiquidObjects();

		fluidCore::levelset* getSolidLevelSet();
		fluidCore::levelset* getLiquidLevelSet();

		void rebuildLiquidLevelSet(vector<fluidCore::particle*>& particles);

	private:

		void addParticle(const vec3& pos, const geomtype& type, const float& thickness, const float& scale,
						 vector<fluidCore::particle*>& particles);

		vector<objCore::objContainer*> solidObjects;
		vector<objCore::objContainer*> liquidObjects;

		fluidCore::levelset* solidLevelSet;
		fluidCore::levelset* liquidLevelSet;

};
}

#endif
