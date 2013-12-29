// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sceneloader.hpp
// Simple JSON loader for fluid sim setups

#ifndef SCENELOADER_HPP
#define SCENELOADER_HPP

#include <json/json.h>
#include "../utilities/utilities.h"
#include "scene.hpp"

using namespace std;
using namespace glm;

namespace sceneCore {
//====================================
// Class Declarations
//====================================

class sceneloader {
	public:
		sceneloader(string filename);
		~sceneloader();

		scene* getScene();
		float getDensity();
		vec3 getDimensions();

	private:
		void loadSettings(const Json::Value& jsonsettings);
		void loadSim(const Json::Value& jsonsettings);
		void loadBox(const Json::Value& jsoncube);
		void loadSphere(const Json::Value& jsonsphere);
		void loadObj(const Json::Value& jsonobj);

		scene* s;
		vec3 dimensions;
		float density;
		string relativePath;

};
}

#endif
