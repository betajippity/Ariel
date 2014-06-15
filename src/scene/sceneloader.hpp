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

namespace sceneCore {
//====================================
// Class Declarations
//====================================

class sceneloader {
	public:
		sceneloader(std::string filename);
		~sceneloader();

		scene* getScene();
		float getDensity();
		glm::vec3 getDimensions();

		glm::vec3 camera_rotate;
		glm::vec3 camera_translate;
		float camera_lookat;
		glm::vec2 camera_resolution;
		glm::vec2 camera_fov;
		
	private:
		void loadSettings(const Json::Value& jsonsettings);
		void loadSim(const Json::Value& jsonsettings);
		void loadBox(const Json::Value& jsoncube);
		void loadSphere(const Json::Value& jsonsphere);
		void loadObj(const Json::Value& jsonobj);
		void loadCamera(const Json::Value& jsoncamera);

		scene* s;
		glm::vec3 dimensions;
		float density;
		std::string relativePath;
		std::string imagePath;
		std::string meshPath;
		std::string vdbPath;
		std::string partioPath;
		
};
}

#endif
