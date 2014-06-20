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

class SceneLoader {
	public:
		SceneLoader(const std::string& filename);
		~SceneLoader();

		Scene* GetScene();
		float GetDensity();
		glm::vec3 GetDimensions();

		glm::vec3		m_cameraRotate;
		glm::vec3		m_cameraTranslate;
		float			m_cameraLookat;
		glm::vec2		m_cameraResolution;
		glm::vec2		m_cameraFov;
		
	private:
		void LoadSettings(const Json::Value& jsonsettings);
		void LoadSim(const Json::Value& jsonsettings);
		void LoadBox(const Json::Value& jsoncube);
		void LoadSphere(const Json::Value& jsonsphere);
		void LoadObj(const Json::Value& jsonobj);
		void LoadCamera(const Json::Value& jsoncamera);

		Scene*			m_s;
		glm::vec3		m_dimensions;
		float			m_density;
		std::string		m_relativePath;
		std::string		m_imagePath;
		std::string		m_meshPath;
		std::string		m_vdbPath;
		std::string		m_partioPath;
		
};
}

#endif
