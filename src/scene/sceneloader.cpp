// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sceneloader.cpp
// Implements sceneloader.hpp

#include "sceneloader.hpp"

namespace sceneCore{

SceneLoader::SceneLoader(const std::string& filename){
	std::cout << "Loading scene from " << filename << "..." << std::endl;

	m_cameraRotate = glm::vec3(0);
	m_cameraTranslate = glm::vec3(0);
	m_cameraLookat = 0.0f;
	m_cameraResolution = glm::vec2(1024);
	m_cameraFov = glm::vec2(45.0f);

	//grab relative path
	std::vector<std::string> pathTokens = utilityCore::tokenizeString(filename, "/");
    if(std::strcmp(filename.substr(0,1).c_str(), "/")==0){
        m_relativePath = "/";
    }else{
        m_relativePath = "";
    }
    for(int i=0; i<pathTokens.size()-1; i++){
        m_relativePath = m_relativePath + pathTokens[i] + "/";
    }

	m_s = new Scene();

	//do json stuff
	std::string jsonInput = utilityCore::readFileAsString(filename);
	Json::Value root;
	Json::Reader reader;
	bool parsedSuccess = reader.parse(jsonInput, root, false);

	//If error, report failures and their locations in the document 
	if(!parsedSuccess){
		std::cout << "Error: Failed to parse JSON" << std::endl 
				  << reader.getFormatedErrorMessages() << std::endl;
	}else{
		LoadSettings(root["settings"][0]);
		LoadSim(root["sim"]);
		if(root.isMember("camera")){
			LoadCamera(root["camera"][0]);
		}
		if(root.isMember("forces")){
			LoadForces(root["forces"]);
		}
	}

	m_s->SetPaths(m_imagePath, m_meshPath, m_vdbPath, m_partioPath);

	std::cout << "Loaded scene from " << filename << ".\n" << std::endl;
}

SceneLoader::~SceneLoader(){

}

Scene* SceneLoader::GetScene(){
	return m_s;
}

float SceneLoader::GetDensity(){
	return m_density;
}

glm::vec3 SceneLoader::GetDimensions(){
	return m_dimensions;
}

float SceneLoader::GetStepsize(){
	return m_stepsize;
}

void SceneLoader::LoadForces(const Json::Value& jsonforces){
	unsigned int forceCount = jsonforces.size();
	for(unsigned int i=0; i<forceCount; i++){
		glm::vec3 force;
		force[0] = jsonforces[i]["x"].asFloat();
		force[1] = jsonforces[i]["y"].asFloat();
		force[2] = jsonforces[i]["z"].asFloat();
		m_s->AddExternalForce(force);
	}
}

void SceneLoader::LoadSettings(const Json::Value& jsonsettings){
	m_density = .5f;
	m_dimensions = glm::vec3(32);
	m_imagePath = m_relativePath;
	m_meshPath = m_relativePath;
	m_vdbPath = m_relativePath;
	m_stepsize = 0.005f;

	if(jsonsettings.isMember("density")){
		m_density = jsonsettings["density"].asFloat();
	}
	if(jsonsettings.isMember("stepsize")){
		m_stepsize = jsonsettings["stepsize"].asFloat();
	}
	if(jsonsettings.isMember("dim_x")){
		m_dimensions.x = jsonsettings["dim_x"].asInt();
	}
	if(jsonsettings.isMember("dim_y")){
		m_dimensions.y = jsonsettings["dim_y"].asInt();
	}
	if(jsonsettings.isMember("dim_z")){
		m_dimensions.z = jsonsettings["dim_z"].asInt();
	}
	if(jsonsettings.isMember("image_output")){
		m_imagePath = jsonsettings["image_output"].asString();
	}
	if(jsonsettings.isMember("mesh_output")){
		m_meshPath = jsonsettings["mesh_output"].asString();
	}
	if(jsonsettings.isMember("vdb_output")){
		m_vdbPath = jsonsettings["vdb_output"].asString();
	}
	if(jsonsettings.isMember("partio_output")){
		m_partioPath = jsonsettings["partio_output"].asString();
	}
}

void SceneLoader::LoadCamera(const Json::Value& jsoncamera){
	//load camera rotation
	if(jsoncamera.isMember("rotation")){
		m_cameraRotate[0] = jsoncamera["rotation"]["x"].asFloat();
		m_cameraRotate[1] = jsoncamera["rotation"]["y"].asFloat();
		m_cameraRotate[2] = jsoncamera["rotation"]["z"].asFloat();
	}
	//load camera translation
	if(jsoncamera.isMember("translation")){
		m_cameraTranslate[0] = jsoncamera["translation"]["x"].asFloat();
		m_cameraTranslate[1] = jsoncamera["translation"]["y"].asFloat();
		m_cameraTranslate[2] = jsoncamera["translation"]["z"].asFloat();
	}
	//load camera resolution
	if(jsoncamera.isMember("resolution")){
		m_cameraResolution[0] = jsoncamera["resolution"]["x"].asInt();
		m_cameraResolution[1] = jsoncamera["resolution"]["y"].asInt();	
	}
	//load camera lookat
	if(jsoncamera.isMember("lookat")){
		m_cameraLookat = jsoncamera["lookat"].asFloat();
	}
	//load camera fov
	if(jsoncamera.isMember("fovx")){
		m_cameraFov[0] = jsoncamera["fovx"].asFloat()/2.0f;
		float xscaled = tan(m_cameraFov.x*(PI/180));
	    float yscaled = (xscaled * m_cameraResolution.y)/m_cameraResolution.x;
	    m_cameraFov.y = (atan(xscaled)*180)/PI;
	}
}

void SceneLoader::LoadSim(const Json::Value& jsonsim){
	int jsoncount = jsonsim.size();
	for(int i=0; i<jsoncount; i++){
		Json::Value object = jsonsim[i];
		if(std::strcmp(object["shape"].asString().c_str(), "box")==0){
			LoadBox(object);
		}else if(std::strcmp(object["shape"].asString().c_str(), "sphere")==0){
			LoadSphere(object);
		}else if(std::strcmp(object["shape"].asString().c_str(), "obj")==0){
			LoadObj(object);
		}
	}
}

void SceneLoader::LoadBox(const Json::Value& jsoncube){
	geomCore::CubeGen cubebuilder;
	glm::vec3 point0;
	point0.x = jsoncube["point0_x"].asFloat();
	point0.y = jsoncube["point0_y"].asFloat();
	point0.z = jsoncube["point0_z"].asFloat();
	glm::vec3 point1;
	point1.x = jsoncube["point1_x"].asFloat();
	point1.y = jsoncube["point1_y"].asFloat();
	point1.z = jsoncube["point1_z"].asFloat();
	int startFrame = -1;
	int endFrame = -1;

	if(jsoncube.isMember("startFrame")){
		startFrame = jsoncube["startFrame"].asInt();
	}
	if(jsoncube.isMember("endFrame")){
		endFrame = jsoncube["endFrame"].asInt();
	}

	objCore::Obj* o = new objCore::Obj();
	cubebuilder.Tesselate(o, point0, point1);

	if(std::strcmp(jsoncube["type"].asString().c_str(), "solid")==0){
		m_s->AddSolidObject(o, startFrame, endFrame);
	}else if(std::strcmp(jsoncube["type"].asString().c_str(), "liquid")==0){
		m_s->AddLiquidObject(o, startFrame, endFrame);
	}
}

void SceneLoader::LoadSphere(const Json::Value& jsonsphere){
	geomCore::SphereGen spherebuilder;
	glm::vec3 center;
	center.x = jsonsphere["center_x"].asFloat();
	center.y = jsonsphere["center_y"].asFloat();
	center.z = jsonsphere["center_z"].asFloat();
	float radius;
	radius = jsonsphere["radius"].asFloat();
	int startFrame = -1;
	int endFrame = -1;

	if(jsonsphere.isMember("startFrame")){
		startFrame = jsonsphere["startFrame"].asInt();
	}
	if(jsonsphere.isMember("endFrame")){
		endFrame = jsonsphere["endFrame"].asInt();
	}

	objCore::Obj* o = new objCore::Obj();
	spherebuilder.Tesselate(o, center, radius);

	if(std::strcmp(jsonsphere["type"].asString().c_str(), "solid")==0){
		m_s->AddSolidObject(o, startFrame, endFrame);
	}else if(std::strcmp(jsonsphere["type"].asString().c_str(), "liquid")==0){
		m_s->AddLiquidObject(o, startFrame, endFrame);
	}
}

void SceneLoader::LoadObj(const Json::Value& jsonobj){
	std::string objpath = m_relativePath + jsonobj["file"].asString();
	objCore::Obj* obj = new objCore::Obj(objpath);

	glm::vec3 scale(1.0f,1.0f,1.0f);
	glm::vec3 rotate(0.0f,0.0f,0.0f);
	glm::vec3 translate(0.0f,0.0f,0.0f);
	if(jsonobj.isMember("rotat_x")){ rotate.x = jsonobj["rotat_x"].asFloat(); }
	if(jsonobj.isMember("rotat_y")){ rotate.y = jsonobj["rotat_y"].asFloat(); }
	if(jsonobj.isMember("rotat_z")){ rotate.z = jsonobj["rotat_z"].asFloat(); }
	if(jsonobj.isMember("scale_x")){ scale.x = jsonobj["scale_x"].asFloat(); }
	if(jsonobj.isMember("scale_y")){ scale.y = jsonobj["scale_y"].asFloat(); }
	if(jsonobj.isMember("scale_z")){ scale.z = jsonobj["scale_z"].asFloat(); }
	if(jsonobj.isMember("trans_x")){ translate.x = jsonobj["trans_x"].asFloat(); }
	if(jsonobj.isMember("trans_y")){ translate.y = jsonobj["trans_y"].asFloat(); }
	if(jsonobj.isMember("trans_z")){ translate.z = jsonobj["trans_z"].asFloat(); }

	glm::mat4 transform = utilityCore::buildTransformationMatrix(translate, rotate, scale);
	obj->BakeTransform(transform);

	int startFrame = -1;
	int endFrame = -1;

	if(jsonobj.isMember("startFrame")){
		startFrame = jsonobj["startFrame"].asInt();
	}
	if(jsonobj.isMember("endFrame")){
		endFrame = jsonobj["endFrame"].asInt();
	}

	if(std::strcmp(jsonobj["type"].asString().c_str(), "solid")==0){
		m_s->AddSolidObject(obj, startFrame, endFrame);
	}else if(std::strcmp(jsonobj["type"].asString().c_str(), "liquid")==0){
		m_s->AddLiquidObject(obj, startFrame, endFrame);
	}
}
}
