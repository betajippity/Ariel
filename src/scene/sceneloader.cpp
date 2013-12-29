// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sceneloader.cpp
// Implements sceneloader.hpp

#include "sceneloader.hpp"

using namespace sceneCore;

sceneloader::sceneloader(string filename){
	cout << "Loading scene from " << filename << "..." << endl;

	//grab relative path
	vector<string> pathTokens = utilityCore::tokenizeString(filename, "/");
    if(strcmp(filename.substr(0,1).c_str(), "/")==0){
        relativePath = "/";
    }else{
        relativePath = "";
    }
    for(int i=0; i<pathTokens.size()-1; i++){
        relativePath = relativePath + pathTokens[i] + "/";
    }

	s = new scene();

	//do json stuff
	string jsonInput = utilityCore::readFileAsString(filename);
	Json::Value root;
	Json::Reader reader;
	bool parsedSuccess = reader.parse(jsonInput, root, false);

	//If error, report failures and their locations in the document 
	if(!parsedSuccess){
		cout << "Error: Failed to parse JSON" << endl << reader.getFormatedErrorMessages() << endl;
	}else{
		loadSettings(root["settings"][0]);
		loadSim(root["sim"]);
	}

	cout << "Loaded scene from " << filename << ".\n" << endl;
}

sceneloader::~sceneloader(){

}

scene* sceneloader::getScene(){
	return s;
}

float sceneloader::getDensity(){
	return density;
}

vec3 sceneloader::getDimensions(){
	return dimensions;
}

void sceneloader::loadSettings(const Json::Value& jsonsettings){
	density = .5f;
	dimensions = vec3(32);

	if(jsonsettings.isMember("density")){
		density = jsonsettings["density"].asFloat();
	}
	if(jsonsettings.isMember("dim_x")){
		dimensions.x = jsonsettings["dim_x"].asInt();
	}
	if(jsonsettings.isMember("dim_y")){
		dimensions.y = jsonsettings["dim_y"].asInt();
	}
	if(jsonsettings.isMember("dim_z")){
		dimensions.z = jsonsettings["dim_z"].asInt();
	}
}

void sceneloader::loadSim(const Json::Value& jsonsim){
	int jsoncount = jsonsim.size();
	for(int i=0; i<jsoncount; i++){
		Json::Value object = jsonsim[i];
		if(strcmp(object["shape"].asString().c_str(), "box")==0){
			loadBox(object);
		}else if(strcmp(object["shape"].asString().c_str(), "sphere")==0){
			loadSphere(object);
		}else if(strcmp(object["shape"].asString().c_str(), "obj")==0){
			loadObj(object);
		}
	}
}

void sceneloader::loadBox(const Json::Value& jsoncube){
	geomCore::cube cubebuilder;
	vec3 point0;
	point0.x = jsoncube["point0_x"].asFloat();
	point0.y = jsoncube["point0_y"].asFloat();
	point0.z = jsoncube["point0_z"].asFloat();
	vec3 point1;
	point1.x = jsoncube["point1_x"].asFloat();
	point1.y = jsoncube["point1_y"].asFloat();
	point1.z = jsoncube["point1_z"].asFloat();
	if(strcmp(jsoncube["type"].asString().c_str(), "solid")==0){
		s->addSolidObject(cubebuilder.tesselate(point0, point1));
	}else if(strcmp(jsoncube["type"].asString().c_str(), "liquid")==0){
		s->addLiquidObject(cubebuilder.tesselate(point0, point1));
	}
}

void sceneloader::loadSphere(const Json::Value& jsoncube){
	geomCore::sphere spherebuilder;
	vec3 center;
	center.x = jsoncube["center_x"].asFloat();
	center.y = jsoncube["center_y"].asFloat();
	center.z = jsoncube["center_z"].asFloat();
	float radius;
	radius = jsoncube["radius"].asFloat();
	if(strcmp(jsoncube["type"].asString().c_str(), "solid")==0){
		s->addSolidObject(spherebuilder.tesselate(center, radius));
	}else if(strcmp(jsoncube["type"].asString().c_str(), "liquid")==0){
		s->addLiquidObject(spherebuilder.tesselate(center, radius));
	}
}

void sceneloader::loadObj(const Json::Value& jsonobj){
	string objpath = relativePath + jsonobj["file"].asString();
	objCore::objContainer* objloader = new objCore::objContainer(objpath);

	vec3 scale(1.0f,1.0f,1.0f);
	vec3 rotate(0.0f,0.0f,0.0f);
	vec3 translate(0.0f,0.0f,0.0f);
	if(jsonobj.isMember("rotat_x")){ rotate.x = jsonobj["rotat_x"].asFloat(); }
	if(jsonobj.isMember("rotat_y")){ rotate.y = jsonobj["rotat_y"].asFloat(); }
	if(jsonobj.isMember("rotat_z")){ rotate.z = jsonobj["rotat_z"].asFloat(); }
	if(jsonobj.isMember("scale_x")){ scale.x = jsonobj["scale_x"].asFloat(); }
	if(jsonobj.isMember("scale_y")){ scale.y = jsonobj["scale_y"].asFloat(); }
	if(jsonobj.isMember("scale_z")){ scale.z = jsonobj["scale_z"].asFloat(); }
	if(jsonobj.isMember("trans_x")){ translate.x = jsonobj["trans_x"].asFloat(); }
	if(jsonobj.isMember("trans_y")){ translate.y = jsonobj["trans_y"].asFloat(); }
	if(jsonobj.isMember("trans_z")){ translate.z = jsonobj["trans_z"].asFloat(); }

	mat4 transform = utilityCore::buildTransformationMatrix(translate, rotate, scale);
	objloader->bakeTransform(transform);

	if(strcmp(jsonobj["type"].asString().c_str(), "solid")==0){
		s->addSolidObject(objloader);
	}else if(strcmp(jsonobj["type"].asString().c_str(), "liquid")==0){
		s->addLiquidObject(objloader);
	}
}