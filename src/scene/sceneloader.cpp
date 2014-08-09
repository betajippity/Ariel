// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: sceneloader.cpp
// Implements sceneloader.hpp

#include "sceneloader.hpp"

namespace sceneCore{

SceneLoader::SceneLoader(const std::string& filename){
	std::cout << "Loading scene from " << filename << "...\n" << std::endl;

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
		if(root.isMember("settings")){
			std::cout << "Loading settings..." << std::endl;
			LoadSettings(root["settings"][0]);
		}
		if(root.isMember("camera")){
			std::cout << "Loading camera..." << std::endl;
			LoadCamera(root["camera"][0]);
		}
		if(root.isMember("globalforces")){
			std::cout << "Loading global forces..." << std::endl;
			LoadGlobalForces(root["globalforces"]);
		}
		if(root.isMember("transforms")){
			std::cout << "Loading transforms..." << std::endl;
			unsigned int nodesInGroup = root["transforms"].size();
			m_s->m_geomTransforms.reserve(nodesInGroup);
			for(unsigned int j=0; j<nodesInGroup; j++){
				LoadGeomTransforms(root["transforms"][j]);
			}
		}
		if(root.isMember("meshfiles")){
			std::cout << "Loading meshfiles..." << std::endl;
			unsigned int nodesInGroup = root["meshfiles"].size();
			m_s->m_meshFiles.reserve(nodesInGroup);
			for(unsigned int j=0; j<nodesInGroup; j++){
				LoadMeshFiles(root["meshfiles"][j]);
			}
		}
		if(root.isMember("animatedmeshes")){
			std::cout << "Loading animatedmeshes..." << std::endl;
			unsigned int nodesInGroup = root["animatedmeshes"].size();
			//since we want sequences to be instanceable, the loading method for sequences
			//gets a little complicated. we have to construct interpolatedObj objects for each
			//frame, store in a vector, and then store those vectors for later referencing
			m_animMeshSequences.reserve(nodesInGroup);
			unsigned int interpObjCount = 0;
			for(unsigned int j=0; j<nodesInGroup; j++){
				interpObjCount = interpObjCount + root["animatedmeshes"][j].size() + 1;
			}
			m_s->m_animMeshes.reserve(interpObjCount);
			for(unsigned int j=0; j<nodesInGroup; j++){
				LoadAnimMeshSequences(root["animatedmeshes"][j]);
			}
		}
		if(root.isMember("geoms")){
			std::cout << "Loading geoms..." << std::endl;
			unsigned int nodesInGroup = root["geoms"].size();
			m_s->m_geoms.reserve(nodesInGroup);
			unsigned int meshContainerCount = 0;
			unsigned int animmeshContainerCount = 0;
			for(unsigned int j=0; j<nodesInGroup; j++){
				if(strcmp(root["geoms"][j]["type"].asString().c_str(), "mesh")==0){
					meshContainerCount++;
				}else if(strcmp(root["geoms"][j]["type"].asString().c_str(), "animated_mesh")==0){
					animmeshContainerCount++;
				}
			}
			m_s->m_meshContainers.reserve(meshContainerCount);
			m_s->m_animmeshContainers.reserve(animmeshContainerCount);
			for(unsigned int j=0; j<nodesInGroup; j++){
				LoadGeom(root["geoms"][j]);
			}
		}
		if(root.isMember("sim")){
			std::cout << "Loading sim..." << std::endl;
			unsigned int fluidNodeCount = 0;
			unsigned int solidNodeCount = 0;
			unsigned int simNodeCount = root["sim"].size();
			for(unsigned int j=0; j<simNodeCount; j++){
				if(strcmp(root["sim"][j]["type"].asString().c_str(), "liquid")==0){
					fluidNodeCount++;
				}else if(strcmp(root["sim"][j]["type"].asString().c_str(), "solid")==0){
					solidNodeCount++;
				}
			}
			m_s->m_solids.reserve(solidNodeCount);
			m_s->m_liquids.reserve(fluidNodeCount);
			for(unsigned int j=0; j<simNodeCount; j++){
				LoadSim(root["sim"][j]);
			}
		}	
	}

	m_s->SetPaths(m_imagePath, m_meshPath, m_vdbPath, m_partioPath);

	std::cout << std::endl;
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

void SceneLoader::LoadSim(const Json::Value& jsonsim){
	std::string id = jsonsim["geom"].asString();
	unsigned int geomID = m_linkNames["geom_"+id];
	geomCore::Geom* geomnode = &m_s->m_geoms[geomID];
	if(strcmp(jsonsim["type"].asString().c_str(), "liquid")==0){
		m_s->m_liquids.push_back(geomnode);
	}else if(strcmp(jsonsim["type"].asString().c_str(), "solid")==0){
		m_s->m_solids.push_back(geomnode);
	}
}

void SceneLoader::LoadGeom(const Json::Value& jsongeom){
	if(jsongeom.isMember("id")==false){
		std::cout << "Warning: Couldn't load geom node, missing ID. Skipping...\n" 
				  << std::endl;
	}else{
		std::string id = jsongeom["id"].asString();
		//Check if id already exists
		std::map<std::string, unsigned int>::iterator it = m_linkNames.find("geom_"+id);
		if(it!=m_linkNames.end()){
			std::cout << "Warning: geom node with ID \"" << id
					  << "\" already exists! Skipping...\n" << std::endl;
		}else{
            //first read common properties
            bool prePersist = jsongeom["pre_persist"].asBool();
            bool postPersist = jsongeom["post_persist"].asBool();
            unsigned int frameInterval = jsongeom["frame_interval"].asInt();
            unsigned int frameOffset = jsongeom["frame_offset"].asInt();
            unsigned int transformFrameCount = jsongeom["transform_frames"].size();
            //then take care of type specific stuff
			if(strcmp(jsongeom["type"].asString().c_str(), "mesh")==0){
                //check frame counts
                unsigned int frameCount = glm::min(jsongeom["geom_frames"].size(), 
                                                   transformFrameCount);
                //build transform pointers array
                geomCore::GeomTransform** transforms = new geomCore::GeomTransform*[frameCount];
                for(unsigned int i=0; i<frameCount; i++){
                    std::string transformname = jsongeom["transform_frames"][i].asString();
                    unsigned int transformID = m_linkNames["transform_"+transformname];
                    transforms[i] = &m_s->m_geomTransforms[transformID]; 
                } 
                //build meshFile pointers array
                spaceCore::Bvh<objCore::Obj>** meshes = new spaceCore::Bvh<objCore::Obj>*
                                                            [frameCount];
                for(unsigned int i=0; i<frameCount; i++){
                    std::string meshname = jsongeom["geom_frames"][i].asString();
                    unsigned int meshID = m_linkNames["meshfile_"+meshname];
                    meshes[i] = &m_s->m_meshFiles[meshID];
                } 
                m_s->m_meshContainers.push_back(geomCore::MeshContainer(frameCount, frameOffset,
                                                                        frameInterval, prePersist,
                                                                        postPersist, transforms,
                                                                        meshes)); 
			    unsigned int meshNodeID = m_s->m_meshContainers.size()-1;
                m_s->m_meshContainers[meshNodeID].m_id = meshNodeID;
                m_linkNames["meshinterface_"+id] = meshNodeID;
                m_s->m_geoms.push_back(geomCore::Geom(&m_s->m_meshContainers[meshNodeID]));
            }else if(strcmp(jsongeom["type"].asString().c_str(), "animated_mesh")==0){
                //check frame counts
                std::string animmeshname = jsongeom["anim_sequence"].asString();
                unsigned int animmeshID = m_linkNames["animmesh_"+animmeshname];
                unsigned int frameCount = m_animMeshSequences[animmeshID].size();
                frameCount = glm::min(frameCount, transformFrameCount);
                //build transform pointers array
                geomCore::GeomTransform** transforms = new geomCore::GeomTransform*[frameCount];
                for(unsigned int i=0; i<frameCount; i++){
                    std::string transformname = jsongeom["transform_frames"][i].asString();
                    unsigned int transformID = m_linkNames["transform_"+transformname];
                    transforms[i] = &m_s->m_geomTransforms[transformID]; 
                } 
                //build meshFile pointers array
                spaceCore::Bvh<objCore::InterpolatedObj>** animmeshes = 
                                        new spaceCore::Bvh<objCore::InterpolatedObj>*[frameCount];
                for(unsigned int i=0; i<frameCount; i++){
                    animmeshes[i] = m_animMeshSequences[animmeshID][i];    
                }
                m_s->m_animmeshContainers.push_back(geomCore::AnimatedMeshContainer(
                                                                frameCount, frameOffset,
                                                                frameInterval, prePersist,
                                                                postPersist, transforms,
                                                                animmeshes)); 
			    unsigned int animmeshNodeID = m_s->m_animmeshContainers.size()-1;
                m_s->m_animmeshContainers[animmeshNodeID].m_id = animmeshNodeID;
                m_linkNames["animmeshinterface_"+id] = animmeshNodeID;
                m_s->m_geoms.push_back(geomCore::Geom(&m_s->m_animmeshContainers[animmeshNodeID]));
            }
            unsigned int geomNodeID = m_s->m_geoms.size()-1;
            m_s->m_geoms[geomNodeID].m_id = geomNodeID;
	        m_linkNames["geom_"+id] = geomNodeID;
        }
    }
}

void SceneLoader::LoadAnimMeshSequences(const Json::Value& jsonanimmesh){
	if(jsonanimmesh.isMember("id")==false){
		std::cout << "Warning: Couldn't load animmesh node, missing ID. Skipping...\n" 
				  << std::endl;
	}else{
		std::string id = jsonanimmesh["id"].asString();
		//Check if id already exists
		std::map<std::string, unsigned int>::iterator it = m_linkNames.find("animmesh_"+id);
		if(it!=m_linkNames.end()){
			std::cout << "Warning: animmesh node with ID \"" << id
					  << "\" already exists! Skipping...\n" << std::endl;
		}else{
			m_animMeshSequences.push_back(std::vector<
                                          spaceCore::Bvh<objCore::InterpolatedObj>*>());
			unsigned int nodeNumber = 0;
            nodeNumber = (unsigned int)m_animMeshSequences.size()-1;
            unsigned int frameCount = jsonanimmesh["frames"].size();
            m_animMeshSequences[nodeNumber].reserve(frameCount);
            std::cout << "Creating animmesh with " << frameCount << " frames" << std::endl;
            for(unsigned int i=0; i<frameCount-1; i++){
            	//grab current and next frame IDs to pass to an InterpolatedObj
            	std::string thisframelink = jsonanimmesh["frames"][i].asString();
            	std::string nextframelink = jsonanimmesh["frames"][i+1].asString();
            	unsigned int thisframeID = m_linkNames["meshfile_"+thisframelink];
            	unsigned int nextframeID = m_linkNames["meshfile_"+nextframelink];
            	m_s->m_animMeshes.push_back(spaceCore::Bvh<objCore::InterpolatedObj>());
            	unsigned int animMeshNodeNumber = (unsigned int)m_s->m_animMeshes.size()-1;
            	objCore::InterpolatedObj interpObj(&m_s->m_meshFiles[thisframeID].m_basegeom,
            									   &m_s->m_meshFiles[nextframeID].m_basegeom);
            	m_s->m_animMeshes[animMeshNodeNumber] = interpObj;
            	m_s->m_animMeshes[animMeshNodeNumber].BuildBvh(24);
            	m_s->m_animMeshes[animMeshNodeNumber].m_id = animMeshNodeNumber;
            	m_animMeshSequences[nodeNumber].push_back(&m_s->m_animMeshes[animMeshNodeNumber]);
            }
            //last frame has to be a special case since there is next frame
            std::string thisframelink = jsonanimmesh["frames"][frameCount-1].asString();
            unsigned int thisframeID = m_linkNames["meshfile_"+thisframelink];
            m_s->m_animMeshes.push_back(spaceCore::Bvh<objCore::InterpolatedObj>());
            unsigned int animMeshNodeNumber = (unsigned int)m_s->m_animMeshes.size()-1;
        	objCore::InterpolatedObj interpObj(&m_s->m_meshFiles[thisframeID].m_basegeom,
        									   &m_s->m_meshFiles[thisframeID].m_basegeom);
        	m_s->m_animMeshes[animMeshNodeNumber] = interpObj;
        	m_s->m_animMeshes[animMeshNodeNumber].BuildBvh(24);
        	m_s->m_animMeshes[animMeshNodeNumber].m_id = animMeshNodeNumber;
        	m_animMeshSequences[nodeNumber].push_back(&m_s->m_animMeshes[animMeshNodeNumber]);
        	m_linkNames["animmesh_"+id] = nodeNumber;
		}	
	}
}

void SceneLoader::LoadGeomTransforms(const Json::Value& jsontransforms){
	if(jsontransforms.isMember("id")==false){
		std::cout << "Warning: Couldn't load transform node, missing ID. Skipping...\n" 
				  << std::endl;
	}else{
		std::string id = jsontransforms["id"].asString();
		//Check if id already exists
		std::map<std::string, unsigned int>::iterator it = m_linkNames.find("transform_"+id);
		if(it!=m_linkNames.end()){
			std::cout << "Warning: transform node with ID \"" << id
					  << "\" already exists! Skipping...\n" << std::endl;
		}else{
			glm::vec3 translation(0);
			glm::vec3 rotation(0);
			glm::vec3 scale(1);
			//read available properties into node
			if(jsontransforms.isMember("translation")==true){
				translation.x = jsontransforms["translation"]["x"].asFloat();
				translation.y = jsontransforms["translation"]["y"].asFloat();
				translation.z = jsontransforms["translation"]["z"].asFloat();
			}
			if(jsontransforms.isMember("rotation")==true){
				rotation.x = jsontransforms["rotation"]["x"].asFloat();
				rotation.y = jsontransforms["rotation"]["y"].asFloat();
				rotation.z = jsontransforms["rotation"]["z"].asFloat();
			}
			if(jsontransforms.isMember("scale")==true){
				scale.x = jsontransforms["scale"]["x"].asFloat();
				scale.y = jsontransforms["scale"]["y"].asFloat();
				scale.z = jsontransforms["scale"]["z"].asFloat();
			}
			//push back geomFrame and set up id links	
			m_s->m_geomTransforms.push_back(geomCore::GeomTransform(translation, rotation, scale));
			unsigned int nodeNumber = 0;
            nodeNumber = (unsigned int)m_s->m_geomTransforms.size()-1;
			m_s->m_geomTransforms[nodeNumber].m_id = nodeNumber;
			m_linkNames["transform_"+id] = nodeNumber;
		}
	}
}

void SceneLoader::LoadMeshFiles(const Json::Value& jsonmeshfile){
	if(jsonmeshfile.isMember("id")==false){
		std::cout << "Warning: Couldn't load meshfile node, missing ID. Skipping...\n" 
				  << std::endl;
	}else{
		std::string id = jsonmeshfile["id"].asString();
		//Check if id already exists
		std::map<std::string, unsigned int>::iterator it = m_linkNames.find("meshfile_"+id);
		if(it!=m_linkNames.end()){
			std::cout << "Warning: meshfile node with ID \"" << id
					  << "\" already exists! Skipping...\n" << std::endl;
		}else{
			m_s->m_meshFiles.push_back(spaceCore::Bvh<objCore::Obj>());
			unsigned int nodeNumber = (unsigned int)m_s->m_meshFiles.size()-1;
			//meshfile can either point to obj file or request a mesh generator
			if(jsonmeshfile.isMember("file")==true){	
				std::string filename = jsonmeshfile["file"].asString();
				m_s->m_meshFiles[nodeNumber].m_basegeom.ReadObj(m_relativePath+filename);
			}else if(jsonmeshfile.isMember("mesh_gen")==true){
				std::string gentype = jsonmeshfile["mesh_gen"].asString();
				if(strcmp(gentype.c_str(), "box")==0){
					std::cout << "Ran box mesh gen" << std::endl;
					geomCore::CubeGen cubebuilder;
					glm::vec3 point0;
					point0.x = jsonmeshfile["point0"]["x"].asFloat();
					point0.y = jsonmeshfile["point0"]["y"].asFloat();
					point0.z = jsonmeshfile["point0"]["z"].asFloat();
					glm::vec3 point1;
					point1.x = jsonmeshfile["point1"]["x"].asFloat();
					point1.y = jsonmeshfile["point1"]["y"].asFloat();
					point1.z = jsonmeshfile["point1"]["z"].asFloat();
					cubebuilder.Tesselate(&m_s->m_meshFiles[nodeNumber].m_basegeom, 
										  point0, point1);
				}else if(strcmp(gentype.c_str(), "sphere")==0){
					std::cout << "Ran sphere mesh gen" << std::endl;
					geomCore::SphereGen spherebuilder;
					glm::vec3 center;
					center.x = jsonmeshfile["center"]["x"].asFloat();
					center.y = jsonmeshfile["center"]["y"].asFloat();
					center.z = jsonmeshfile["center"]["z"].asFloat();
					float radius;
					radius = jsonmeshfile["radius"].asFloat();
					spherebuilder.Tesselate(&m_s->m_meshFiles[nodeNumber].m_basegeom, 
											center, radius);
				}
			}
			m_s->m_meshFiles[nodeNumber].BuildBvh(24);
			m_linkNames["meshfile_"+id] = nodeNumber;
			m_s->m_meshFiles[nodeNumber].m_basegeom.m_id = nodeNumber;
			m_s->m_meshFiles[nodeNumber].m_id = nodeNumber;
		}
	}
}

void SceneLoader::LoadGlobalForces(const Json::Value& jsonforces){
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
	if(jsonsettings.isMember("step_size")){
		m_stepsize = jsonsettings["step_size"].asFloat();
	}
	if(jsonsettings.isMember("dim")){
		m_dimensions.x = jsonsettings["dim"]["x"].asInt();
		m_dimensions.y = jsonsettings["dim"]["y"].asInt();
		m_dimensions.z = jsonsettings["dim"]["z"].asInt();
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
}
