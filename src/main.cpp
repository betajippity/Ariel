// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: main.cpp
// Entry point for Ariel

#include <iostream>
#include "grid/particlegrid.hpp"
#include "sim/flip.hpp"
#include "viewer/viewer.hpp"
#include "scene/sceneloader.hpp"

using namespace std;
using namespace glm;

int main(int argc, char** argv){ 

	cout << "" << endl;
	cout << "===================================================" << endl;
	cout << "Ariel: FLIP Fluid Simulator" << endl;
	cout << "Version 0.1.14.24a" << endl;
	cout << "Copyright (C) Yining Karl Li. All rights reserved." << endl;
	cout << "===================================================" << endl;
	cout << "" << endl;

	bool retina = false;
	bool verbose = false;
	gridtype type = RAW;
	string scenefile = "";

	for(int i=1; i<argc; i++){
		string header; string data;
		istringstream liness(argv[i]);
		getline(liness, header, '='); getline(liness, data, '=');
		if(strcmp(header.c_str(), "-retina")==0){
			retina = true;
			cout << "Framebuffer dumps set to Retina mode..." << endl;
		}else if(strcmp(header.c_str(), "-verbose")==0){
			verbose = true;
			cout << "Verbose mode activated..." << endl;
		}else if(strcmp(header.c_str(), "-scene")==0){
			scenefile = data;
		}else if(strcmp(header.c_str(), "-mode")==0){
			if(strcmp(data.c_str(), "vdb")==0){
				type = VDB;
			}
		}
	}

	if(type==VDB){
		cout << "Running in VDB grid mode...\n" << endl;
	}else if(type==RAW){
		cout << "Running in RAW grid mode...\n" << endl;
	}

	if(strcmp(scenefile.c_str(), "")==0){
		cout << "Error: no scene specified! Use -scene=[file]\n" << endl;
		exit(EXIT_FAILURE);
	} 

	sceneCore::sceneloader* sloader = new sceneCore::sceneloader(scenefile);

    fluidCore::flipsim* f = new fluidCore::flipsim(sloader->getDimensions(), sloader->getScene(), 
    											   sloader->getDensity(), type, verbose);

    viewerCore::viewer* glview = new viewerCore::viewer();
    glview->load(f, retina, sloader->camera_resolution, sloader->camera_rotate, 
    			 sloader->camera_translate, sloader->camera_fov, sloader->camera_lookat);
    glview->launch();

}
