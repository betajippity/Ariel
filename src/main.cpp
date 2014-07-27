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
	cout << "Version 0.1.14.24b" << endl;
	cout << "Copyright (C) Yining Karl Li. All rights reserved." << endl;
	cout << "===================================================" << endl;
	cout << "" << endl;

	bool retina = false;
	bool verbose = false;
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
		}
	}

	if(strcmp(scenefile.c_str(), "")==0){
		cout << "Error: no scene specified! Use -scene=[file]\n" << endl;
		exit(EXIT_FAILURE);
	} 

	sceneCore::SceneLoader* sloader = new sceneCore::SceneLoader(scenefile);

    fluidCore::FlipSim* f = new fluidCore::FlipSim(sloader->GetDimensions(), sloader->GetDensity(), 
    											   sloader->GetStepsize(), sloader->GetScene(), 
    											   verbose);

    viewerCore::Viewer* glview = new viewerCore::Viewer();
    glview->Load(f, retina, sloader->m_cameraResolution, sloader->m_cameraRotate, 
    			 sloader->m_cameraTranslate, sloader->m_cameraFov, sloader->m_cameraLookat);
    glview->Launch();

}
