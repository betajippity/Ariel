// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: main.cpp
// Entry point for Ariel

#include <iostream>
#include "grid/particlegrid.hpp"
#include "sim/flip.hpp"
#include "viewer/viewer.hpp"

using namespace std;
using namespace glm;

int main(int argc, char** argv){ 

	cout << "" << endl;
	cout << "===================================================" << endl;
	cout << "Ariel: FLIP Fluid Simulator" << endl;
	cout << "Version 0.1.13.52a" << endl;
	cout << "Copyright (C) Yining Karl Li. All rights reserved." << endl;
	cout << "===================================================" << endl;
	cout << "" << endl;

	bool retina = false;
	bool verbose = false;

	for(int i=1; i<argc; i++){
		string header; string data;
		istringstream liness(argv[i]);
		getline(liness, header, '='); getline(liness, data, '=');
		if(strcmp(header.c_str(), "-retina")==0){
			retina = true;
			cout << "Framebuffer dumps set to Retina mode\n" << endl;
		}else if(strcmp(header.c_str(), "-verbose")==0){
			verbose = true;
			cout << "Verbose mode activated\n" << endl;
		}
	}

	vec3 dimensions = vec3(32,32,32);
	float density = .5f;

	geomCore::cube cubebuilder;
	sceneCore::scene* scene = new sceneCore::scene();

	//set up dam-break test scene
	vec3 p0 = vec3(6.4f, 1.0f, 6.4f);
	vec3 p1 = vec3(12.8f, 12.8f, 25.6f);
	scene->addLiquidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(1.0f, 1.0f, 1.0f);
	p1 = vec3(32.0f-1.0f, 1.92f, 32.0f-1.0f);
	scene->addLiquidObject(cubebuilder.tesselate(p0, p1));

	// //set up walls
	p0 = vec3(0.0f, 0.0f, 0.0f); // left wall
	p1 = vec3(1.0f, dimensions.y, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(dimensions.x-1.0f, 0.0f, 0.0f); // right wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, 0.0f, 0.0f); // floor wall
	p1 = vec3(dimensions.x, 1.0f, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, dimensions.y-1.0f, 0.0f); // ceiling wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, 0.0f, 0.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, 1.0f);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, 0.0f, dimensions.z-1.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	geomCore::sphere spherebuilder;
	// scene->addSolidObject(spherebuilder.tesselate(dimensions/2.0f, 10.0f));

    fluidCore::flipsim* f = new fluidCore::flipsim(dimensions, scene, density, RAW, verbose);

    viewerCore::viewer* glview = new viewerCore::viewer();
    glview->load(f, retina);
    glview->launch();

}
