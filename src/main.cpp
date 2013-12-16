// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: main.cpp
// Entry point for Kai

#include <iostream>
#include "grid/particlegrid.hpp"
#include "sim/flip.hpp"
#include "viewer/viewer.hpp"

using namespace std;
using namespace glm;

int main(int argc, char** argv){ 

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
	// scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, 0.0f, 0.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, 1.0f);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	p0 = vec3(0.0f, 0.0f, dimensions.z-1.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	scene->addSolidObject(cubebuilder.tesselate(p0, p1));

	geomCore::sphere spherebuilder;
	// scene->addLiquidObject(spherebuilder.tesselate(dimensions/2.0f, 10.0f));

    fluidCore::particlegrid* test = new fluidCore::particlegrid(dimensions.x,
    															dimensions.y,
    															dimensions.z);

    fluidCore::flipsim* f = new fluidCore::flipsim(dimensions, scene, density);

    viewerCore::viewer* glview = new viewerCore::viewer();
    glview->load(f);
    glview->launch();

}
