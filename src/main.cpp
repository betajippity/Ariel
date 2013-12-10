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
	float density = 0.5f;

	//set up dam-break test scene
	vec3 p0 = vec3(6.4f, 1.0f, 6.4f);
	vec3 p1 = vec3(12.8f, 12.8f, 25.6f);
	geomCore::box* fluid1 = new geomCore::box(p0, p1, FLUID);

	p0 = vec3(1.0f, 1.0f, 1.0f);
	p1 = vec3(32.0f-1.0f, 1.92f, 32.0f-1.0f);
	geomCore::box* fluid2 = new geomCore::box(p0, p1, FLUID);

	sceneCore::scene* dambreak = new sceneCore::scene();
	dambreak->addGeom(fluid1);
	dambreak->addGeom(fluid2);


    fluidCore::particlegrid* test = new fluidCore::particlegrid(dimensions.x,
    															dimensions.y,
    															dimensions.z);

    fluidCore::flipsim* f = new fluidCore::flipsim(dimensions, dambreak, density);

    viewerCore::viewer* glview = new viewerCore::viewer();
    glview->load(f);
    glview->launch();

}
