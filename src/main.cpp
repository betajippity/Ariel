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

	int n = 32;
	vec3 dimensions = vec3(n);
	float density = 0.5f;
	float thickness = 1.0/n;

	float scale = (float)n/32.0f;

	//set up dam-break test scene
	vec3 p0 = vec3(0.2f/scale,
				   thickness/scale,
				   0.2f/scale);
	vec3 p1 = vec3(0.4f/scale,
				   0.4f/scale,
				   0.8f/scale);
	geomCore::box* fluid1 = new geomCore::box(p0, p1, FLUID);

	p0 = vec3(thickness/scale,
			  thickness/scale,
			  thickness/scale);
	p1 = vec3((1.0f-thickness)/scale,
			   0.06f/scale,
			  (1.0f-thickness)/scale);
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
