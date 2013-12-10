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

	//set up dam-break test scene
	vec3 p0 = vec3(6.4f, 1.0f, 6.4f);
	vec3 p1 = vec3(12.8f, 12.8f, 25.6f);
	geomCore::box* fluid1 = new geomCore::box(p0, p1, FLUID);

	p0 = vec3(1.0f, 1.0f, 1.0f);
	p1 = vec3(32.0f-1.0f, 1.92f, 32.0f-1.0f);
	geomCore::box* fluid2 = new geomCore::box(p0, p1, FLUID);

	//set up walls
	p0 = vec3(0.0f, 0.0f, 0.0f); // left wall
	p1 = vec3(1.0f, dimensions.y, dimensions.z);
	geomCore::box* leftwall = new geomCore::box(p0, p1, SOLID);

	p0 = vec3(dimensions.x-1.0f, 0.0f, 0.0f); // right wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	geomCore::box* rightwall = new geomCore::box(p0, p1, SOLID);

	p0 = vec3(0.0f, 0.0f, 0.0f); // floor wall
	p1 = vec3(dimensions.x, 1.0f, dimensions.z);
	geomCore::box* floorwall = new geomCore::box(p0, p1, SOLID);

	p0 = vec3(0.0f, dimensions.y-1.0f, 0.0f); // ceiling wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	geomCore::box* ceilingwall = new geomCore::box(p0, p1, SOLID);

	p0 = vec3(0.0f, 0.0f, 0.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, 1.0f);
	geomCore::box* frontwall = new geomCore::box(p0, p1, SOLID);

	p0 = vec3(0.0f, 0.0f, dimensions.z-1.0f); // front wall
	p1 = vec3(dimensions.x, dimensions.y, dimensions.z);
	geomCore::box* backwall = new geomCore::box(p0, p1, SOLID);

	sceneCore::scene* dambreak = new sceneCore::scene();
	dambreak->addGeom(fluid1);
	dambreak->addGeom(fluid2);
	dambreak->addGeom(leftwall);
	dambreak->addGeom(rightwall);
	dambreak->addGeom(floorwall);
	dambreak->addGeom(ceilingwall);
	dambreak->addGeom(frontwall);
	dambreak->addGeom(backwall);

    fluidCore::particlegrid* test = new fluidCore::particlegrid(dimensions.x,
    															dimensions.y,
    															dimensions.z);

    fluidCore::flipsim* f = new fluidCore::flipsim(dimensions, dambreak, density);

    viewerCore::viewer* glview = new viewerCore::viewer();
    glview->load(f);
    glview->launch();

}
