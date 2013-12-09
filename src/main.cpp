// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: main.cpp
// Entry point for Kai

#include <iostream>
#include "grid/particlegrid.hpp"
#include "sim/flip.hpp"

using namespace std;
using namespace glm;

int main(int argc, char** argv){ 

    fluidCore::particlegrid* test = new fluidCore::particlegrid(32, 32, 32);

    fluidCore::flipsim* f = new fluidCore::flipsim(vec3(32), 0.5f);

    //lol i don't do anything yet...

}
