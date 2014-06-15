// ObjCore2.6: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: objContainer.hpp
// A simple obj storage/loading class. This library superscedes ObjCore 1.0. 

#ifndef OBJCONTAINER_HPP
#define OBJCONTAINER_HPP

#include <glm/glm.hpp>
#include <sstream>
#include <fstream>
#include <iostream>
#include "obj.inl"
#include "../../utilities/utilities.h"

namespace objCore {
//====================================
// Class Declarations
//====================================

//Class wrapper for obj struct and related externed functions
class objContainer {
  public:
    //Initializers
    objContainer(std::string filename);
    objContainer(obj* o);
    ~objContainer();

    void bakeTransform(glm::mat4 transform);

   	//Getters
   	obj* getObj();
    bool writeObj(std::string filename); //returns true if successful, false otherwise
    
    void keepObj(bool keep); //if keep is true, destructor won't clear obj

  private:
    void preload(std::string filename);
    void load(std::string filename); 

    std::ifstream fp_in;
    obj* mesh;
    bool keep;
};
}

#endif
