// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: viewer.hpp
// GL viewer based on standard Takua viewer framework

#ifndef VIEWER_HPP
#define VIEWER_HPP

#define GLEW_STATIC

#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include <thread> 
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <map>
#include "../sim/flip.hpp"

namespace viewerCore {

struct vboData{
	GLuint vboID;
	GLuint cboID;
	int size;
	GLenum type;
	std::string key;
};

//Used just for tracking OpenGL viewport camera position/keeping in sync with rendercam
struct glCamera{
	glm::vec2 mouseOld;
	glm::vec3 rotate;
	glm::vec3 translate;
	glm::vec2 fov;
	float lookat;
	int currentKey;
	int currentMouseClick;
	float rotateSpeed;
	float zoomSpeed;
	float panSpeed;

	//Initializer
	glCamera(): mouseOld(glm::vec2(0.0f,0.0f)), 
				rotate(glm::vec3(0.0f,0.0f,0.0f)), 
				translate(glm::vec3(0.0f,0.0f,0.0f)),
				lookat(0.0f),
				fov(45.0f),
				currentKey(0),
				currentMouseClick(0),
				rotateSpeed(0.2f),
				zoomSpeed(0.1f),
				panSpeed(0.1f){};
};

//====================================
// Class Declarations
//====================================

class viewer{
	public:
		viewer();
		~viewer();

		bool launch();
		void load(fluidCore::flipsim* sim, bool retina);
		void load(fluidCore::flipsim* sim, bool retina, glm::vec2 resolution, glm::vec3 camrotate, 
				  glm::vec3 camtranslate, glm::vec2 camfov, float camlookat);
	private:
		//Initialize stuff
		bool init();

		//Sim thread stuff
		void simLoopThread();

		//Main draw functions
		void mainLoop();
		void updateInputs();

		//VBO stuff
		vboData createVBO(vboData data, float* vertices, int vertexcount, float* colors,
                          int colorcount, GLenum type, std::string key);
		vboData createVBOFromObj(objCore::objContainer* o, glm::vec4 color, std::string key);

		void saveFrame();

		//Interface callbacks
		static void errorCallback(int error, const char* description);		
		static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

		//Cached data
		bool loaded;
		bool runrender;

		glm::vec2 resolution;
		GLFWwindow* window;
		std::vector<vboData> vbos;
		std::map<std::string, int> vbokeys;
		std::map<std::string, glm::vec2> frameranges;

		glCamera cam;

		std::vector<fluidCore::particle*>* particles;

    	fluidCore::flipsim* sim;
    	bool siminitialized;
    	bool drawobjects;
    	bool drawInvalid;
    	
    	unsigned char* bitmapData;
    	bool dumpFramebuffer;
    	bool dumpReady;
    	bool pause;
    	tbb::mutex framebufferWriteLock;

    	int framebufferScale;

    	bool dumpVDB;
    	bool dumpOBJ;
    	bool dumpPARTIO;
};
}

#endif
