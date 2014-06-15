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

enum vbotype{QUADS, TRIANGLES, LINES, POINTS};

using namespace std;
using namespace glm;

namespace viewerCore {

struct vboData{
	GLuint vboID;
	GLuint cboID;
	int size;
	vbotype type;
	string key;
};

//Used just for tracking OpenGL viewport camera position/keeping in sync with rendercam
struct glCamera{
	vec2 mouseOld;
	vec3 rotate;
	vec3 translate;
	vec2 fov;
	float lookat;
	int currentKey;
	int currentMouseClick;
	float rotateSpeed;
	float zoomSpeed;
	float panSpeed;

	//Initializer
	glCamera(): mouseOld(vec2(0.0f,0.0f)), 
				rotate(vec3(0.0f,0.0f,0.0f)), 
				translate(vec3(0.0f,0.0f,0.0f)),
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
		void load(fluidCore::flipsim* sim, bool retina, vec2 resolution, vec3 camrotate, 
				  vec3 camtranslate, vec2 camfov, float camlookat);
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
                          int colorcount, vbotype type, string key);
		vboData createVBOFromObj(objCore::objContainer* o, vec4 color, string key);

		void saveFrame();

		//Interface callbacks
		static void errorCallback(int error, const char* description);		
		static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

		//Cached data
		bool loaded;
		bool runrender;

		vec2 resolution;
		GLFWwindow* window;
		vector<vboData> vbos;
		map<string, int> vbokeys;
		map<string, vec2> frameranges;

		glCamera cam;

		vector<fluidCore::particle*>* particles;

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
