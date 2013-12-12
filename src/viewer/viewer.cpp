// Kai: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: viewer.cpp
// Implements viewer.hpp

#include "viewer.hpp"
#include "../utilities/utilities.h"
#include <sstream>

using namespace viewerCore;

viewer::viewer(){
    loaded = false;
}

viewer::~viewer(){

}

void viewer::load(fluidCore::flipsim* sim){
    resolution = vec2(1000, 1000);

    cam.zoomSpeed = 1.0f;
    cam.panSpeed = 0.2f;

    loaded = true;

    frame = 0;
    // flip3D::init(frame);
    // particles = flip3D::getParticles();
    newframe = 0;

    recordWidth = 1000;
    recordHeight = 1000;
    bitmapData = new unsigned char[3 * recordWidth * recordHeight];

    this->sim = sim;
    siminitialized = false;
}

//returns true if viewer launches and closes successfully, otherwise, returns false
bool viewer::launch(){
    if(loaded==true){
        if(init()==true){


            omp_set_nested(true);
            #pragma omp parallel
            {
                #pragma omp master
                {
                    mainLoop();
                }
                #pragma omp single
                {
                    if(frame==0){
                        sim->init();
                        siminitialized = true;
                    }
                    while(1){
                        // flip3D::simulateStep();
                        particles = sim->getParticles();
                        frame++;
                    }
                }
            }

            return true;
        }else{
            cout << "Error: GL initialization failed.\n" << endl;
            return false;
        }
    }else{
        cout << "Error: No sim loaded!\n" << endl;
        return false;
    } 
}

//====================================
// Draw/Interaction Loop
//====================================

void viewer::updateInputs(){
    double x; double y;
    glfwGetCursorPos(window, &x, &y);
    vec2 d;
    d.x = float(x-cam.mouseOld.x);
    d.y = float(y-cam.mouseOld.y);
    cam.mouseOld.x = x;
    cam.mouseOld.y = y;
    if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == 1 || 
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == 1 ||
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == 1){

        bool doCamera = false;

        if(glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS || 
           glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS){
            doCamera = true;
        }
        if(doCamera==true){
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == 1){
                cam.rotate.x += d.y * cam.rotateSpeed;
                cam.rotate.y += d.x * cam.rotateSpeed;
            }
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == 1){
                cam.translate.z += d.y * cam.zoomSpeed;
            }
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == 1){
                cam.translate.x += d.x * cam.panSpeed;
                cam.translate.y -= d.y * cam.panSpeed;
            } 
        }else{
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == 1){
                //mouseclick event goes here
            }
        }
    }
}

void viewer::mainLoop(){
    while (!glfwWindowShouldClose(window)){

        if(siminitialized){
            vboData data;
            vector<float> vertexData;
            int psize = particles->size();

            vec3 gridSize = sim->getDimensions();
            float maxd = glm::max(glm::max(gridSize.x, gridSize.z), gridSize.y);

            for(int j=0; j<psize; j++){
                if(particles->operator[](j)->type==FLUID){
                    vertexData.push_back(particles->operator[](j)->p[0]*maxd);
                    vertexData.push_back(particles->operator[](j)->p[1]*maxd);
                    vertexData.push_back(particles->operator[](j)->p[2]*maxd);
                }
            }
            glDeleteBuffers(1, &data.vboID);
            data.color = vbos[vbokeys["fluid"]].color;
            string key = "fluid";
            data = createVBO(data, (float*)&vertexData[0], vertexData.size(), POINTS, key);
            vertexData.clear();
            vbos[vbokeys["fluid"]] = data;

            vboData data2;
            vector<float> vertexData2;
            for(int j=0; j<psize; j++){
                if(particles->operator[](j)->type==SOLID){
                    vertexData2.push_back(particles->operator[](j)->p[0]*maxd-(gridSize.x/2.0f));
                    vertexData2.push_back(particles->operator[](j)->p[1]*maxd-0.4f);
                    vertexData2.push_back(particles->operator[](j)->p[2]*maxd-(gridSize.z/2.0f));
                }
            }
            data2.color = vbos[vbokeys["walls"]].color;
            key = "walls";
            data2 = createVBO(data2, (float*)&vertexData2[0], vertexData2.size(), POINTS, key);
            vertexData2.clear();
            vbos[vbokeys["walls"]] = data2;
        }

        glClearColor(0.325, 0.325, 0.325, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glEnable (GL_BLEND);
        glBlendFunc (GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);

        glPushMatrix();
            glTranslatef(cam.translate.x, cam.translate.y, cam.translate.z);
            // apply the current rotation
            glRotatef(cam.rotate.x, 1, 0, 0);
            glRotatef(cam.rotate.y, 0, 1, 0);
            glRotatef(cam.rotate.z, 0, 0, 1);
            
            for(int i=0; i<vbos.size(); i++){
                glPushMatrix();
                glBindBuffer(GL_ARRAY_BUFFER, vbos[i].vboID);
                glVertexPointer(3, GL_FLOAT, 0, NULL);
                glEnableClientState(GL_VERTEX_ARRAY);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glColor4f(vbos[i].color.x, vbos[i].color.y, vbos[i].color.z, 0.5f);

                vec3 res = sim->getDimensions();
                glTranslatef(-res.x/2, 0, -res.y/2);

                if(i==vbokeys["boundingbox"]){
                    glLineWidth(2.0f);
                }

                if(vbos[i].type==QUADS){
                    glDrawArrays(GL_QUADS, 0, vbos[i].size/3);
                }else if(vbos[i].type==TRIANGLES){
                    glDrawArrays(GL_TRIANGLES, 0, vbos[i].size/3);
                }else if(vbos[i].type==LINES){
                    glDrawArrays(GL_LINES, 0, vbos[i].size/3);
                }else if(vbos[i].type==POINTS){
                    glPointSize(5.0f);
                    glDrawArrays(GL_POINTS, 0, vbos[i].size/3);
                }
                glDisableClientState(GL_VERTEX_ARRAY);
                glPopMatrix();
            }

            //draw unit axis
            glLineWidth(2.0f);
            glBegin(GL_LINES);
                glColor4f(1.0f, 0.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(2.0f, 0.0f, 0.0f);
                glColor4f(0.0f, 1.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 2.0f, 0.0f);
                glColor4f(0.0f, 0.0f, 1.0f, 0.0f);
                glVertex3f(0.0f, 0.0f, 0.0f);
                glVertex3f(0.0f, 0.0f, 2.0f);
            glEnd();
            glLineWidth(1.0f);

        glPopMatrix();

        glfwSwapBuffers(window);
        glfwPollEvents();
        updateInputs();

        // for (int i=0; i<recordHeight; i++) 
        // {
        //     glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
        //         bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
        // }

    }
    glfwDestroyWindow(window);
    glfwTerminate();
}

//====================================
// Init Stuff
//====================================

bool viewer::init(){
    //Camera setup stuff
    vec2 fov = vec2(45.0f, 45.0f);
    cam.translate = vec3(0.0f,0.0f,-30.0f);

    //Window setup stuff
    glfwSetErrorCallback(errorCallback);
    if (!glfwInit()){
        return false;
    }

    window = glfwCreateWindow(resolution.x, resolution.y, "Kai: now with 100% more VDB!", NULL, NULL);
    if (!window){
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, keyCallback);   
    glewExperimental = GL_TRUE;
    if(glewInit()!=GLEW_OK){
        return false;   
    }

    //camera stuff
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    vec2 xBounds;
    vec2 yBounds;
    utilityCore::fovToPerspective(fov.x, 1, 1, xBounds, yBounds); 
    glFrustum(xBounds[0], xBounds[1], yBounds[0], yBounds[1], 1, 10000000);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_MODELVIEW);

    //dummy buffer for particles
    vboData data;
    vector<float> vertexData;
    data.color = vec3(0,0,1);
    string key = "fluid";
    data = createVBO(data, (float*)&vertexData[0], vertexData.size(), POINTS, key);
    vertexData.clear();
    vbos.push_back(data);
    vbokeys["fluid"] = vbos.size()-1;

    data.color = vec3(1,0,0);
    key = "walls";
    data = createVBO(data, (float*)&vertexData[0], vertexData.size(), POINTS, key);
    vertexData.clear();
    vbos.push_back(data);
    vbokeys["walls"] = vbos.size()-1;

    //buffer for sim bounding box
    data.color = vec3(.2,.2,.2);
    key = "boundingbox";
    vector<vec3> boxData;
    vec3 res = sim->getDimensions();

    boxData.push_back(vec3(res.x, res.y, res.z));   boxData.push_back(vec3(0.0f, res.y, res.z));
    boxData.push_back(vec3(0.0f, res.y, res.z));    boxData.push_back(vec3(0.0f, 0.0f, res.z));
    boxData.push_back(vec3(0.0f, 0.0f, res.z));     boxData.push_back(vec3(res.x, 0.0f, res.z));
    boxData.push_back(vec3(res.x, 0.0f, res.z));    boxData.push_back(vec3(res.x, res.y, res.z));

    boxData.push_back(vec3(res.x, res.y, 0.0f));    boxData.push_back(vec3(0.0f, res.y, 0.0f));
    boxData.push_back(vec3(0.0f, res.y, 0.0f));     boxData.push_back(vec3(0.0f, 0.0f, 0.0f));
    boxData.push_back(vec3(0.0f, 0.0f, 0.0f));      boxData.push_back(vec3(res.x, 0.0f, 0.0f));
    boxData.push_back(vec3(res.x, 0.0f, 0.0f));     boxData.push_back(vec3(res.x, res.y, 0.0f));

    boxData.push_back(vec3(res.x, res.y, res.z));   boxData.push_back(vec3(res.x, res.y, 0.0f));
    boxData.push_back(vec3(0.0f, res.y, res.z));    boxData.push_back(vec3(0.0f, res.y, 0.0f));
    boxData.push_back(vec3(0.0f, 0.0f, res.z));     boxData.push_back(vec3(0.0f, 0.0f, 0.0f));
    boxData.push_back(vec3(res.x, 0.0f, res.z));    boxData.push_back(vec3(res.x, 0.0f, 0.0f));

    data = createVBO(data, (float*)&boxData[0], boxData.size()*3, LINES, key);
    boxData.clear();
    vbos.push_back(data);
    vbokeys["boundingbox"] = vbos.size()-1;

    int numberOfSolidObjects = sim->getScene()->getSolidObjects().size();
    vector<objCore::objContainer*> solids = sim->getScene()->getSolidObjects();
    for(int i=0; i<numberOfSolidObjects; i++){
        vboData objectdata;
        objectdata.color = vec3(1, 0, 0);
        key = "solid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(solids[i], objectdata.color, key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
    }

    int numberOfLiquidObjects = sim->getScene()->getLiquidObjects().size();
    vector<objCore::objContainer*> liquids = sim->getScene()->getLiquidObjects();
    for(int i=0; i<numberOfLiquidObjects; i++){
        vboData objectdata;
        objectdata.color = vec3(0, 0, 1);
        key = "liquid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(liquids[i], objectdata.color, key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
    }

    return true;
}

vboData viewer::createVBO(vboData data, float* vertices, int numberOfVertices, vbotype type, string key){
    data.size = numberOfVertices;
    glGenBuffers(1, &data.vboID);
    glBindBuffer(GL_ARRAY_BUFFER, data.vboID);
    glBufferData(GL_ARRAY_BUFFER, data.size*sizeof(float), vertices, GL_STATIC_DRAW);
    data.type = type;
    data.key = key;
    return data;
}

vboData viewer::createVBOFromObj(objCore::objContainer* o, vec3 color, string key){
    objCore::obj* oData = o->getObj();
    vboData data;
    vector<vec3> vertexData;
    for(int i=0; i<oData->numberOfPolys; i++){
        vec4 f = oData->polyVertexIndices[i];
        vec3 p0 = oData->vertices[int(f[0])-1];
        vec3 p1 = oData->vertices[int(f[1])-1];
        vec3 p2 = oData->vertices[int(f[2])-1];
        vertexData.push_back(p0);
        vertexData.push_back(p1);
        vertexData.push_back(p2);      
        if(f[3]>0){
            vec3 p3 = oData->vertices[int(f[3])-1];
            vertexData.push_back(p0);
            vertexData.push_back(p3);
            vertexData.push_back(p2);  
        }
    }
    data.color = color;
    data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, TRIANGLES, key);
    vertexData.clear();
    return data;
}

//====================================
// Interaction Callbacks
//====================================

void viewer::errorCallback(int error, const char* description){
    fputs(description, stderr);
}

void viewer::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods){
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
        glfwSetWindowShouldClose(window, GL_TRUE);
        exit(EXIT_SUCCESS);
    }
}
