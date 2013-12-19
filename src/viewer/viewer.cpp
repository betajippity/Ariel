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

    drawobjects = true;
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
                        particles = sim->getParticles();
                        siminitialized = true;
                    }
                    
                    while(1){
                        sim->step();
                        particles = sim->getParticles();
                        
                    //     frame++;
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

void viewer::mainLoop(){
    while (!glfwWindowShouldClose(window)){

        if(siminitialized){
            vboData data;
            vector<vec3> vertexData;
            vector<vec4> colorData;
            int psize = particles->size();

            vec3 gridSize = sim->getDimensions();
            float maxd = glm::max(glm::max(gridSize.x, gridSize.z), gridSize.y);

            for(int j=0; j<psize; j++){
                if(particles->operator[](j)->type==FLUID){
                    vertexData.push_back(particles->operator[](j)->p*maxd);
                    float c = length(particles->operator[](j)->u)/3.0f;
                    colorData.push_back(vec4(c,c,1,0));
                }
            }
            glDeleteBuffers(1, &data.vboID);
            glDeleteBuffers(1, &data.cboID);
            string key = "fluid";
            data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0], 
                             colorData.size()*4, POINTS, key);
            vertexData.clear();
            colorData.clear();
            vbos[vbokeys["fluid"]] = data;
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
                glBindBuffer(GL_ARRAY_BUFFER, vbos[i].cboID);
                glColorPointer(4, GL_FLOAT, 0, NULL);
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

                vec3 res = sim->getDimensions();
                glTranslatef(-res.x/2, 0, -res.y/2);

                if(i==vbokeys["boundingbox"]){
                    glLineWidth(2.0f);
                }

                if(vbos[i].type==QUADS){
                    if(!(i!=vbokeys["boundingbox"] && drawobjects==false)){
                        glDrawArrays(GL_QUADS, 0, vbos[i].size/3);
                    }
                }else if(vbos[i].type==TRIANGLES){
                    if(!(i!=vbokeys["boundingbox"] && drawobjects==false)){
                        glDrawArrays(GL_TRIANGLES, 0, vbos[i].size/3);
                    }
                }else if(vbos[i].type==LINES){
                    glDrawArrays(GL_LINES, 0, vbos[i].size/3);
                }else if(vbos[i].type==POINTS){
                    glPointSize(5.0f);
                    glDrawArrays(GL_POINTS, 0, vbos[i].size/3);
                }
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
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
    if(glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_Q){
            drawobjects = !drawobjects;
            cam.currentKey = GLFW_KEY_Q;
        }
    }else{
        cam.currentKey = 0;
    }
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
    vector<float> colorData;
    string key = "fluid";
    data = createVBO(data, (float*)&vertexData[0], vertexData.size(), (float*)&colorData[0], 
                     colorData.size(), POINTS, key);
    vertexData.clear();
    vbos.push_back(data);
    vbokeys["fluid"] = vbos.size()-1;

    //buffer for sim bounding box
    key = "boundingbox";
    vec3 res = sim->getDimensions();
    geomCore::cube cubebuilder;
    data = createVBOFromObj(cubebuilder.tesselate(vec3(0), res), vec4(.2,.2,.2,0), key);
    vbos.push_back(data);
    vbokeys["boundingbox"] = vbos.size()-1;

    int numberOfSolidObjects = sim->getScene()->getSolidObjects().size();
    vector<objCore::objContainer*> solids = sim->getScene()->getSolidObjects();
    for(int i=0; i<numberOfSolidObjects; i++){
        vboData objectdata;
        key = "solid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(solids[i], vec4(1,0,0,.75), key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
    }

    int numberOfLiquidObjects = sim->getScene()->getLiquidObjects().size();
    vector<objCore::objContainer*> liquids = sim->getScene()->getLiquidObjects();
    for(int i=0; i<numberOfLiquidObjects; i++){
        vboData objectdata;
        key = "liquid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(liquids[i], vec4(0,0,1,.75), key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
    }

    return true;
}

vboData viewer::createVBO(vboData data, float* vertices, int vertexcount, float* colors,
                          int colorcount, vbotype type, string key){
    data.size = vertexcount;
    glGenBuffers(1, &data.vboID);
    glGenBuffers(1, &data.cboID);
    //bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, data.vboID);
    glBufferData(GL_ARRAY_BUFFER, data.size*sizeof(float), vertices, GL_STATIC_DRAW);
    //bind cbo
    glBindBuffer(GL_ARRAY_BUFFER, data.cboID);
    glBufferData(GL_ARRAY_BUFFER, colorcount*sizeof(float), colors, GL_STATIC_DRAW);

    data.type = type;
    data.key = key;
    return data;
}

vboData viewer::createVBOFromObj(objCore::objContainer* o, vec4 color, string key){
    objCore::obj* oData = o->getObj();
    vboData data;
    vector<vec3> vertexData;
    vector<vec4> colorData;

    vec4 fcheck = oData->polyVertexIndices[0];
    if(int(fcheck[3])-1>0){
        for(int i=0; i<oData->numberOfPolys; i++){
            vec4 f = oData->polyVertexIndices[i];
            vec3 p0 = oData->vertices[int(f[0])-1];
            vec3 p1 = oData->vertices[int(f[1])-1];
            vec3 p2 = oData->vertices[int(f[2])-1];
            vec3 p3 = oData->vertices[int(f[3])-1];
            if(int(f[3])-1<0){
                p3 = p0;
            }
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);
            vertexData.push_back(p3);      
        }
        int vcount = vertexData.size();
        for(int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, QUADS, key);
        colorData.clear();
    }else{
        for(int i=0; i<oData->numberOfPolys; i++){
            vec4 f = oData->polyVertexIndices[i];
            vec3 p0 = oData->vertices[int(f[0])-1];
            vec3 p1 = oData->vertices[int(f[1])-1];
            vec3 p2 = oData->vertices[int(f[2])-1];
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);  
            if(int(f[3])-1>=0){
                vec3 p3 = oData->vertices[int(f[3])-1];
                vertexData.push_back(p3);
                vertexData.push_back(p1);
                vertexData.push_back(p2); 
            }    
        }
        int vcount = vertexData.size();
        for(int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, TRIANGLES, key); 
        colorData.clear();
    }
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
