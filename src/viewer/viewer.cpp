// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: viewer.cpp
// Implements viewer.hpp

#include <stb_image/stb_image_write.h>
#include <sstream>
#include "viewer.hpp"
#include "../utilities/utilities.h"

namespace viewerCore{

viewer::viewer(){
    loaded = false;
}

viewer::~viewer(){

}

void viewer::load(fluidCore::flipsim* sim, bool retina){
    load(sim, retina, glm::vec2(1024), glm::vec3(0), glm::vec3(0,0,30), glm::vec2(45.0f), 30.0f);
}

void viewer::load(fluidCore::flipsim* sim, bool retina, glm::vec2 resolution, 
                  glm::vec3 camrotate, glm::vec3 camtranslate, glm::vec2 camfov, float camlookat){
    this->resolution = resolution;

    cam.zoomSpeed = 0.1f;
    cam.panSpeed = 0.1f;
    cam.translate = camtranslate;
    cam.rotate = camrotate;
    cam.lookat = camlookat;
    cam.fov = camfov;

    loaded = true;

    this->sim = sim;
    siminitialized = false;

    drawobjects = true;
    drawInvalid = false;

    dumpFramebuffer = false;
    dumpReady = false;

    dumpVDB = false;
    dumpOBJ = false;

    if(retina){
        framebufferScale = 2;
    }else{
        framebufferScale = 1;
    }

    bitmapData = new unsigned char[3 * (int)resolution.x*framebufferScale * 
                                       (int)resolution.y*framebufferScale];

    pause = false;
}

void viewer::simLoopThread(){
    if(sim->frame==0){
        sim->init();
        particles = sim->getParticles();
        siminitialized = true;
    }
    while(1){
        if(!pause){
            sim->step(dumpVDB, dumpOBJ, dumpPARTIO);
            particles = sim->getParticles();
            if(dumpFramebuffer && dumpReady){
                framebufferWriteLock.lock();
                {
                    saveFrame();
                }
                framebufferWriteLock.unlock();
            }
        }
        siminitialized = true;
    }
}

//returns true if viewer launches and closes successfully, otherwise, returns false
bool viewer::launch(){
    if(loaded==true){
        if(init()==true){
            std::thread simThread([=](){
                simLoopThread();
            });
            mainLoop();
            return true;
        }else{
            std::cout << "Error: GL initialization failed.\n" << std::endl;
            return false;
        }
    }else{
        std::cout << "Error: No sim loaded!\n" << std::endl;
        return false;
    } 
}

void viewer::saveFrame(){
    std::string filename = sim->getScene()->imagePath;
    std::string frameString = utilityCore::padString(4, 
                                                     utilityCore::convertIntToString(sim->frame));
    utilityCore::replaceString(filename, ".png", "."+frameString+".png");

    char anim_filename[2048];
    sprintf(anim_filename, "%s", (char*)filename.c_str());
    stbi_write_png(anim_filename, (int)resolution.x*framebufferScale, 
                                  (int)resolution.y*framebufferScale, 3, bitmapData, 
                                  (int)resolution.x*framebufferScale * 3);
}

//====================================
// Draw/Interaction Loop
//====================================

void viewer::mainLoop(){
    while (!glfwWindowShouldClose(window)){

        if(siminitialized){
            vboData data = vbos[vbokeys["fluid"]];
            std::vector<glm::vec3> vertexData;
            std::vector<glm::vec4> colorData;
            int psize = particles->size();

            glm::vec3 gridSize = sim->getDimensions();
            vertexData.reserve(psize);
            colorData.reserve(psize);
            float maxd = glm::max(glm::max(gridSize.x, gridSize.z), gridSize.y);

            for(int j=0; j<psize; j++){
                if(particles->operator[](j)->type==FLUID){
                    if(!particles->operator[](j)->invalid || 
                       (particles->operator[](j)->invalid && drawInvalid)){
                        vertexData.push_back(particles->operator[](j)->p*maxd);
                        float c = glm::length(particles->operator[](j)->u)/3.0f;
                        c = glm::max(c, 
                                     1.0f*glm::max((.7f-particles->operator[](j)->density),0.0f));
                        bool invalid = particles->operator[](j)->invalid;
                        if(invalid){
                            colorData.push_back(glm::vec4(1,0,0,0));
                        }else{
                            colorData.push_back(glm::vec4(c,c,1,0));
                        }
                    }
                }
            }
            glDeleteBuffers(1, &data.vboID);
            glDeleteBuffers(1, &data.cboID);
            std::string key = "fluid";
            data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, 
                             (float*)&colorData[0], colorData.size()*4, GL_POINTS, key);
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
            // apply the current rotation
            glRotatef(-cam.rotate.x, 1, 0, 0);
            glRotatef(-cam.rotate.y, 0, 1, 0);
            glRotatef(-cam.rotate.z, 0, 0, 1);
            glTranslatef(-cam.translate.x, -cam.translate.y, -cam.translate.z);
            
            for(int i=0; i<vbos.size(); i++){
                glPushMatrix();
                glBindBuffer(GL_ARRAY_BUFFER, vbos[i].vboID);
                glVertexPointer(3, GL_FLOAT, 0, NULL);
                glBindBuffer(GL_ARRAY_BUFFER, vbos[i].cboID);
                glColorPointer(4, GL_FLOAT, 0, NULL);
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

                glm::vec3 res = sim->getDimensions();
                glTranslatef(-res.x/2, 0, -res.y/2);

                if(i==vbokeys["boundingbox"]){
                    glLineWidth(2.0f);
                }

                bool skipDraw = false;
                if(vbos[i].type==GL_QUADS || vbos[i].type==GL_TRIANGLES){
                    glm::vec2 frame = frameranges[vbos[i].key];
                    if(!((frame[0]<0 && frame[1]<0) || (frame[0]<=sim->frame 
                        && sim->frame<=frame[1]))){
                        skipDraw = true;
                    }
                }

                if(vbos[i].type==GL_QUADS){
                    if(!(i!=vbokeys["boundingbox"] && drawobjects==false) && skipDraw==false){
                        glDrawArrays(GL_QUADS, 0, vbos[i].size/3);
                    }
                }else if(vbos[i].type==GL_TRIANGLES){
                    if(!(i!=vbokeys["boundingbox"] && drawobjects==false) && skipDraw==false){
                        glDrawArrays(GL_TRIANGLES, 0, vbos[i].size/3);
                    }
                }else if(vbos[i].type==GL_LINES){
                    glDrawArrays(GL_LINES, 0, vbos[i].size/3);
                }else if(vbos[i].type==GL_POINTS){
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

        if(dumpFramebuffer==true){
            framebufferWriteLock.lock();
            {
                for (unsigned int i=0; i<resolution.y*framebufferScale; ++i){
                    glReadPixels(0,i,(int)resolution.x*framebufferScale, 1, GL_RGB, 
                                 GL_UNSIGNED_BYTE, 
                                 bitmapData + ((int)resolution.x*framebufferScale * 3 * 
                                 (((int)resolution.y*framebufferScale-1)-i)));
                }
            }
            framebufferWriteLock.unlock();
            dumpReady = true; 
        }
    }
    glfwDestroyWindow(window);
    glfwTerminate();
}

void viewer::updateInputs(){
    double x; double y;
    glfwGetCursorPos(window, &x, &y);
    glm::vec2 d;
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
                glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), cam.rotate, 
                                                                     glm::vec3(1,1,1));
                glm::vec3 view = glm::normalize(glm::vec3(m*glm::vec4(0,0,-1,0)));
                glm::vec3 lookatPoint = cam.translate + view*cam.lookat;

                cam.rotate.x += -d.y * cam.rotateSpeed;
                cam.rotate.y += -d.x * cam.rotateSpeed;

                m = utilityCore::buildTransformationMatrix(glm::vec3(0), cam.rotate, 
                                                           glm::vec3(1,1,1));
                view = glm::normalize(glm::vec3(m*glm::vec4(0,0,1,0)));

                cam.translate = lookatPoint +  view*cam.lookat;
            }
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == 1){
                glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), cam.rotate, 
                                                                     glm::vec3(1,1,1));
                glm::vec3 view = glm::normalize(glm::vec3(m*glm::vec4(0,0,-1,0)));

                glm::vec3 lookatPoint = cam.translate + view*cam.lookat;
                cam.translate = cam.translate + (d.y * view * cam.zoomSpeed);
                cam.lookat = glm::length(cam.translate-lookatPoint);
            }
            if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == 1){
                glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), cam.rotate, 
                                                                     glm::vec3(1,1,1));
                glm::vec3 up = glm::normalize(glm::vec3(m*glm::vec4(0,1,0,0)));
                glm::vec3 right = glm::normalize(glm::vec3(m*glm::vec4(-1,0,0,0)));

                cam.translate = cam.translate + (d.x * right * cam.panSpeed);
                cam.translate = cam.translate + (d.y * up * cam.panSpeed);
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
    }else if(glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_R){
            dumpFramebuffer = !dumpFramebuffer;
            dumpReady = false;
            cam.currentKey = GLFW_KEY_R;
            if(dumpFramebuffer){
                std::cout << "\nFramebuffer recording ON.\n" << std::endl;
            }else{
                std::cout << "\nFramebuffer recording OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_P){
            pause = !pause;
            cam.currentKey = GLFW_KEY_P;
            if(pause){
                std::cout << "\nSimulation paused.\n" << std::endl; 
            }
        }
    }else if(glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_V){
            dumpVDB = !dumpVDB;
            cam.currentKey = GLFW_KEY_V;
            if(dumpVDB){
                std::cout << "\nVDB Export ON.\n" << std::endl;
            }else{
                std::cout << "\nVDB Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_O){
            dumpOBJ = !dumpOBJ;
            cam.currentKey = GLFW_KEY_O;
            if(dumpOBJ){
                std::cout << "\nOBJ Export ON.\n" << std::endl;
            }else{
                std::cout << "\nOBJ Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_G){
            dumpPARTIO = !dumpPARTIO;
            cam.currentKey = GLFW_KEY_G;
            if(dumpPARTIO){
                std::cout << "\nPARTIO Export ON.\n" << std::endl;
            }else{
                std::cout << "\nPARTIO Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS){
        if(cam.currentKey!=GLFW_KEY_I){
            drawInvalid = !drawInvalid;
            cam.currentKey = GLFW_KEY_I;
            if(drawInvalid){
                std::cout << "\nDraw out of bound particles ON.\n" << std::endl;
            }else{
                std::cout << "\nDraw out of bound particles OFF.\n" << std::endl;
            }
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
    glm::vec2 fov = cam.fov;

    //Window setup stuff
    glfwSetErrorCallback(errorCallback);
    if (!glfwInit()){
        return false;
    }

    window = glfwCreateWindow(resolution.x, resolution.y, "Ariel: now with 100% more FLIP!", 
                              NULL, NULL);
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
    glm::vec2 xBounds;
    glm::vec2 yBounds;
    utilityCore::fovToPerspective(fov.x, (float)resolution.x/resolution.y, 1, xBounds, yBounds); 
    glFrustum(xBounds[0], xBounds[1], yBounds[0], yBounds[1], 1, 10000000);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_MODELVIEW);

    //dummy buffer for particles
    vboData data;
    std::vector<glm::vec3> vertexData;
    vertexData.push_back(glm::vec3(0,0,0));
    std::vector<glm::vec4> colorData;
    colorData.push_back(glm::vec4(0,0,0,0));
    std::string key = "fluid";
    data = createVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0], 
                     colorData.size()*4, GL_POINTS, key);
    vertexData.clear();
    vbos.push_back(data);
    vbokeys["fluid"] = vbos.size()-1;

    //buffer for sim bounding box
    key = "boundingbox";
    glm::vec3 res = sim->getDimensions();
    geomCore::cube cubebuilder;
    data = createVBOFromObj(cubebuilder.tesselate(glm::vec3(0), res), glm::vec4(.2,.2,.2,0), key);
    vbos.push_back(data);
    vbokeys["boundingbox"] = vbos.size()-1;

    int numberOfSolidObjects = sim->getScene()->getSolidObjects().size();
    std::vector<objCore::objContainer*> solids = sim->getScene()->getSolidObjects();
    for(int i=0; i<numberOfSolidObjects; i++){
        vboData objectdata;
        key = "solid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(solids[i], glm::vec4(1,0,0,.75), key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
        frameranges[key] = sim->getScene()->getSolidFrameRange(i);
    }

    int numberOfLiquidObjects = sim->getScene()->getLiquidObjects().size();
    std::vector<objCore::objContainer*> liquids = sim->getScene()->getLiquidObjects();
    for(int i=0; i<numberOfLiquidObjects; i++){
        vboData objectdata;
        key = "liquid"+utilityCore::convertIntToString(i);
        objectdata = createVBOFromObj(liquids[i], glm::vec4(0,0,1,.75), key);
        vbos.push_back(objectdata);
        vbokeys[key] = vbos.size()-1;
        frameranges[key] = sim->getScene()->getLiquidFrameRange(i);
    }

    return true;
}

vboData viewer::createVBO(vboData data, float* vertices, int vertexcount, float* colors,
                          int colorcount, GLenum type, std::string key){
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

vboData viewer::createVBOFromObj(objCore::objContainer* o, glm::vec4 color, std::string key){
    objCore::obj* oData = o->getObj();
    vboData data;
    std::vector<glm::vec3> vertexData;
    std::vector<glm::vec4> colorData;

    glm::vec4 fcheck = oData->polyVertexIndices[0];
    if(int(fcheck[3])-1>0){
        for(int i=0; i<oData->numberOfPolys; i++){
            glm::vec4 f = oData->polyVertexIndices[i];
            glm::vec3 p0 = oData->vertices[int(f[0])-1];
            glm::vec3 p1 = oData->vertices[int(f[1])-1];
            glm::vec3 p2 = oData->vertices[int(f[2])-1];
            glm::vec3 p3 = oData->vertices[int(f[3])-1];
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
                         colorData.size()*4, GL_QUADS, key);
        colorData.clear();
    }else{
        for(int i=0; i<oData->numberOfPolys; i++){
            glm::vec4 f = oData->polyVertexIndices[i];
            glm::vec3 p0 = oData->vertices[int(f[0])-1];
            glm::vec3 p1 = oData->vertices[int(f[1])-1];
            glm::vec3 p2 = oData->vertices[int(f[2])-1];
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);  
            if(int(f[3])-1>=0){
                glm::vec3 p3 = oData->vertices[int(f[3])-1];
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
                         colorData.size()*4, GL_TRIANGLES, key); 
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
}

