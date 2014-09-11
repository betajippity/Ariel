// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: viewer.cpp
// Implements viewer.hpp

#include <stb_image/stb_image_write.h>
#include <sstream>
#include "viewer.hpp"
#include "../utilities/utilities.h"
#include "../camera/cameralist.hpp"

namespace viewerCore{

Viewer::Viewer(){
    m_loaded = false;
}

Viewer::~Viewer(){

}

void Viewer::Load(fluidCore::FlipSim* sim, const bool& retina){
    Load(sim, retina, glm::vec2(1024), glm::vec3(0), glm::vec3(0,0,30), glm::vec2(45.0f), 30.0f);
}

void Viewer::Load(fluidCore::FlipSim* sim, const bool& retina, const glm::vec2& resolution, 
                  const glm::vec3& camrotate, const glm::vec3& camtranslate, 
				  const glm::vec2& camfov, const float& camlookat){
    m_resolution = resolution;

    m_cam.m_zoomSpeed = 0.1f;
    m_cam.m_panSpeed = 0.1f;
    m_cam.m_translate = camtranslate;
    m_cam.m_rotate = camrotate;
    m_cam.m_lookat = camlookat;
    m_cam.m_fov = camfov;

    m_loaded = true;
    m_currentFrame = 0;

    m_sim = sim;
    m_siminitialized = false;

    m_drawobjects = true;
    m_drawInvalid = false;

    m_dumpFramebuffer = false;
    m_dumpReady = false;

    m_dumpVDB = false;
    m_dumpOBJ = false;

    if(retina){
        m_framebufferScale = 2;
    }else{
        m_framebufferScale = 1;
    }

    m_bitmapData = new unsigned char[3 * (int)resolution.x*m_framebufferScale * 
                                       (int)resolution.y*m_framebufferScale];

    m_pause = false;
}

void Viewer::SimLoopThread(){
    if(m_sim->m_frame==0){
        m_sim->Init();
        m_particles = m_sim->GetParticles();
        m_siminitialized = true;
    }
    while(1){
        if(!m_pause){
            m_sim->Step(m_dumpVDB, m_dumpOBJ, m_dumpPARTIO);
            m_particles = m_sim->GetParticles();
            if(m_dumpFramebuffer && m_dumpReady){
                m_framebufferWriteLock.lock();
                {
                    SaveFrame();
                }
                m_framebufferWriteLock.unlock();
            }
        }
        m_siminitialized = true;
    }
}

//returns true if viewer launches and closes successfully, otherwise, returns false
bool Viewer::Launch(){
    if(m_loaded==true){
        if(Init()==true){
            std::thread simThread([=](){
                SimLoopThread();
            });
            MainLoop();
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

void Viewer::SaveFrame(){
    std::string filename = m_sim->GetScene()->m_imagePath;
    std::string frameString = utilityCore::padString(4, 
                                                utilityCore::convertIntToString(m_sim->m_frame));
    std::cout << filename << std::endl;
    utilityCore::replaceString(filename, ".png", "."+frameString+".png");

    char anim_filename[2048];
    sprintf(anim_filename, "%s", (char*)filename.c_str());
    stbi_write_png(anim_filename, (int)m_resolution.x*m_framebufferScale, 
                                  (int)m_resolution.y*m_framebufferScale, 3, m_bitmapData, 
                                  (int)m_resolution.x*m_framebufferScale * 3);
}

//====================================
// Draw/Interaction Loop
//====================================

void Viewer::UpdateParticles(){
    //rebuild particles
    if(m_siminitialized){

        m_sim->GetScene()->m_particleLock.lock();

        //dummy buffer for particles
        VboData data;
        std::vector<glm::vec3> vertexData;
        std::vector<glm::vec4> colorData;
        unsigned int psize = m_particles->size();

        glm::vec3 gridSize = m_sim->GetDimensions();
        vertexData.reserve(psize);
        colorData.reserve(psize);
        float maxd = glm::max(glm::max(gridSize.x, gridSize.z), gridSize.y);

        unsigned int lpsize = m_sim->GetScene()->GetLiquidParticleCount();
        // lpsize = m_particles->size();
        for(unsigned int j=0; j<lpsize; j++){
            // if(m_particles->operator[](j)->m_type==FLUID){
                if(!m_particles->operator[](j)->m_invalid || 
                   (m_particles->operator[](j)->m_invalid && m_drawInvalid)){
                    vertexData.push_back(m_particles->operator[](j)->m_p*maxd);
                    float c = glm::length(m_particles->operator[](j)->m_u)/3.0f;
                    c = glm::max(c, 
                                 1.0f*glm::max((.7f-m_particles->operator[](j)->m_density),
                                 0.0f));
                    bool invalid = m_particles->operator[](j)->m_invalid;

                    if(invalid){
                        colorData.push_back(glm::vec4(1,1,0,0));
                    }else if(m_particles->operator[](j)->m_type==SOLID){
                        colorData.push_back(glm::vec4(1,0,0,0));
                    }else{
                        /*if(m_particles->operator[](j)->m_temp2){
                            colorData.push_back(glm::vec4(0,1,0,0));
                        }else if(m_particles->operator[](j)->m_temp){
                            colorData.push_back(glm::vec4(1,0,0,0));
                        }else{*/
                            colorData.push_back(glm::vec4(c,c,1,0));
                        //}
                    }
                }
            // }
        }
        glDeleteBuffers(1, &data.m_vboID);
        glDeleteBuffers(1, &data.m_cboID);
        std::string key = "fluid";
        data = CreateVBO(data, (float*)&vertexData[0], vertexData.size()*3, 
                         (float*)&colorData[0], colorData.size()*4, GL_POINTS, key);
        glm::mat4 m;
        for(int x=0; x<4; x++){
            for(int y=0; y<4; y++){
                data.m_transform[x][y] = m[x][y];
            }
        }
        vertexData.clear();
        colorData.clear();
        m_vbos.push_back(data);
        m_vbokeys["fluid"] = m_vbos.size()-1;

        m_sim->GetScene()->m_particleLock.unlock();
    }
}

void Viewer::MainLoop(){
    while (!glfwWindowShouldClose(m_window)){

        if(m_currentFrame==0 && m_siminitialized==true){
            UpdateParticles();
        }

        //check if frame has incremented; if yes, rebuild vbos
        if(m_currentFrame!=m_sim->m_frame){
            m_vbos.clear();
            m_currentFrame = m_sim->m_frame;
            UpdateParticles();
            UpdateMeshes();
        }

        glClearColor(0.325, 0.325, 0.325, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glEnable (GL_BLEND);
        glBlendFunc (GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);

        glPushMatrix();
            // apply the current rotation
            glRotatef(-m_cam.m_rotate.x, 1, 0, 0);
            glRotatef(-m_cam.m_rotate.y, 0, 1, 0);
            glRotatef(-m_cam.m_rotate.z, 0, 0, 1);
            glTranslatef(-m_cam.m_translate.x, -m_cam.m_translate.y, -m_cam.m_translate.z);
            
            glm::vec3 res = m_sim->GetDimensions();

            for(unsigned int i=0; i<m_vbos.size(); i++){
                glPushMatrix();       
                    glTranslatef(-res.x/2, 0, -res.y/2);
                    glMultMatrixf(m_vbos[i].m_transform[0]);
                    glBindBuffer(GL_ARRAY_BUFFER, m_vbos[i].m_vboID);
                    glVertexPointer(3, GL_FLOAT, 0, NULL);
                    glBindBuffer(GL_ARRAY_BUFFER, m_vbos[i].m_cboID);
                    glColorPointer(4, GL_FLOAT, 0, NULL);
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glEnableClientState(GL_COLOR_ARRAY);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    
                    if(i==m_vbokeys["boundingbox"]){
                        glLineWidth(2.0f);
                    }

                    if(m_vbos[i].m_type==GL_QUADS){
                        if(m_drawobjects==true){
                            glDrawArrays(GL_QUADS, 0, m_vbos[i].m_size/3);
                        }else{
                            if(i==m_vbokeys["boundingbox"]){
                                glDrawArrays(GL_QUADS, 0, m_vbos[i].m_size/3);
                            }
                        }
                    }else if(m_vbos[i].m_type==GL_TRIANGLES){
                        if(m_drawobjects==true){
                            glDrawArrays(GL_TRIANGLES, 0, m_vbos[i].m_size/3);
                        }else{
                            if(i==m_vbokeys["boundingbox"]){
                                glDrawArrays(GL_TRIANGLES, 0, m_vbos[i].m_size/3);
                            }
                        }
                    }else if(m_vbos[i].m_type==GL_LINES){
                        glDrawArrays(GL_LINES, 0, m_vbos[i].m_size/3);
                    }else if(m_vbos[i].m_type==GL_POINTS){
                        glPointSize(5.0f);
                        glDrawArrays(GL_POINTS, 0, m_vbos[i].m_size/3);
                    }
                    glDisableClientState(GL_VERTEX_ARRAY);
                    glDisableClientState(GL_COLOR_ARRAY);

                glPopMatrix();
            }

            /*glPushMatrix();
                glTranslatef(-res.x/2, 0, -res.y/2);
                //draw rays
                for(unsigned int i=0; i<m_rays.size(); i++){
                    glm::vec3 origin = m_rays[i].m_origin;
                    glm::vec3 distantPoint = m_rayendpoints[i];
                    glLineWidth(1.0f);
                    glBegin(GL_LINES);
                        glColor4f(0.0f, 1.0f, 1.0f, 0.0f);
                        glVertex3f(origin.x, origin.y, origin.z);
                        glVertex3f(distantPoint.x, distantPoint.y, distantPoint.z); 
                    glEnd();
                }
            glPopMatrix();*/

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

        glfwSwapBuffers(m_window);
        glfwPollEvents();
        UpdateInputs();

        if(m_dumpFramebuffer==true){
            m_framebufferWriteLock.lock();
            {
                for(unsigned int i=0; i<m_resolution.y*m_framebufferScale; ++i){
                    glReadPixels(0,i,(int)m_resolution.x*m_framebufferScale, 1, GL_RGB, 
                                 GL_UNSIGNED_BYTE, 
                                 m_bitmapData + ((int)m_resolution.x*m_framebufferScale * 3 * 
                                 (((int)m_resolution.y*m_framebufferScale-1)-i)));
                }
            }
            m_framebufferWriteLock.unlock();
            m_dumpReady = true; 
        }
    }
    glfwDestroyWindow(m_window);
    glfwTerminate();
}

void Viewer::UpdateInputs(){
    double x; double y;
    glfwGetCursorPos(m_window, &x, &y);
    glm::vec2 d;
    d.x = float(x-m_cam.m_mouseOld.x);
    d.y = float(y-m_cam.m_mouseOld.y);
    m_cam.m_mouseOld.x = x;
    m_cam.m_mouseOld.y = y;
    if(glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_LEFT) == 1 || 
        glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_RIGHT) == 1 ||
        glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_MIDDLE) == 1){

        bool doCamera = false;

        if(glfwGetKey(m_window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS || 
           glfwGetKey(m_window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS){
            doCamera = true;
        }
        if(doCamera==true){
            if(glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_LEFT) == 1){
                if(glfwGetKey(m_window, GLFW_KEY_Z) == GLFW_PRESS){
                    glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), m_cam.m_rotate,
                                                                         glm::vec3(1,1,1));
                    glm::vec3 up = glm::normalize(glm::vec3(m*glm::vec4(0,1,0,0)));
                    glm::vec3 right = glm::normalize(glm::vec3(m*glm::vec4(-1,0,0,0)));
                    
                    m_cam.m_translate = m_cam.m_translate + (d.x * right * m_cam.m_panSpeed);
                    m_cam.m_translate = m_cam.m_translate + (d.y * up * m_cam.m_panSpeed);
                }else{
                    glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), m_cam.m_rotate,
                                                                         glm::vec3(1,1,1));
                    glm::vec3 view = glm::normalize(glm::vec3(m*glm::vec4(0,0,-1,0)));
                    glm::vec3 lookatPoint = m_cam.m_translate + view*m_cam.m_lookat;

                    m_cam.m_rotate.x += -d.y * m_cam.m_rotateSpeed;
                    m_cam.m_rotate.y += -d.x * m_cam.m_rotateSpeed;

                    m = utilityCore::buildTransformationMatrix(glm::vec3(0), m_cam.m_rotate, 
                                                               glm::vec3(1,1,1));
                    view = glm::normalize(glm::vec3(m*glm::vec4(0,0,1,0)));

                    m_cam.m_translate = lookatPoint +  view*m_cam.m_lookat;
                }
            }
            if(glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_RIGHT) == 1){
                glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), m_cam.m_rotate, 
                                                                     glm::vec3(1,1,1));
                glm::vec3 view = glm::normalize(glm::vec3(m*glm::vec4(0,0,-1,0)));

                glm::vec3 lookatPoint = m_cam.m_translate + view*m_cam.m_lookat;
                m_cam.m_translate = m_cam.m_translate + (d.y * view * m_cam.m_zoomSpeed);
                m_cam.m_lookat = glm::length(m_cam.m_translate-lookatPoint);
            }
            if(glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_MIDDLE) == 1){
                glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0), m_cam.m_rotate, 
                                                                     glm::vec3(1,1,1));
                glm::vec3 up = glm::normalize(glm::vec3(m*glm::vec4(0,1,0,0)));
                glm::vec3 right = glm::normalize(glm::vec3(m*glm::vec4(-1,0,0,0)));

                m_cam.m_translate = m_cam.m_translate + (d.x * right * m_cam.m_panSpeed);
                m_cam.m_translate = m_cam.m_translate + (d.y * up * m_cam.m_panSpeed);
            } 
        }else{
            if(glfwGetMouseButton(m_window, GLFW_MOUSE_BUTTON_LEFT) == 1){
                glm::vec3 res = m_sim->GetDimensions();
                glm::vec3 translate = m_cam.m_translate + glm::vec3(res.x/2, 0, res.y/2);
                cameraCore::CameraFrame* camFrame = new cameraCore::CameraFrame(translate, 
                                                    m_cam.m_rotate, 0.001f, 
                                                    100000000.0f, m_cam.m_lookat);
                cameraCore::PerspectiveCamera perspCam(m_resolution, m_cam.m_fov, 
                                                       1, 1, 1, &camFrame);
                m_rays.clear();
                m_rayendpoints.clear();
                float frame = (float)m_sim->m_frame;

                rayCore::Ray r = perspCam.Raycast(glm::vec2(x,y), glm::vec4(0), 
                                                  frame); 
                float distance = 1000000000.0f;
                std::vector<geomCore::Geom*> liquids = m_sim->GetScene()->GetSolidGeoms();
                unsigned int liquidCount = liquids.size();
                for(unsigned int i=0; i<liquidCount; i++){
                    spaceCore::HitCountTraverseAccumulator result(r.m_origin);
                    geomCore::Geom* geom = liquids[i];
                    geom->Intersect(r, result);
                    // std::cout << i << " " << result.m_numberOfHits << std::endl;
                    if(result.m_intersection.m_hit==true){
                        float hitDistance = glm::length(r.m_origin - result.m_intersection.m_point);
                        distance = glm::min(hitDistance, distance);
                    }
                }
                m_rayendpoints.push_back(r.GetPointAlongRay(distance));
                m_rays.push_back(r);

                delete camFrame;
            }
        }
    }
    if(glfwGetKey(m_window, GLFW_KEY_Q) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_Q){
            m_drawobjects = !m_drawobjects;
            m_cam.m_currentKey = GLFW_KEY_Q;
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_R) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_R){
            m_dumpFramebuffer = !m_dumpFramebuffer;
            m_dumpReady = false;
            m_cam.m_currentKey = GLFW_KEY_R;
            if(m_dumpFramebuffer){
                std::cout << "\nFramebuffer recording ON.\n" << std::endl;
            }else{
                std::cout << "\nFramebuffer recording OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_P) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_P){
            m_pause = !m_pause;
            m_cam.m_currentKey = GLFW_KEY_P;
            if(m_pause){
                std::cout << "\nSimulation paused.\n" << std::endl; 
            }
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_V) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_V){
            m_dumpVDB = !m_dumpVDB;
            m_cam.m_currentKey = GLFW_KEY_V;
            if(m_dumpVDB){
                std::cout << "\nVDB Export ON.\n" << std::endl;
            }else{
                std::cout << "\nVDB Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_O) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_O){
            m_dumpOBJ = !m_dumpOBJ;
            m_cam.m_currentKey = GLFW_KEY_O;
            if(m_dumpOBJ){
                std::cout << "\nOBJ Export ON.\n" << std::endl;
            }else{
                std::cout << "\nOBJ Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_G) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_G){
            m_dumpPARTIO = !m_dumpPARTIO;
            m_cam.m_currentKey = GLFW_KEY_G;
            if(m_dumpPARTIO){
                std::cout << "\nPARTIO Export ON.\n" << std::endl;
            }else{
                std::cout << "\nPARTIO Export OFF.\n" << std::endl;
            }
        }
    }else if(glfwGetKey(m_window, GLFW_KEY_I) == GLFW_PRESS){
        if(m_cam.m_currentKey!=GLFW_KEY_I){
            m_drawInvalid = !m_drawInvalid;
            m_cam.m_currentKey = GLFW_KEY_I;
            if(m_drawInvalid){
                std::cout << "\nDraw out of bound particles ON.\n" << std::endl;
            }else{
                std::cout << "\nDraw out of bound particles OFF.\n" << std::endl;
            }
        }
    }else{
        m_cam.m_currentKey = 0;
    }
}

void Viewer::UpdateMeshes(){
    //buffer for sim bounding box
    std::string key = "boundingbox";
    glm::vec3 res = m_sim->GetDimensions();
    geomCore::CubeGen cubebuilder;
    objCore::Obj* simboundingbox = new objCore::Obj();
    cubebuilder.Tesselate(simboundingbox, glm::vec3(0), res);
    VboData data = CreateVBOFromObj(simboundingbox, glm::vec4(.2,.2,.2,0), key);
    glm::mat4 m;
    for(int x=0; x<4; x++){
        for(int y=0; y<4; y++){
            data.m_transform[x][y] = m[x][y];
        }
    }
    m_vbos.push_back(data);
    m_vbokeys["boundingbox"] = m_vbos.size()-1;

    std::vector<geomCore::Geom*> solids = m_sim->GetScene()->GetSolidGeoms();
    unsigned int numberOfSolidObjects = solids.size();
    for(unsigned int i=0; i<numberOfSolidObjects; i++){
        key = "solid_"+utilityCore::convertIntToString(i);
        if(solids[i]->GetType()==MESH){
            if(CreateVBOFromMeshContainer(dynamic_cast<geomCore::MeshContainer*>
                                          (solids[i]->m_geom), (float)m_currentFrame, 
                                          glm::vec4(1,0,0,.75), key, data)==true){
                m_vbos.push_back(data);
                m_vbokeys[key] = m_vbos.size()-1;
            }
        }else if(solids[i]->GetType()==ANIMMESH){
            if(CreateVBOFromAnimmeshContainer(dynamic_cast<geomCore::AnimatedMeshContainer*>
                                              (solids[i]->m_geom), (float)m_currentFrame, 
                                              glm::vec4(1,0,0,.75), key, data)==true){
                m_vbos.push_back(data);
                m_vbokeys[key] = m_vbos.size()-1;
            }
        }
    }

    std::vector<geomCore::Geom*> liquids = m_sim->GetScene()->GetLiquidGeoms();
    unsigned int numberOfLiquidObjects = liquids.size();
    for(unsigned int i=0; i<numberOfLiquidObjects; i++){
        key = "liquid_"+utilityCore::convertIntToString(i);
        if(liquids[i]->GetType()==MESH){
            if(CreateVBOFromMeshContainer(dynamic_cast<geomCore::MeshContainer*>
                                          (liquids[i]->m_geom), (float)m_currentFrame, 
                                          glm::vec4(0,0,1,.75), key, data)==true){
                m_vbos.push_back(data);
                m_vbokeys[key] = m_vbos.size()-1;
            }
        }else if(liquids[i]->GetType()==ANIMMESH){
            if(CreateVBOFromAnimmeshContainer(dynamic_cast<geomCore::AnimatedMeshContainer*>
                                              (liquids[i]->m_geom), (float)m_currentFrame, 
                                              glm::vec4(0,0,1,.75), key, data)==true){
                m_vbos.push_back(data);
                m_vbokeys[key] = m_vbos.size()-1;
            }
        }
    }
}

//====================================
// Init Stuff
//====================================

bool Viewer::Init(){
    //Camera setup stuff
    glm::vec2 fov = m_cam.m_fov;

    //Window setup stuff
    glfwSetErrorCallback(ErrorCallback);
    if(!glfwInit()){
        return false;
    }

    m_window = glfwCreateWindow(m_resolution.x, m_resolution.y, "Ariel: now with 100% more FLIP!", 
                              NULL, NULL);
    if(!m_window){
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(m_window);
    glfwSetKeyCallback(m_window, KeyCallback);   
    glewExperimental = GL_TRUE;
    if(glewInit()!=GLEW_OK){
        return false;   
    }

    //camera stuff
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glm::vec2 xBounds;
    glm::vec2 yBounds;
    utilityCore::fovToPerspective(fov.x, (float)m_resolution.x/m_resolution.y, 1, xBounds, yBounds); 
    glFrustum(xBounds[0], xBounds[1], yBounds[0], yBounds[1], 1, 10000000);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_MODELVIEW);

    UpdateMeshes();

    return true;
}

VboData Viewer::CreateVBO(VboData& data, float* vertices, const unsigned int& vertexcount, 
						  float* colors, const unsigned int& colorcount, const GLenum& type, 
						  const std::string& key){
    data.m_size = vertexcount;
    glGenBuffers(1, &data.m_vboID);
    glGenBuffers(1, &data.m_cboID);
    //bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, data.m_vboID);
    glBufferData(GL_ARRAY_BUFFER, data.m_size*sizeof(float), vertices, GL_STATIC_DRAW);
    //bind cbo
    glBindBuffer(GL_ARRAY_BUFFER, data.m_cboID);
    glBufferData(GL_ARRAY_BUFFER, colorcount*sizeof(float), colors, GL_STATIC_DRAW);

    data.m_type = type;
    data.m_key = key;
    return data;
}

bool Viewer::CreateVBOFromMeshContainer(geomCore::MeshContainer* o, const float& frame, 
                                        const glm::vec4& color, const std::string& key,
                                        VboData& data){
    glm::mat4 transform;
    glm::mat4 inversetransform;
    if(o->GetTransforms(frame, transform, inversetransform)==false){
        return false;
    }
    data = CreateVBOFromObj(&o->GetMeshFrame(frame)->m_basegeom, color, key);
    for(int x=0; x<4; x++){
        for(int y=0; y<4; y++){
            data.m_transform[x][y] = transform[x][y];
        }
    }
    return true;
}

bool Viewer::CreateVBOFromAnimmeshContainer(geomCore::AnimatedMeshContainer* o, const float& frame,
                                            const glm::vec4& color, const std::string& key,
                                            VboData& data){
    glm::mat4 transform;
    glm::mat4 inversetransform;
    if(o->GetTransforms(frame, transform, inversetransform)==false){
        return false;
    }
    //get objs that we need and interp weight
    objCore::InterpolatedObj* io = &o->GetMeshFrame(frame)->m_basegeom;
    objCore::Obj* o0 = io->m_obj0;
    objCore::Obj* o1 = io->m_obj1;
    float lerpWeight = o->GetInterpolationWeight(frame);

    //build vertex buffers
    std::vector<glm::vec3> vertexData;
    std::vector<glm::vec4> colorData;

    glm::uvec4 fcheck = o0->m_polyVertexIndices[0];
    if(int(fcheck[3])-1>0){
        for(unsigned int i=0; i<o0->m_numberOfPolys; i++){
            glm::uvec4 f = o0->m_polyVertexIndices[i];
            glm::vec3 p0 = o0->m_vertices[int(f[0])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[0])-1] * lerpWeight;
            glm::vec3 p1 = o0->m_vertices[int(f[1])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[1])-1] * lerpWeight;
            glm::vec3 p2 = o0->m_vertices[int(f[2])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[2])-1] * lerpWeight;
            glm::vec3 p3 = o0->m_vertices[int(f[3])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[3])-1] * lerpWeight;
            if(f[3]-1<0){
                p3 = p0;
            }
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);
            vertexData.push_back(p3);      
        }
        unsigned int vcount = vertexData.size();
        for(unsigned int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = CreateVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, GL_QUADS, key);
        colorData.clear();
    }else{
        for(unsigned int i=0; i<o0->m_numberOfPolys; i++){
            glm::uvec4 f = o0->m_polyVertexIndices[i];
            glm::vec3 p0 = o0->m_vertices[int(f[0])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[0])-1] * lerpWeight;
            glm::vec3 p1 = o0->m_vertices[int(f[1])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[1])-1] * lerpWeight;
            glm::vec3 p2 = o0->m_vertices[int(f[2])-1] * (1.0f-lerpWeight) +
                           o1->m_vertices[int(f[2])-1] * lerpWeight;
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);  
            if(int(f[3])-1>=0){
                glm::vec3 p3 = o0->m_vertices[int(f[3])-1] * (1.0f-lerpWeight) +
                               o1->m_vertices[int(f[3])-1] * lerpWeight;
                vertexData.push_back(p3);
                vertexData.push_back(p1);
                vertexData.push_back(p2); 
            }    
        }
        unsigned int vcount = vertexData.size();
        for(unsigned int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = CreateVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, GL_TRIANGLES, key); 
        colorData.clear();
    }
    vertexData.clear();
    //build transform
    for(int x=0; x<4; x++){
        for(int y=0; y<4; y++){
            data.m_transform[x][y] = transform[x][y];
        }
    }
    return true;
}

VboData Viewer::CreateVBOFromObj(objCore::Obj* o, const glm::vec4& color, const std::string& key){
    VboData data;
    std::vector<glm::vec3> vertexData;
    std::vector<glm::vec4> colorData;

    glm::uvec4 fcheck = o->m_polyVertexIndices[0];
    if(int(fcheck[3])-1>0){
        for(unsigned int i=0; i<o->m_numberOfPolys; i++){
            glm::uvec4 f = o->m_polyVertexIndices[i];
            glm::vec3 p0 = o->m_vertices[int(f[0])-1];
            glm::vec3 p1 = o->m_vertices[int(f[1])-1];
            glm::vec3 p2 = o->m_vertices[int(f[2])-1];
            glm::vec3 p3 = o->m_vertices[int(f[3])-1];
            if(f[3]-1<0){
                p3 = p0;
            }
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);
            vertexData.push_back(p3);      
        }
        unsigned int vcount = vertexData.size();
        for(unsigned int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = CreateVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, GL_QUADS, key);
        colorData.clear();
    }else{
        for(unsigned int i=0; i<o->m_numberOfPolys; i++){
            glm::uvec4 f = o->m_polyVertexIndices[i];
            glm::vec3 p0 = o->m_vertices[int(f[0])-1];
            glm::vec3 p1 = o->m_vertices[int(f[1])-1];
            glm::vec3 p2 = o->m_vertices[int(f[2])-1];
            vertexData.push_back(p0);
            vertexData.push_back(p1);
            vertexData.push_back(p2);  
            if(int(f[3])-1>=0){
                glm::vec3 p3 = o->m_vertices[f[3]-1];
                vertexData.push_back(p3);
                vertexData.push_back(p1);
                vertexData.push_back(p2); 
            }    
        }
        unsigned int vcount = vertexData.size();
        for(unsigned int i=0; i<vcount; i++){
            colorData.push_back(color);
        }
        data = CreateVBO(data, (float*)&vertexData[0], vertexData.size()*3, (float*)&colorData[0],
                         colorData.size()*4, GL_TRIANGLES, key); 
        colorData.clear();
    }
    vertexData.clear();
    return data;
}

//====================================
// Interaction Callbacks
//====================================

void Viewer::ErrorCallback(int error, const char* description){
    fputs(description, stderr);
}

void Viewer::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods){
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
        glfwSetWindowShouldClose(window, GL_TRUE);
        exit(EXIT_SUCCESS);
    }
}
}

