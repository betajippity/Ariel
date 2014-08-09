// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: camera.cpp
// Implements camera.hpp

#include "camera.hpp"
#include "../utilities/utilities.h"

namespace cameraCore {

//====================================
// Camera Class
//====================================

HOST DEVICE Camera::Camera(){
    SetContents(NULL);
}

HOST DEVICE Camera::Camera(CameraInterface* cam){
    SetContents(cam);
}

HOST DEVICE Camera::~Camera(){

}

HOST DEVICE void Camera::SetContents(CameraInterface* cam){
    m_cam = cam;
}

HOST DEVICE float Camera::GetAperture(const unsigned int& frame){
    return m_cam->GetAperture(frame);
}

HOST DEVICE float Camera::GetFocal(const unsigned int& frame){
    return m_cam->GetFocal(frame);
}

HOST DEVICE glm::vec3 Camera::GetTranslation(const unsigned int& frame){
    return m_cam->GetTranslation(frame);
}

HOST DEVICE glm::vec3 Camera::GetRotation(const unsigned int& frame){
    return m_cam->GetRotation(frame);
}

HOST DEVICE glm::vec3 Camera::GetView(const unsigned int& frame){
    return m_cam->GetView(frame);
}

HOST DEVICE glm::vec3 Camera::GetUp(const unsigned int& frame){
    return m_cam->GetUp(frame);
}

HOST DEVICE float Camera::GetLookat(const unsigned int& frame){
    return m_cam->GetLookat(frame);
}

HOST DEVICE rayCore::Ray Camera::Raycast(const glm::vec2& xy, const glm::vec4& randomSample,
                                         float& frame){
    return m_cam->Raycast(xy, randomSample, frame);
}

HOST DEVICE glm::uvec2 Camera::GetResolution(){
    return m_cam->GetResolution();
}
    
HOST DEVICE glm::vec2 Camera::GetFOV(){
    return m_cam->GetFOV();
}
    
HOST DEVICE unsigned int Camera::GetIterations(){
    return m_cam->GetIterations();
}
    
HOST DEVICE unsigned int Camera::GetTraceDepth(){
    return m_cam->GetTraceDepth();
}
    
//====================================
// CameraFrame Class
//====================================

CameraFrame::CameraFrame(){
    SetContents(glm::vec3(0), glm::vec3(0), 1.0f, 1.0f, 0.0f);
}

CameraFrame::CameraFrame(const glm::vec3& translation, const glm::vec3& rotation,
                         const float& aperture, const float& focal, const float& lookat){
    SetContents(translation, rotation, aperture, focal, lookat);
}

CameraFrame::CameraFrame(const glm::vec3& up, const glm::vec3& view, const glm::vec3& position,
                         const float& aperture, const float& focal, const float& lookat){
    SetContents(up, view, position, aperture, focal, lookat);
}

CameraFrame::~CameraFrame(){

}

void CameraFrame::SetContents(const glm::vec3& up, const glm::vec3& view, const glm::vec3& position,
                              const float& aperture, const float& focal, const float& lookat){
    glm::vec3 rotate = utilityCore::calculateKabschRotation(view, up, glm::vec3(0,0,0),
                                                            glm::vec3(0,0,-1), glm::vec3(0,1,0),
                                                            glm::vec3(0,0,0));
    SetContents(position, rotate, aperture, focal, lookat);
}

void CameraFrame::SetContents(const glm::vec3& translation, const glm::vec3& rotation,
                              const float& aperture, const float& focal, const float& lookat){
    m_translation = translation;
    m_rotation = rotation;
    m_aperture = aperture;
    m_focal = focal;
    m_lookat = lookat;
    glm::mat4 m = utilityCore::buildTransformationMatrix(glm::vec3(0.0f), m_rotation,
                                                         glm::vec3(1.0f,1.0f,1.0f));
    m_view = glm::vec3(utilityCore::multiply(m, glm::vec4(0.0f,0.0f,-1.0f,1.0f)));
    m_up = glm::vec3(utilityCore::multiply(m, glm::vec4(0.0f,1.0f,0.0f,1.0f)));
}
}
