// TAKUA Render: Physically Based Renderer
// Written by Yining Karl Li
//
// File: camera.hpp
// Basic perspective camera, also serves as base camera class

#ifndef CAMERA_HPP
#define CAMERA_HPP

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define SHARED __shared__
#else
#define HOST
#define DEVICE
#define SHARED
#endif

#include <glm/glm.hpp>
#include "../ray/ray.hpp"

enum CameraType{PERSP, FISHEYE};

namespace cameraCore {

//====================================
// Class Declarations
//====================================

class Camera;
class CameraFrame;
class CameraInterface;

class Camera {
    public:
        HOST DEVICE Camera();
        HOST DEVICE Camera(CameraInterface* cam);
        HOST DEVICE ~Camera();

        HOST DEVICE void SetContents(CameraInterface* cam);
        HOST DEVICE float GetAperture(const unsigned int& frame);
        HOST DEVICE float GetFocal(const unsigned int& frame);
        HOST DEVICE glm::vec3 GetTranslation(const unsigned int& frame);
        HOST DEVICE glm::vec3 GetRotation(const unsigned int& frame);
        HOST DEVICE glm::vec3 GetView(const unsigned int& frame);
        HOST DEVICE glm::vec3 GetUp(const unsigned int& frame);
        HOST DEVICE float GetLookat(const unsigned int& frame);
        HOST DEVICE glm::uvec2 GetResolution();
        HOST DEVICE glm::vec2 GetFOV();
        HOST DEVICE unsigned int GetIterations();
        HOST DEVICE unsigned int GetTraceDepth();
        HOST DEVICE rayCore::Ray Raycast(const glm::vec2& xy, const glm::vec4& randomSample, 
                                         float& frame);

        CameraInterface*    m_cam;
        unsigned int        m_id;
};

class CameraInterface {
    public:
        HOST DEVICE CameraInterface(){};
        HOST DEVICE ~CameraInterface(){};

        HOST DEVICE virtual rayCore::Ray Raycast(const glm::vec2& xy, const glm::vec4& sample,
                                                 float& frame) = 0;
        HOST DEVICE virtual float GetAperture(const unsigned int& frame) = 0;
        HOST DEVICE virtual float GetFocal(const unsigned int& frame) = 0;
        HOST DEVICE virtual glm::vec3 GetTranslation(const unsigned int& frame) = 0;
        HOST DEVICE virtual glm::vec3 GetRotation(const unsigned int& frame) = 0;
        HOST DEVICE virtual glm::vec3 GetView(const unsigned int& frame) = 0;
        HOST DEVICE virtual glm::vec3 GetUp(const unsigned int& frame) = 0;
        HOST DEVICE virtual float GetLookat(const unsigned int& frame) = 0;
        HOST DEVICE virtual unsigned int GetID() = 0;
        HOST DEVICE virtual CameraType GetType() = 0;
        HOST DEVICE virtual glm::uvec2 GetResolution() = 0;
        HOST DEVICE virtual glm::vec2 GetFOV() = 0;
        HOST DEVICE virtual unsigned int GetIterations() = 0;
        HOST DEVICE virtual unsigned int GetTraceDepth() = 0;
};

class CameraFrame {
    public:
        CameraFrame();
        CameraFrame(const glm::vec3& translation, const glm::vec3& rotation, 
                    const float& aperture, 
                    const float& focal, const float& lookat);
        CameraFrame(const glm::vec3& up, const glm::vec3& view, const glm::vec3& position,
                    const float& aperture, const float& focal, const float& lookat);
        ~CameraFrame();

        void SetContents(const glm::vec3& translation, const glm::vec3& rotation,
                         const float& aperture, const float& focal, const float& lookat);
        void SetContents(const glm::vec3& up, const glm::vec3& view, const glm::vec3& position, 
                         const float& aperture, const float& focal, const float& lookat);

        glm::vec3       m_translation;
        glm::vec3       m_rotation;
        glm::vec3       m_view;
        glm::vec3       m_up;
        float           m_aperture;
        float           m_focal;
        float           m_lookat;
        unsigned int    m_id;
};
}

#endif
