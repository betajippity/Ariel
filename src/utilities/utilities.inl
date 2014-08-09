// UtilityCore: A utility library. Part of the TAKUA Render project.
// Written by Yining Karl Li
// Version 0.5.13.39a
//  
// File: utilities.inl
// A collection/kitchen sink of generally useful functions

#ifndef UTILITIES_INL
#define UTILITIES_INL

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

#define GLM_FORCE_RADIANS

#include <glm/gtc/matrix_transform.inl>
#include <glm/glm.hpp>
#include <iostream>
#include <sys/timeb.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <rmsd/rmsd.h>

//====================================
// Math stuff
//====================================

HOST DEVICE float utilityCore::toRadian(float degree){
    return degree*PI/180.0f;
}

HOST DEVICE float utilityCore::toDegree(float radian){
    return radian*180.0f/PI;
}

float utilityCore::clamp(float f, float min, float max){
    if(f<min){
        return min;
    }else if(f>max){
        return max;
    }else{
        return f;
    }
}

glm::vec3 utilityCore::clampRGB(glm::vec3 color){
    if(color[0]<0){
        color[0]=0;
    }else if(color[0]>255){
        color[0]=255;
    }
    if(color[1]<0){
        color[1]=0;
    }else if(color[1]>255){
        color[1]=255;
    }
    if(color[2]<0){
        color[2]=0;
    }else if(color[2]>255){
        color[2]=255;
    }
    return color;
}

HOST DEVICE bool utilityCore::epsilonCheck(float a, float b){
    if(glm::abs(glm::abs(a)-glm::abs(b))<EPSILON){
        return true;
    }else{
        return false;
    }
}

//====================================
// String wrangling stuff
//====================================

std::string utilityCore::padString(int length, std::string str){
    int strlength = str.length();
    std::string pad = "";
    for(int i=strlength; i<length; i++){
        pad = pad + "0";
    }
    return pad + str;
}

bool utilityCore::replaceString(std::string& str, const std::string& from, const std::string& to){
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string utilityCore::convertIntToString(int number){
    std::stringstream ss;
    ss << number;
    return ss.str();
}

std::vector<std::string> utilityCore::tokenizeString(std::string str, std::string separator){
    std::vector<std::string> results;
    char * cstr, *p;
    std::string strt = str;
    cstr = new char[strt.size()+1];
    std::strcpy (cstr, strt.c_str());
    p=std::strtok (cstr, separator.c_str());
    while (p!=NULL){
        results.push_back(p);
        p=strtok(NULL, separator.c_str());
    }
    delete [] cstr;
    delete [] p;
    return results;
}

std::vector<std::string> utilityCore::tokenizeStringByAllWhitespace(std::string str){
    std::stringstream strstr(str);
    std::istream_iterator<std::string> it(strstr);
    std::istream_iterator<std::string> end;
    std::vector<std::string> results(it, end);
    return results;
}

std::string utilityCore::getLastNCharactersOfString(std::string s, int n){
    return s.substr(s.length()-n, n);
}

std::string utilityCore::getFirstNCharactersOfString(std::string s, int n){
    return s.substr(0, n);
}

//====================================
// Time stuff
//====================================

int utilityCore::getMilliseconds(){
    timeb tb;
    ftime( &tb );
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}

int utilityCore::compareMilliseconds(int referenceTime){
    int elapsedMilliseconds = getMilliseconds() - referenceTime;
    if ( elapsedMilliseconds < 0 ){
        elapsedMilliseconds += 0x100000 * 1000;
    }
    return elapsedMilliseconds;
}

//====================================
// Matrix stuff
//====================================

glm::vec3 utilityCore::calculateKabschRotation(glm::vec3 mov0, glm::vec3 mov1, glm::vec3 mov2,
                                               glm::vec3 ref0, glm::vec3 ref1, glm::vec3 ref2){
    double rotatedList[3][3];
    rotatedList[0][0] = mov0.x; rotatedList[0][1] = mov0.y; rotatedList[0][2] = mov0.z;
    rotatedList[1][0] = mov1.x; rotatedList[1][1] = mov1.y; rotatedList[1][2] = mov1.z;
    rotatedList[2][0] = mov2.x; rotatedList[2][1] = mov2.y; rotatedList[2][2] = mov2.z;
    double referenceList[3][3];
    referenceList[0][0] = ref0.x; referenceList[0][1] = ref0.y; referenceList[0][2] = ref0.z;
    referenceList[1][0] = ref1.x; referenceList[1][1] = ref1.y; referenceList[1][2] = ref1.z;
    referenceList[2][0] = ref2.x; referenceList[2][1] = ref2.y; referenceList[2][2] = ref2.z;
    int listCount = 3;
    glm::vec3 movcom = (mov0+mov1+mov2)/3.0f;
    double mov_com[3];
    mov_com[0] = movcom.x; mov_com[1] = movcom.y; mov_com[2] = movcom.z;
    glm::vec3 refcom = (ref0+ref1+ref2)/3.0f;
    glm::vec3 movtoref = refcom - movcom;
    double mov_to_ref[3];
    mov_to_ref[0] = movtoref[0]; mov_to_ref[1] = movtoref[1]; mov_to_ref[2] = movtoref[2];
    double U[3][3];
    double rmsd;

    calculate_rotation_rmsd(referenceList, rotatedList, listCount, mov_com, mov_to_ref, U, &rmsd);

    float x = atan2( U[1][2], U[2][2]  ) * 180.0f/3.1415926f;
    float y = atan2(-U[0][2], std::sqrt(U[1][2]*U[1][2] + U[2][2]*U[2][2]))  * 180.0f/3.1415926f;
    float z = atan2( U[0][1], U[0][0]) * 180.0f/3.1415926f;

    return glm::vec3(x-180.0f,180.0f-y,180.0f-z);
}

glm::mat4 utilityCore::buildTransformationMatrix(glm::vec3 translation, glm::vec3 rotation,
                                                 glm::vec3 scale){
    glm::mat4 translationMat = buildTranslation(translation);
    glm::mat4 rotationMat = buildRotation(toRadian(rotation.z), glm::vec3(0,0,1));
    rotationMat = rotationMat*buildRotation(toRadian(rotation.y), glm::vec3(0,1,0));
    rotationMat = rotationMat*buildRotation(toRadian(rotation.x), glm::vec3(1,0,0));
    glm::mat4 scaleMat = buildScale(scale);
    return translationMat*rotationMat*scaleMat;
}

glm::mat4 utilityCore::buildInverseTransformationMatrix(glm::vec3 translation, glm::vec3 rotation, 
                                                        glm::vec3 scale){
    glm::mat4 translationMat = buildTranslation(-translation);
    glm::mat4 rotationMat = buildRotation(toRadian(-rotation.x), glm::vec3(1,0,0));
    rotationMat = rotationMat*buildRotation(toRadian(-rotation.y), glm::vec3(0,1,0));
    rotationMat = rotationMat*buildRotation(toRadian(-rotation.z), glm::vec3(0,0,1));
    glm::mat4 scaleMat = buildScale(glm::vec3(1.0f)/scale);
    return scaleMat*rotationMat*translationMat;
}

HOST DEVICE glm::vec4 utilityCore::multiply(glm::mat4 m, glm::vec4 v){
    glm::vec4 r;
    r.x = m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z + m[3][0] * v.w;
    r.y = m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z + m[3][1] * v.w;
    r.z = m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z + m[3][2] * v.w;
    r.w = m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] * v.w;
    return r;
}

//this duplicates GLM functionality, but is necessary to work on CUDA
HOST DEVICE glm::mat4 utilityCore::buildTranslation(glm::vec3 translation){
    glm::mat4 m = glm::mat4();
    m[3][0] = translation[0];
    m[3][1] = translation[1];
    m[3][2] = translation[2];
    return m;
}

HOST DEVICE glm::mat4 utilityCore::buildRotation(float radian, glm::vec3 axis){
    axis = axis/glm::length(axis);
    float a = radian;
    float c = glm::cos(a);
    float s = glm::sin(a);
    glm::mat4 m = glm::mat4();
    m[0][0] = c + (1.0f - c)      * axis.x     * axis.x;
    m[0][1] = (1.0f - c) * axis.x * axis.y + s * axis.z;
    m[0][2] = (1.0f - c) * axis.x * axis.z - s * axis.y;
    m[0][3] = 0.0f;
    m[1][0] = (1.0f - c) * axis.y * axis.x - s * axis.z;
    m[1][1] = c + (1.0f - c) * axis.y * axis.y;
    m[1][2] = (1.0f - c) * axis.y * axis.z + s * axis.x;
    m[1][3] = 0.0f;
    m[2][0] = (1.0f - c) * axis.z * axis.x + s * axis.y;
    m[2][1] = (1.0f - c) * axis.z * axis.y - s * axis.x;
    m[2][2] = c + (1.0f - c) * axis.z * axis.z;
    m[2][3] = 0.0f;
    return m;
}

HOST DEVICE glm::mat4 utilityCore::buildScale(glm::vec3 scale){
    glm::mat4 m = glm::mat4();
    m[0][0] = scale[0];
    m[1][1] = scale[1];
    m[2][2] = scale[2];
    return m;
}

//====================================
// GLM Printers
//====================================

void utilityCore::printMat4(glm::mat4 m){
    std::cout << m[0][0] << " " << m[1][0] << " " << m[2][0] << " " << m[3][0] << std::endl;
    std::cout << m[0][1] << " " << m[1][1] << " " << m[2][1] << " " << m[3][1] << std::endl;
    std::cout << m[0][2] << " " << m[1][2] << " " << m[2][2] << " " << m[3][2] << std::endl;
    std::cout << m[0][3] << " " << m[1][3] << " " << m[2][3] << " " << m[3][3] << std::endl;
}

void utilityCore::printVec4(glm::vec4 m){
    std::cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << std::endl;
}

void utilityCore::printVec3(glm::vec3 m){
    std::cout << m[0] << " " << m[1] << " " << m[2] << std::endl;
}

//====================================
// GL Stuff
//====================================

void utilityCore::fovToPerspective(float fovy, float aspect, float zNear, glm::vec2& xBounds, 
                                   glm::vec2& yBounds){
    yBounds[1] = zNear * tan(fovy*(float)PI/360.0f);
    yBounds[0] = -yBounds[1];
    xBounds[0] = yBounds[0]*aspect;
    xBounds[1] = yBounds[1]*aspect;
}

//====================================
// IO Stuff
//====================================

std::string utilityCore::readFileAsString(std::string filename){
    std::ifstream t(filename.c_str());
    std::string str;
    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);
    str.assign((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    return str;
}

std::string utilityCore::getRelativePath(std::string path){
    std::string relativePath;
    std::vector<std::string> pathTokens = utilityCore::tokenizeString(path, "/");
    for(int i=0; i<pathTokens.size()-1; i++){
        relativePath = relativePath + pathTokens[i] + "/";
    }
    return relativePath;
}

#endif
