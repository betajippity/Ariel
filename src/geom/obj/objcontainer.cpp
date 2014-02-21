// ObjCore2.6: An (improved) obj mesh wrangling library. Part of TAKUA Render.
// Written by Yining Karl Li
//
// File: objContainer.cpp
// Implements objContainer.hpp

#include <vector>
#include "objcontainer.hpp"
#include <omp.h>

using namespace objCore;

objContainer::objContainer(string filename){
    mesh = new obj;
    preload(filename);
    load(filename);
    keep = false;
}

objContainer::objContainer(obj* o){
    mesh = o;
    keep = false;
}

objContainer::~objContainer(){
    if(!keep){
        clearObj(mesh);
    }
    delete mesh;
}

void objContainer::load(string filename){
    mesh->vertices = new vec3[mesh->numberOfVertices];
    mesh->normals = new vec3[mesh->numberOfNormals];
    mesh->uvs = new vec2[mesh->numberOfUVs];
    mesh->polyVertexIndices = new vec4[mesh->numberOfPolys];
    mesh->polyNormalIndices = new vec4[mesh->numberOfPolys];
    mesh->polyUVIndices = new vec4[mesh->numberOfPolys];

    int currentVertex = 0;
    int currentNormal = 0;
    int currentUV = 0;
    int currentFace = 0;

    //Load loop
    ifstream fp_in;
    char* fname = (char*)filename.c_str();
    fp_in.open(fname);
    if(fp_in.is_open()){
        while(fp_in.good()){
            string line;
            getline(fp_in, line);
            if(!line.empty()){
                vector<string> tokens = utilityCore::tokenizeString(line, " ");
                if(tokens.size()>1){
                    if(tokens[0].compare("v")==0){
                        mesh->vertices[currentVertex] = vec3(atof(tokens[1].c_str()),
                                                             atof(tokens[2].c_str()),
                                                             atof(tokens[3].c_str()));
                        currentVertex++;
                    }else if(tokens[0].compare("vn")==0){
                        mesh->normals[currentNormal] = vec3(atof(tokens[1].c_str()),
                                                            atof(tokens[2].c_str()),
                                                            atof(tokens[3].c_str()));
                        currentNormal++;
                    }else if(tokens[0].compare("vt")==0){
                        mesh->uvs[currentUV] = vec2(atof(tokens[1].c_str()), 
                                                    atof(tokens[2].c_str()));
                        currentUV++;
                    }else if(tokens[0].compare("f")==0){
                        vec4 vertices(-1,-1,-1,-1);
                        vec4 normals(-1,-1,-1,-1);
                        vec4 uvs(-1,-1,-1,-1);
                        int loops = std::min(4, (int)tokens.size()-1)+1;
                        for(int i=1; i<loops; i++){
                            vector<string> faceindices = utilityCore::tokenizeString(tokens[i], 
                                                                                     "/");
                            vertices[i-1] = atoi(faceindices[0].c_str());
                            if(faceindices.size()==2 && mesh->numberOfNormals>1){
                                normals[i-1] = atoi(faceindices[1].c_str());
                            }else if(faceindices.size()==2 && mesh->numberOfUVs>1){
                                uvs[i-1] = atoi(faceindices[1].c_str());
                            }else if(faceindices.size()==3){
                                uvs[i-1] = atoi(faceindices[1].c_str());
                                normals[i-1] = atoi(faceindices[2].c_str());
                            }
                        }
                        mesh->polyVertexIndices[currentFace] = vertices;
                        mesh->polyNormalIndices[currentFace] = normals;
                        mesh->polyUVIndices[currentFace] = uvs;
                        currentFace++;
                    }
                }               
            }
        }
    } 
    fp_in.close();

    if(mesh->numberOfUVs==0){
        cout << "No UVs found, creating default UVs..." << endl;
        delete [] mesh->uvs;
        mesh->uvs = new vec2[4];
        mesh->uvs[0] = vec2(0,0);
        mesh->uvs[1] = vec2(0,1);
        mesh->uvs[2] = vec2(1,1);
        mesh->uvs[3] = vec2(1,0);
        mesh->numberOfUVs = 4;
        for(int i=0; i<mesh->numberOfPolys; i++){
            mesh->polyUVIndices[i] = vec4(1,2,3,4);
        }
    }

    if(mesh->numberOfNormals==0){
        cout << "No normals found, creating default normals..." << endl;
        delete [] mesh->normals;
        mesh->numberOfNormals = mesh->numberOfPolys;
        mesh->normals = new vec3[mesh->numberOfPolys];
        for(int i=0; i<mesh->numberOfPolys; i++){
            vec4 vertexIndices = mesh->polyVertexIndices[i];
            vec3 v0 = mesh->vertices[(int)vertexIndices[0]-1];
            vec3 v1 = mesh->vertices[(int)vertexIndices[1]-1];
            vec3 v2 = mesh->vertices[(int)vertexIndices[2]-1];
            vec3 n0 = v0-v1;
            vec3 n1 = v2-v1;
            utilityCore::printVec3(n0);
            mesh->normals[i] = normalize(cross(n1,n0));
            mesh->polyNormalIndices[i] = vec4(i+1,i+1,i+1,i+1);
        }
    }
}

bool objContainer::writeObj(string filename){
    ofstream outputFile(filename.c_str());
    if(outputFile.is_open()){
        //write out vertices
        for(int i=0; i<mesh->numberOfVertices; i++){
            vec3 v = mesh->vertices[i];
            outputFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        //write out uvs
        for(int i=0; i<mesh->numberOfUVs; i++){
            vec2 uv = mesh->uvs[i];
            outputFile << "vt " << uv.x << " " << uv.y << "\n";
        }
        //write out normals
        for(int i=0; i<mesh->numberOfNormals; i++){
            vec3 n = mesh->normals[i];
            outputFile << "vn " << n.x << " " << n.y << " " << n.z << "\n";
        }
        //Write out faces
        if(mesh->numberOfUVs==0 && mesh->numberOfNormals==0){
            for(int i=0; i<mesh->numberOfPolys; i++){
                vec4 fv = mesh->polyVertexIndices[i];
                outputFile << "f ";
                outputFile << (int)fv.x << " " << (int)fv.y << " " << (int)fv.z << " ";
                if(fv.w>=0){
                    outputFile << (int)fv.w;
                }
                outputFile << "\n";
            }   
        }else{
            for(int i=0; i<mesh->numberOfPolys; i++){
                vec4 fv = mesh->polyVertexIndices[i];
                vec4 fn = mesh->polyNormalIndices[i];
                vec4 fuv = mesh->polyUVIndices[i];
                outputFile << "f ";
                outputFile << (int)fv.x << "/" << (int)fuv.x << "/" << (int)fn.x << " ";
                outputFile << (int)fv.y << "/" << (int)fuv.y << "/" << (int)fn.y << " ";
                outputFile << (int)fv.z << "/" << (int)fuv.z << "/" << (int)fn.z << " ";
                if(fv.w>=0){
                    outputFile << (int)fv.w << "/" << (int)fuv.w << "/" << (int)fn.w << " ";
                }
                outputFile << "\n";
            }            
        }
        outputFile.close();
        cout << "Wrote obj file to " << filename << endl;
        return true;
    }else{
        cout << "Error: Unable to write to " << filename << endl;
        return false;
    }
}

void objContainer::bakeTransform(mat4 transform){
    #pragma parallel for
    for(int i=0; i<mesh->numberOfVertices; i++){
        mesh->vertices[i] = vec3(transform * vec4(mesh->vertices[i], 1.0f));
    }
    #pragma parallel for
    for(int i=0; i<mesh->numberOfNormals; i++){
        mesh->normals[i] = normalize(vec3(transform * vec4(mesh->normals[i], 0.0f)));
    }
}

void objContainer::keepObj(bool keep){
    this->keep = keep;
}

void objContainer::preload(string filename){
	ifstream fp_in;
	int vertexCount = 0;
    int normalCount = 0;
    int uvCount = 0;
    int faceCount = 0;
	char* fname = (char*)filename.c_str();
    fp_in.open(fname);
    if(fp_in.is_open()){
        while(fp_in.good()){
            string line;
            getline(fp_in, line);
            if(!line.empty()){
                string header = utilityCore::getFirstNCharactersOfString(line, 2);
                if(header.compare("v ")==0){
                    vertexCount++;
                }else if(header.compare("vt")==0){
                    uvCount++;
                }else if(header.compare("vn")==0){
                    normalCount++;
                }else if(header.compare("f ")==0){
                    faceCount++;
                }
            }
        }
    } 
    fp_in.close();
    mesh->numberOfVertices = vertexCount;
    mesh->numberOfNormals = normalCount;
    mesh->numberOfUVs = uvCount;
    mesh->numberOfPolys = faceCount;
}

obj* objContainer::getObj(){
	return mesh;
}