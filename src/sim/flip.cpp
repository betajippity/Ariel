// Ariel: FLIP Fluid Simulator
// Written by Yining Karl Li
//
// File: flip.cpp
// Implements the FLIP sim

#include "flip.hpp"
#include "../math/kernels.inl"
#include "particlegridoperations.inl"
#include "particleresampler.inl"
#include "solver.inl"

namespace fluidCore{

flipsim::flipsim(const glm::vec3& maxres, sceneCore::Scene* s, const float& density, 
				 const bool& verbose){
	dimensions = maxres;	
	pgrid = new ParticleGrid(maxres);
	mgrid = createMacgrid(maxres);
	mgrid_previous = createMacgrid(maxres);
	max_density = 0.0f;
	this->density = density;
	scene = s;
	frame = 0;
	stepsize = 0.005f;
	subcell = 1;
	picflipratio = .95f;
	densitythreshold = 0.04f;
	this->verbose = verbose;
}

flipsim::~flipsim(){
	delete pgrid;
	unsigned int particlecount = particles.size();
	for(unsigned int i=0; i<particlecount; i++){
		delete particles[i];
	}
	particles.clear();
	clearMacgrid(mgrid);
}

void flipsim::init(){
	//We need to figure out maximum particle pressure, 
	//so we generate a bunch of temporary particles
	//inside of a known area, sort them back onto the underlying grid, and calculate the density
	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = density/maxd;
 	//generate temp particles
	for(unsigned int i = 0; i < 10; i++){  				//FOR_EACH_CELL
		for(unsigned int j = 0; j < 10; j++){ 
			for(unsigned int k = 0; k < 10; k++){ 
				particle* p = new particle;
				p->p = (glm::vec3(i,j,k) + glm::vec3(0.5f))*h;
				p->type = FLUID;
				p->mass = 1.0f;
				particles.push_back(p);
			}
		}
	}
	pgrid->Sort(particles);
	max_density = 1.0f;
	computeDensity(); 
	max_density = 0.0f;
	//sum densities across particles
	for(unsigned int n=0; n<particles.size(); n++) {
		particle *p = particles[n];
		max_density = glm::max(max_density,p->density);
		delete p;
	}
	particles.clear();

	scene->BuildLevelSets(frame);

	//Generate particles and sort
	scene->GenerateParticles(particles, dimensions, density, pgrid, 0);
	pgrid->Sort(particles);
	pgrid->MarkCellTypes(particles, mgrid.A, density);

	//Remove fluid particles that are stuck in walls
	for(std::vector<particle *>::iterator iter=particles.begin(); 
		iter!=particles.end();) { //NONCHECKED
		particle &p = **iter;
		if(p.type == SOLID){
			iter++;
			continue;
		}
		unsigned int i = glm::min(maxd-1,p.p.x*maxd);
		unsigned int j = glm::min(maxd-1,p.p.y*maxd);
		unsigned int k = glm::min(maxd-1,p.p.z*maxd);
		if( mgrid.A->GetCell(i,j,k) == SOLID ) {
			delete *iter;
			iter = particles.erase(iter);
		} else {
			iter ++;
		}
	}
}

void flipsim::step(bool saveVDB, bool saveOBJ, bool savePARTIO){
	frame++;	
	std::cout << "Simulating Step: " << frame << "..." << std::endl;
	
	scene->BuildLevelSets(frame);
	scene->GenerateParticles(particles, dimensions, density, pgrid, frame);

	pgrid->Sort(particles);
	computeDensity();
	applyExternalForces(); 
	splatParticlesToMACGrid(pgrid, particles, &mgrid);
	pgrid->MarkCellTypes(particles, mgrid.A, density);
	storePreviousGrid();
	enforceBoundaryVelocity(&mgrid);
	project();
	enforceBoundaryVelocity(&mgrid);
	extrapolateVelocity();
	subtractPreviousGrid();
	solvePicFlip();
	advectParticles();
	
	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = density/maxd;
	resampleParticles(pgrid, particles, stepsize, h, dimensions);

	unsigned int particlecount = particles.size();

	//mark particles as inside walls or out of bounds
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particlecount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int p=r.begin(); p!=r.end(); ++p){	
				particles[p]->invalid = false;
				glm::vec3 t = particles[p]->p*maxd;
				if(t.x>dimensions.x || t.y>dimensions.y || t.z>dimensions.z){
					particles[p]->invalid = true;
				}
				if(t.x<0 || t.y<0 || t.z<0){
					particles[p]->invalid = true;
				}
				if(mgrid.A->GetCell(t)==SOLID){
					particles[p]->invalid = true;
				}
			}
		}
	);

	//Remove fluid particles that only valid in this frame
	for(std::vector<particle *>::iterator iter=particles.begin(); iter!=particles.end();) {
		particle &p = **iter;
		if( p.temp ) {
			delete *iter;
			iter = particles.erase(iter);
		} else {
			iter ++;
		}
	}

	//Attempt to push particles in walls out 
	particlecount = particles.size();
	std::vector<glm::vec3> stuckPositions;
	std::vector<particle*> stuckParticles;

	for(unsigned int p=0; p<particlecount; p++){
		if(particles[p]->invalid && particles[p]->type == FLUID){
			stuckParticles.push_back(particles[p]);
			stuckPositions.push_back(particles[p]->p*maxd);
		}
	}

	scene->ProjectPointsToSolidSurface(stuckPositions);

	particlecount = stuckPositions.size();
	for(unsigned int p=0; p<particlecount; p++){
		if(glm::length(stuckPositions[p] - stuckParticles[p]->p*maxd)>0.0001f){
			float penaltyForce = 10.0f;
			glm::vec3 vdir = stuckPositions[p] - stuckParticles[p]->p*maxd;
			stuckParticles[p]->p = stuckPositions[p]/maxd;
			stuckParticles[p]->u = vdir*penaltyForce;
		}
	}

	if(saveVDB || saveOBJ || savePARTIO){
		scene->ExportParticles(particles, maxd, frame, saveVDB, saveOBJ, savePARTIO);
	}
}

void flipsim::advectParticles(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;
	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	unsigned int particleCount = particles.size();

	//update positions
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				if(particles[i]->type == FLUID){
					glm::vec3 velocity = interpolateVelocity(particles[i]->p, &mgrid);
					particles[i]->p += stepsize*velocity;
				}
			}
		}
	);
	pgrid->Sort(particles); //sort

	//apply constraints for outer walls of sim
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int p0=r.begin(); p0!=r.end(); ++p0){	
				float r = 1.0f/maxd;
				if( particles[p0]->type == FLUID ) {
					particles[p0]->p = glm::max(glm::vec3(r),glm::min(glm::vec3(1.0f-r),
												particles[p0]->p));
				}

				particle *p = particles[p0];
				if(p->type == FLUID){
					unsigned int i = glm::min(x-1.0f,p->p.x*maxd);
					unsigned int j = glm::min(y-1.0f,p->p.y*maxd);
					unsigned int k = glm::min(z-1.0f,p->p.z*maxd);			
					std::vector<particle*> neighbors = pgrid->GetCellNeighbors(glm::vec3(i,j,k), 
																			   glm::vec3(1));
					for(int p1=0; p1<neighbors.size(); p1++){
						particle* np = neighbors[p1];
						float re = 1.5f*density/maxd;
						if(np->type == SOLID){
							float dist = glm::length(p->p-np->p); //check this later
							if(dist<re){
								glm::vec3 normal = np->n;
								if(glm::length(normal)<0.0000001f && dist){
									normal = glm::normalize(p->p - np->p);
								}
								p->p += (re-dist)*normal;
								p->u -= glm::dot(p->u, normal) * normal;
							}
						}
					}
				}
			}
		}
	);

	//remove stuck particles
	/*std::vector<unsigned int> repositionList;
	for(unsigned int n=0; n<particles.size(); n++){
		particle& p = *particles[n];
		bool reposition = false;
		if(p.type==FLUID){
			unsigned int i = (unsigned int)glm::min(x-1.0f,0.0f);
			unsigned int j = (unsigned int)glm::min(y-1.0f,0.0f);
			unsigned int k = (unsigned int)glm::min(z-1.0f,0.0f);

			//reposition particles stuck on walls
			if(mgrid.A->getCell(i,j,k)){
				reposition = true;
			}

			i = glm::min(x-3.0f, glm::max(2.0f, p.p.x*maxd));
        	j = glm::min(y-3.0f, glm::max(2.0f, p.p.y*maxd));
       		k = glm::min(z-3.0f, glm::max(2.0f, p.p.z*maxd));
       		if(p.density < densitythreshold && mgrid.A->getCell(i,glm::max(0,j-1),k) == SOLID ||
       										   mgrid.A->getCell(i,glm::min(y-1,j+1),k) == SOLID){
       			reposition = true;
       		}
		}
		if(reposition==true){
			repositionList.push_back(n);
		}
	}*/
}

void flipsim::solvePicFlip(){
	int particleCount = particles.size();

	//store copy of current velocities for later
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				particles[i]->t = particles[i]->u;
			}
		}
	);

	splatMACGridToParticles(particles, &mgrid_previous);

	//set FLIP velocity
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				particles[i]->t = particles[i]->u + particles[i]->t;
			}
		}
	);

	//set PIC velocity
	splatMACGridToParticles(particles, &mgrid);

	//combine PIC and FLIP
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particleCount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				particles[i]->u = (1.0f-picflipratio)*particles[i]->u + 
								  picflipratio*particles[i]->t;
			}
		}
	);
}

void flipsim::storePreviousGrid(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;
	//for every x face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		mgrid_previous.u_x->SetCell(i,j,k,mgrid.u_x->GetCell(i,j,k));
			    	}
			    }
			}
		}
	);
	//for every y face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y+1; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		mgrid_previous.u_y->SetCell(i,j,k,mgrid.u_y->GetCell(i,j,k));
			    	}
			    }
			}
		}
	);
	//for every z face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z+1; ++k){
			    		mgrid_previous.u_z->SetCell(i,j,k,mgrid.u_z->GetCell(i,j,k));
			    	}
			    }
			}
		}
	);
}

void flipsim::subtractPreviousGrid(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;
	//for every x face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		float subx = mgrid.u_x->GetCell(i,j,k) - 
			    					 mgrid_previous.u_x->GetCell(i,j,k);
			    		mgrid_previous.u_x->SetCell(i,j,k,subx);
			    	}
			    }
			}
		}
	);
	//for every y face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y+1; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		float suby = mgrid.u_y->GetCell(i,j,k) - 
			    					 mgrid_previous.u_y->GetCell(i,j,k);
			    		mgrid_previous.u_y->SetCell(i,j,k,suby);
			    	}
			    }
			}
		}
	);
	//for every z face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z+1; ++k){
			    		float subz = mgrid.u_z->GetCell(i,j,k) - 
			    					 mgrid_previous.u_z->GetCell(i,j,k);
			    		mgrid_previous.u_z->SetCell(i,j,k,subz);
			    	}
			    }
			}
		}
	);
}

void flipsim::project(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = 1.0f/maxd; //cell width

	//compute divergence per cell
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
				for(unsigned int j = 0; j < y; ++j){
					for(unsigned int k = 0; k < z; ++k){
						float divergence = (mgrid.u_x->GetCell(i+1, j, k) - 
											mgrid.u_x->GetCell(i, j, k) + 
										    mgrid.u_y->GetCell(i, j+1, k) - 
										    mgrid.u_y->GetCell(i, j, k) +
										    mgrid.u_z->GetCell(i, j, k+1) - 
										    mgrid.u_z->GetCell(i, j, k)
										    ) / h;
						mgrid.D->SetCell(i,j,k,divergence);
					}
				}
			}
		}
	);

	//compute internal level set for liquid surface
	pgrid->BuildSDF(mgrid, density);
	
	Solve(mgrid, subcell, verbose);

	if(verbose){
		std::cout << " " << std::endl; //TODO: no more stupid formatting hacks like this to std::out
	}

	//subtract pressure gradient
	subtractPressureGradient();
}

void flipsim::extrapolateVelocity(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;

	Grid<int>** mark = new Grid<int>*[3];
	Grid<int>** wallmark = new Grid<int>*[3];
	for(unsigned int i=0; i<3; i++){
		mark[i] = new Grid<int>(dimensions, 0);
		wallmark[i] = new Grid<int>(dimensions, 0);
	}

	//initalize temp grids with values
	//for every x face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
						mark[0]->SetCell(i,j,k, (i>0 && mgrid.A->GetCell(i-1,j,k)==FLUID) || 
											    (i<x && mgrid.A->GetCell(i,j,k)==FLUID));
						wallmark[0]->SetCell(i,j,k, (i<=0 || mgrid.A->GetCell(i-1,j,k)==SOLID) && 
												    (i>=x || mgrid.A->GetCell(i,j,k)==SOLID));
			    	}
			    }
			}
		}
	);
	//for every y face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
			  	for(unsigned int j = 0; j < y+1; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
						mark[1]->SetCell(i,j,k, (j>0 && mgrid.A->GetCell(i,j-1,k)==FLUID) || 
											    (j<y && mgrid.A->GetCell(i,j,k)==FLUID));
						wallmark[1]->SetCell(i,j,k, (j<=0 || mgrid.A->GetCell(i,j-1,k)==SOLID) && 
												    (j>=y || mgrid.A->GetCell(i,j,k)==SOLID));
			    	}
			    }
			}
		}
	);
	//for every z face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){		
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z+1; ++k){
						mark[2]->SetCell(i,j,k, (k>0 && mgrid.A->GetCell(i,j,k-1)==FLUID) || 
											    (k<z && mgrid.A->GetCell(i,j,k)==FLUID));
						wallmark[2]->SetCell(i,j,k, (k<=0 || mgrid.A->GetCell(i,j,k-1)==SOLID) && 
												    (k>=z || mgrid.A->GetCell(i,j,k)==SOLID));
			    	}
			    }
			}
		}
	);

	//extrapolate
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y+1; ++j){ 
			    	for(unsigned int k = 0; k < z+1; ++k){
						for(unsigned int n=0; n<3; ++n){
							if(n!=0 && i>x-1){ continue; };
							if(n!=1 && j>y-1){ continue; };
							if(n!=2 && k>z-1){ continue; };
							if(!mark[n]->GetCell(i,j,k) && wallmark[n]->GetCell(i,j,k)){
								unsigned int wsum = 0;
								float sum = 0.0f;
								glm::vec3 q[6] = { glm::vec3(i-1,j,k), glm::vec3(i+1,j,k), 
												   glm::vec3(i,j-1,k), glm::vec3(i,j+1,k), 
												   glm::vec3(i,j,k-1), glm::vec3(i,j,k+1) };
								for(unsigned int qk=0; qk<6; ++qk){
									if(q[qk][0]>=0 && q[qk][0]<x+(n==0) && q[qk][1]>=0 && 
									   q[qk][1]<y+(n==1) && q[qk][2]>=0 && q[qk][2]<z+(n==2) ) {
										if(mark[n]->GetCell(q[qk][0],q[qk][1],q[qk][2])){
											wsum ++;
											if(n==0){
												sum += mgrid.u_x->GetCell(q[qk][0],q[qk][1],
																		  q[qk][2]);
											}else if(n==1){
												sum += mgrid.u_y->GetCell(q[qk][0],q[qk][1],
																		  q[qk][2]);
											}else if(n==2){
												sum += mgrid.u_z->GetCell(q[qk][0],q[qk][1],
																		  q[qk][2]);
											}
										}
									}
								}
								if(wsum){
									if(n==0){
										mgrid.u_x->SetCell(i,j,k,sum/wsum);
									}else if(n==1){
										mgrid.u_y->SetCell(i,j,k,sum/wsum);
									}else if(n==2){
										mgrid.u_z->SetCell(i,j,k,sum/wsum);
									}
								}
							}
						}
			    	}
			    }
			}
		}
	);
	for(unsigned int i=0; i<3; i++){
		delete mark[i];
		delete wallmark[i];
	}
	delete mark;
	delete wallmark;
}

void flipsim::subtractPressureGradient(){
	unsigned int x = (unsigned int)dimensions.x; unsigned int y = (unsigned int)dimensions.y; 
	unsigned int z = (unsigned int)dimensions.z;

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);
	float h = 1.0f/maxd; //cell width

	//for every x face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x+1),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){ 
			    	for(unsigned int k = 0; k < z; ++k){
			    		if(i>0 && i<x){
							float pf = mgrid.P->GetCell(i,j,k);
							float pb = mgrid.P->GetCell(i-1,j,k);
							if(subcell && mgrid.L->GetCell(i,j,k) * 
							   mgrid.L->GetCell(i-1,j,k) < 0.0f){
								if(mgrid.L->GetCell(i,j,k)<0.0f){
									pf = mgrid.P->GetCell(i,j,k);
								}else{
									pf = mgrid.L->GetCell(i,j,k)/
										 glm::min(1.0e-3f,mgrid.L->GetCell(i-1,j,k))*
										 		  mgrid.P->GetCell(i-1,j,k);
								}
								if(mgrid.L->GetCell(i-1,j,k)<0.0f){
									pb = mgrid.P->GetCell(i-1,j,k);
								}else{
									pb = mgrid.L->GetCell(i-1,j,k)/
										 glm::min(1.0e-6f,mgrid.L->GetCell(i,j,k))*
										 		  mgrid.P->GetCell(i,j,k);
								}				
							}
							float xval = mgrid.u_x->GetCell(i,j,k);
							xval -= (pf-pb)/h;
							mgrid.u_x->SetCell(i,j,k,xval);
						}
			    	}
			    }
			} 
		}
	);
	//for every y face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y+1; ++j){
			    	for(unsigned int k = 0; k < z; ++k){
			    		if(j>0 && j<y){
							float pf = mgrid.P->GetCell(i,j,k);
							float pb = mgrid.P->GetCell(i,j-1,k);   
							if(subcell && mgrid.L->GetCell(i,j,k) * 
							   mgrid.L->GetCell(i,j-1,k) < 0.0f){
								if(mgrid.L->GetCell(i,j,k)<0.0f){
									pf = mgrid.P->GetCell(i,j,k);
								}else{
									pf = mgrid.L->GetCell(i,j,k)/
										 glm::min(1.0e-3f,mgrid.L->GetCell(i,j-1,k))*
										 		  mgrid.P->GetCell(i,j-1,k);
								}
								if(mgrid.L->GetCell(i,j-1,k)<0.0f){
									pb = mgrid.P->GetCell(i,j-1,k);
								}else{
									pb = mgrid.L->GetCell(i,j-1,k)/
										 glm::min(1.0e-6f,mgrid.L->GetCell(i,j,k))*
										 		  mgrid.P->GetCell(i,j,k);
								}	
							}
							float yval = mgrid.u_y->GetCell(i,j,k);
							yval -= (pf-pb)/h;
							mgrid.u_y->SetCell(i,j,k,yval);
			    		} 		
			    	} 
			    }
			}
		}
	);
	//for every z face
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,x),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
			  	for(unsigned int j = 0; j < y; ++j){
			    	for(unsigned int k = 0; k < z+1; ++k){
			    		if(k>0 && k<z){
							float pf = mgrid.P->GetCell(i,j,k);
							float pb = mgrid.P->GetCell(i,j,k-1);
							if(subcell && mgrid.L->GetCell(i,j,k) * 
							   mgrid.L->GetCell(i,j,k-1) < 0.0f){
								if(mgrid.L->GetCell(i,j,k)<0.0f){
									pf = mgrid.P->GetCell(i,j,k);
								}else{
									pf = mgrid.L->GetCell(i,j,k)/
										 glm::min(1.0e-3f,mgrid.L->GetCell(i,j,k-1))*
										 		  mgrid.P->GetCell(i,j,k-1);
								}
								if(mgrid.L->GetCell(i,j,k-1)<0.0f){
									pb = mgrid.P->GetCell(i,j,k-1);
								}else{
									pb = mgrid.L->GetCell(i,j,k-1)/
										 glm::min(1.0e-6f,mgrid.L->GetCell(i,j,k))*
										 		  mgrid.P->GetCell(i,j,k);
								}	
							}
							float zval = mgrid.u_z->GetCell(i,j,k);
							zval -= (pf-pb)/h;
							mgrid.u_z->SetCell(i,j,k,zval);
			    		} 		
			    	} 
			    }
			}
		}
	);
}

void flipsim::applyExternalForces(){
	glm::vec3 gravity = glm::vec3(0,-9.8f, 0); //for now, just gravity
	unsigned int particlecount = particles.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particlecount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				particles[i]->u += gravity*stepsize;
			}
		}
	);
}

void flipsim::computeDensity(){

	float maxd = glm::max(glm::max(dimensions.x, dimensions.z), dimensions.y);

	unsigned int particlecount = particles.size();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0,particlecount),
		[=](const tbb::blocked_range<unsigned int>& r){
			for(unsigned int i=r.begin(); i!=r.end(); ++i){	
				//Find neighbours
				if(particles[i]->type==SOLID){
					particles[i]->density = 1.0f;
				}else{
					glm::vec3 position = particles[i]->p;

					position.x = (int)glm::max(0.0f,glm::min((int)maxd-1.0f,(int)maxd*position.x));
					position.y = (int)glm::max(0.0f,glm::min((int)maxd-1.0f,(int)maxd*position.y));
					position.z = (int)glm::max(0.0f,glm::min((int)maxd-1.0f,(int)maxd*position.z));
					std::vector<particle *> neighbors;
					neighbors = pgrid->GetCellNeighbors(position, glm::vec3(1));
					float weightsum = 0.0f;
					unsigned int neighborscount = neighbors.size();
					for(unsigned int m=0; m<neighborscount; m++){
						if(neighbors[m]->type!=SOLID){
							float sqd = mathCore::sqrlength(neighbors[m]->p, particles[i]->p);
							//TODO: figure out a better density smooth approx than density/maxd
							float weight = neighbors[m]->mass * 
										   mathCore::smooth(sqd, 4.0f*density/maxd);
							weightsum = weightsum + weight;
						}
					}
					particles[i]->density = weightsum/max_density;
				}
			}
		}
	);
}

bool flipsim::isCellFluid(const int& x, const int& y, const int& z){
	if(scene->GetLiquidLevelSet()->getCell(x,y,z)<0.0f &&
	   scene->GetSolidLevelSet()->getCell(x,y,z)>=0.0f){
		return true;
	}else{
		return false;
	}
}

std::vector<particle*>* flipsim::getParticles(){
	return &particles;
}

glm::vec3 flipsim::getDimensions(){
	return dimensions;
}

sceneCore::Scene* flipsim::getScene(){
	return scene;	
}

fliptask::fliptask(flipsim* sim, bool dumpVDB, bool dumpOBJ, bool dumpPARTIO){
	m_sim = sim;
	m_dumpPARTIO = dumpPARTIO;
	m_dumpOBJ = dumpOBJ;
	m_dumpVDB = dumpVDB;
}

tbb::task* fliptask::execute(){
	m_sim->step(m_dumpVDB, m_dumpOBJ, m_dumpPARTIO);
	return NULL;
}
}
