#ifndef __INTEGRATOR_H
#define __INTEGRATOR_H

//#define DIY_USE_SPDLOG 1

#include <vtkDataSet.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <vtkStreamTracer.h>

#include "util.h"
#include "options.h"
#include <string>
#include <deque>

using namespace std;
//using namespace std::deque;

class Integrator
{
public:
    enum Status {OK, OUT, ERROR};
    //vector<Particle> outParticles;
    //list<Particle> outParticles;
    deque<Particle> outParticles;
    int maxParticlesPerRound;
    vector<vector<float> > *vec;
   // list<Particle> particles;
    //Particle* outParticles; 
    Integrator(double _h, bool _recordSteps, vector<vector<float> > *_vec):h(_h), recordSteps(_recordSteps), vec(_vec), status(ERROR) {}
    virtual ~Integrator() {}

    virtual void Init() = 0;
    virtual int Push(Particle &p, const options::Options &opts)
    {

        //Copy maxParticlesperRound so vtkm can access it 
        maxParticlesPerRound = opts.maxParticlesPerRound;
        //Cap the number of advection steps taken.
        int nSteps = opts.maxAdvectSteps - p.nSteps;
        // There are 2 maxSteps:
        //    1- per particle
        if (opts.maxStepsPerParticle > 0)
            nSteps = min(nSteps, opts.maxStepsPerParticle);
       //   2- per DIY round 
       else if (opts.maxStepsPerRound > 0)
            nSteps = min(nSteps, opts.maxStepsPerRound);
     
           
        return Push(p, nSteps);
    }
    virtual int Push(Particle &p, int nSteps) = 0;
    //virtual void setParticles(list<Particle> newParticles)
    virtual void setParticles(deque<Particle> newParticles)
   {


         ;//cerr<<"INTERG newParticles "<<newParticles.size()<<"\n";

   }
/*
    virtual int PushParticles(list<Particle> particles, const options::Options &opts)
    {
       //Copy maxParticlesperRound so vtkm can access it 
        maxParticlesPerRound = opts.maxParticlesPerRound;
        //Cap the number of advection steps taken.
        int nSteps = 0;//opts.maxAdvectSteps - p.nSteps;
        // There are 2 maxSteps:
        //    1- per particle
        if (opts.maxStepsPerParticle > 0)
            nSteps = min(nSteps, opts.maxStepsPerParticle);
       //   2- per DIY round 
       else if (opts.maxStepsPerRound > 0)
            nSteps = min(nSteps, opts.maxStepsPerRound);

           return 0; // PushParticles(particles, nSteps);
    }*/
//    virtual int PushParticles(list<Particle> particles, int nSteps) { return 0; } 
 /* 
   virtual int Push(list<Particle>& seedsArr,  const options::Options &opts)
    {

       int nSteps = opts.maxAdvectSteps;
     //  return Push(seedsArr, nSteps); 
      
    }
  */ 
   //virtual int Push(list<Particle>& seedsArr, int nSteps);
    string StatusStr() const
    {
        if (status == OK) return string("OK");
        else if (status == OUT) return string("OUT");
        else if (status == ERROR) return string("ERROR");
        else return string("Unknown");
    }
    
    bool recordSteps;
    double h;
    Status status;
    map<pair<int,int>, vector<float> > samples;    
};

#endif //__INTEGRATOR_H
