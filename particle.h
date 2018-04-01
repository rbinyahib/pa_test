#ifndef __PARTICLE_H
#define __PARTICLE_H

#define DIY_USE_SPDLOG 1

#include <iostream>
#include "vec.h"
#include <list>

class Particle : public Vec
{
public:
    Particle() : Vec(), id(-1), seq(0), nSteps(0), status(0),blockID(0), blockCount(0), isMsg(0), requestingID(-1) {}
    Particle(const double *c, int _id) : Vec(c), id(_id), nSteps(0), seq(0), status(0), blockID(0), blockCount(0), isMsg(0), requestingID(-1) {}
    Particle(const Vec &_v, int _id) : Vec(_v), id(_id), nSteps(0), seq(0), status(0), blockID(0), blockCount(0), isMsg(0), requestingID(-1){}
    Particle(const Particle &p) : Vec(p.coords), nSteps(p.nSteps), id(p.id), seq(p.seq), status(p.status), blockID(p.blockID), blockCount(p.blockCount), isMsg(p.isMsg), requestingID(p.requestingID) {}

    int nSteps;
    int id, seq;
    int blockID;
    unsigned long long int  status;
    int blockCount;
    int isMsg;
    int requestingID;
    //std::map<int ,int> blockCount;
/*
 enum ParticleStatus
{
  STATUS_OK = 0x0001,
  TERMINATED = 0x0002,
  EXITED_SPATIAL_BOUNDARY = 0x0004,
  EXITED_TEMPORAL_BOUNDARY = 0x0008,
  STATUS_ERROR = 0x0010
};
*/
};

std::ostream &operator<<(std::ostream &os, Particle const &p)
{
    os<<"Particle id "<<p.id;
//<<" seq "<<p.seq<<" #Steps "<<p.nSteps<<" Coord ["<<p.coords[0]<<" "<<p.coords[1]<<" "<<p.coords[2]<<"] blockID "<<p.blockID<<" isMsg "<<p.isMsg<<" requestingID "<<p.requestingID;
    return os;
}

#if 0
//Serialization of Particle.
namespace diy
{
    template<>
    struct Serialization<Particle>
    {
        static
        void save(diy::BinaryBuffer &bb, const Particle &x)
        {
            diy::Serialization<int>::save(bb, x.id);
            diy::Serialization<int>::save(bb, x.seq);
            diy::Serialization<int>::save(bb, x.nSteps);
            diy::Serialization<int>::save(bb, x.blockID);
            diy::Serialization<int>::save(bb, x.isMsg);
            diy::Serialization<int>::save(bb, x.requestingID);
            diy::Serialization<unsigned long long int>::save(bb, x.status);
            diy::Serialization<double>::save(bb, x.coords[0]);
            diy::Serialization<double>::save(bb, x.coords[1]);
            diy::Serialization<double>::save(bb, x.coords[2]);
        }

        static
        void load(diy::BinaryBuffer &bb, Particle &x)
        {
            diy::Serialization<int>::load(bb, x.id);
            diy::Serialization<int>::load(bb, x.seq);
            diy::Serialization<int>::load(bb, x.nSteps);
            diy::Serialization<int>::load(bb, x.blockID);
            diy::Serialization<int>::load(bb, x.isMsg);
            diy::Serialization<int>::load(bb, x.requestingID);
            diy::Serialization<unsigned long long int>::load(bb, x.status);
            diy::Serialization<double>::load(bb, x.coords[0]);
            diy::Serialization<double>::load(bb, x.coords[1]);
            diy::Serialization<double>::load(bb, x.coords[2]);
        }
    };
}
#endif

#endif //__PARTICLE_H
