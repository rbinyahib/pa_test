#ifndef __OPTIONS_H
#define __OPTIONS_H

#define DIY_USE_SPDLOG 1

#include <stdlib.h>
#include <iostream>

using namespace std;

/*
struct blockTable
{

   int blockID;
   float BlockBox[6];

};
*/
//blockTable* blockBooking;
//std::vector<blockTable> blockBooking;

struct Bounds
{

  float min[3];
  float max[3];


};

namespace options {
    
enum SeedMethod{UNIFORM, RANDOM, RANDOM_BLOCK, RANDOM_BOX, TEST_SEED};
    
class Options
{
public:

    Options();
    ~Options() {}
    
    void Parse(int argc, char **argv);
    bool Validate();

    int numDims;
    int numThreads;
    int numDIYBlocks;
    int numDataBlocks;
    int numBlocksInMemory;
    int numLoadedBlocks;
    int maxRounds;
    int maxParticlesPerRound;
    int maxStepsPerParticle;
    int maxStepsPerRound;
    SeedMethod seedMethod;
    float uniformRate[3];
    int numSeeds;
    float seedBBox[6];
    string filename;
    string foldername;
    int fileSize[3];
    int resampleSize[3];
    int blockDivisions[3];
    bool dumpData;
    unsigned int randSeed;
    int maxAdvectSteps;
    float testSeed[3];
    float advectionStep;
    bool recordParticlePaths;
    bool vtkmON;
    bool resample;
    bool POS;
    bool DOS;
    bool WOR;
};

Options::Options()
{
    numDims = 3;
    numThreads = 1;
    numDIYBlocks = 1;
    numDataBlocks = 1;
    numBlocksInMemory = -1; 
    numLoadedBlocks = 4;
    maxRounds = 0;
    maxParticlesPerRound = 0;
    maxStepsPerParticle = 0;
    maxStepsPerRound = 0;
    maxAdvectSteps = 10000;
    seedMethod = UNIFORM;
    uniformRate[0] = uniformRate[1] = uniformRate[2] = 1;
    numSeeds = 1000;
    fileSize[0] = fileSize[1] = fileSize[2] = -1;
    resampleSize[0] = resampleSize[1] = resampleSize[2] = -1;
    blockDivisions[0] = blockDivisions[1] = blockDivisions[2] = 0;
    dumpData = false;
    advectionStep = 0.1;
    //first several digits of pi make the perfect seed.
    randSeed = 314159265;
    recordParticlePaths = false;
    vtkmON = false;
    POS = false;
    DOS = true;
    WOR = false;
    resample = false;
}

bool
Options::Validate()
{

    //numDataBlocks = numDIYBlocks;
    if (filename.size() == 0)
        throw "Error: no input file specified.";
    if (fileSize[0] < 1 || fileSize[1] < 1 || fileSize[2] < 1)
        throw "Error: input file SIZE not specified";
 //   if (!(blockDivisions[0] == 0 && blockDivisions[1] == 0 && blockDivisions[2] == 0))
   //     if (blockDivisions[0]*blockDivisions[1]*blockDivisions[2] != numDataBlocks)
     //       throw "Error: Block divisions wrong for number of blocks specified";
    if (resampleSize[0] < 1 || resampleSize[1] < 1 || resampleSize[2] < 1)
        for (int i = 0; i < 3; i++) resampleSize[i] = fileSize[i];

    
    //Seed the random number generator.
    srand(randSeed);
    //blockBooking.resize(numDataBlocks);
   // blockBooking = new blockTable[numDataBlocks];
    return true;
}

void
Options::Parse(int argc, char **argv)
{
    for (int i = 0; i < argc; i++)
    {
        string arg(argv[i]);
        if (arg == "--filename")
           {
               filename = argv[++i];
               //std::cerr<<" filename "<<filename<<"\n";
            }
        else if (arg == "--file-dim")
        {
            fileSize[0] = atoi(argv[++i]);
            fileSize[1] = atoi(argv[++i]);
            fileSize[2] = atoi(argv[++i]);
        }
        else if (arg == "--resample-dim")
        {
            resample = true;
            resampleSize[0] = atoi(argv[++i]);
            resampleSize[1] = atoi(argv[++i]);
            resampleSize[2] = atoi(argv[++i]);
        }
        else if (arg == "--folder-time")
        {
               foldername = argv[++i];

         }
        else if(arg == "--numDataBlocks")
        {

             numDataBlocks = atoi(argv[++i]);
        }
        else if (arg == "--block-divisions")
        {
            blockDivisions[0] = atoi(argv[++i]);
            blockDivisions[1] = atoi(argv[++i]);
            blockDivisions[2] = atoi(argv[++i]);
      
           //cerr<<"blockDivisions "<<blockDivisions[0]<<", "<<blockDivisions[1]<<", "<<blockDivisions[2]<<"\n";
        }        
        else if (arg == "--blocks")
            numDIYBlocks = atoi(argv[++i]);
        else if (arg == "--blocks-in-memory")
            numBlocksInMemory = atoi(argv[++i]);
        else if (arg == "--load-blocks")
              numLoadedBlocks = atoi(argv[++i]);
        else if (arg == "--threads")
            numThreads = atoi(argv[++i]);
        else if (arg == "--max-rounds")
            maxRounds = atoi(argv[++i]);
        else if (arg == "--max-particles-round")
            maxParticlesPerRound = atoi(argv[++i]);
        else if (arg == "--max-steps-round")
            maxStepsPerRound = atoi(argv[++i]);
        else if (arg == "--max-advect-steps")
            maxAdvectSteps = atoi(argv[++i]);
        else if (arg == "--advect-step")
            advectionStep = atof(argv[++i]);
        else if (arg == "--seed-uniform")
        {
            seedMethod = UNIFORM;
            uniformRate[0] = atof(argv[++i]);
            uniformRate[1] = atof(argv[++i]);
            uniformRate[2] = atof(argv[++i]);
        }
        else if (arg == "--seed-random-whole")
        {
            seedMethod = RANDOM;
            numSeeds = atoi(argv[++i]);
        }
        else if (arg == "--seed-random-block")
        {
            seedMethod = RANDOM_BLOCK;
            numSeeds = atoi(argv[++i]);
        }
        else if (arg == "--seed-random-box")
        {
            seedMethod = RANDOM_BOX;
            numSeeds = atoi(argv[++i]);
            for (int j = 0; j < 6; j++)
                seedBBox[j] = atof(argv[++i]);
        }
        else if (arg == "--test-seed")
        {
            seedMethod = TEST_SEED;
            testSeed[0] = atof(argv[++i]);
            testSeed[1] = atof(argv[++i]);
            testSeed[2] = atof(argv[++i]);
            numSeeds = 1;
        }
        else if (arg == "--dump-data")
            dumpData = true;
        else if (arg == "--random-seed")
            randSeed = atoi(argv[++i]);
        else if (arg == "--random-seed-time")
            randSeed = time(NULL);
        else if (arg == "--record-particle-path")
            recordParticlePaths = true;
        else if (arg == "--vtkm-on")
           {
                vtkmON = true;
                  cerr<<"***** VTK-m is ON\n";
           }
        else if (arg == "--POS")
          {
               POS = true;
               DOS = false;
               cerr<<"***** Parallelize over Seeds\n";

            }
        else if(arg == "--WOR")
         {  WOR = true;
            DOS = false;
           cerr<<"**** Work Requesting \n";
         } 
 }

    try
    {
        Validate();
    }
    catch (const char *s)
    {
        cerr<<s<<endl;
        exit(-1);
    }
}
}

#endif //__OPTIONS_H
