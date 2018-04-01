
#include <mpi.h>
#include <iostream>
#include <fstream>
//#include <pcomm/avtParICAlgorithm.h>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include "particle.h"
#include "util.h"
#include "stats.h"
#include "vtkmIntegrator.h"
#include "vtkIntegrator.h"
#include "rk4Integrator.h"
#include "options.h"

#include <iostream>
#include <fstream>
#include  "TimingsManager.h"

using namespace std;

/*
struct Bounds
{

  float min[3];
  float max[3];


};
*/

struct blockTable
{

   int blockID;
   //float BlockBox[6];
   Bounds BlockBox;
  // int count = 0;
};

struct DataBlock
{
   int blockID;
   vector<vector<float> > vel;
   Bounds bounds;

};
//-----------------------------
enum ParticleStatus
{
/*
  STATUS_OK = 0x0001,
  TERMINATED = 0x0002,
  EXITED_SPATIAL_BOUNDARY = 0x0004,
  EXITED_TEMPORAL_BOUNDARY = 0x0008,
  STATUS_ERROR = 0x0010
*/
  STATUS_OK = 1, // 1
  STATUS_TERMINATED = 1 << 1,
  AT_SPATIAL_BOUNDARY = 1 << 2,
  AT_TEMPORAL_BOUNDARY = 1 << 3,
  EXITED_SPATIAL_BOUNDARY = 1 << 4, // 16
  EXITED_TEMPORAL_BOUNDARY = 1 << 5, // 32
  STATUS_ERROR = 1 << 6
};

//Returns 1 if it matches the passed stat and 0 if it does not 
bool CheckBit(const int stat, const ParticleStatus& b)
{
    return ( stat & b) != 0;
}


Bounds box;
static inline float
rand01() { return (float)rand() / (RAND_MAX+1.0f);}

static inline void
randomV(float *v)
{
    v[0] = rand01();
    v[1] = rand01();
    v[2] = rand01();
}

static void
randPoint(Particle &p, const Bounds &box)
{
    float dx = box.max[0]-box.min[0] -1;
    float dy = box.max[1]-box.min[1] -1;
    float dz = box.max[2]-box.min[2] -1;

    p.coords[0] = box.min[0] + rand01()*dx;
    p.coords[1] = box.min[1] + rand01()*dy;
    p.coords[2] = box.min[2] + rand01()*dz;
   
    //p.x = box.min[0] + rand01()*dx;
    //p.y = box.min[1] + rand01()*dy;
    //p.z = box.min[2] + rand01()*dz;
    Particle pa;
}

int maxSteps = 100;
float h = 0.01;

static const int TERMINATED = 0;
static const int BOUNDARY = 1;

static const int MSG_TERMINATED = 1000;
static const int MSG_NEEDPARTICLE = 2000;
//list<Particle> active, boundary, terminated, inActive;
deque<Particle> active, boundary, terminated, inActive,communicate;
std::vector<vtkm::Vec<vtkm::Float32, 3> > vtkmVec;
//Abhishek
vectorFeild* vtkmVec2;
int nVecs;
vector<vector<float> > vel;
vector< vector<int> > BlocksTracking;
vector<blockTable> blockBooking;
deque<DataBlock> LoadedBounds;
int* startInd;
int* endInd;
struct TIMEINFO startTime, endTime;
double elapsed_secsAdvection = 0;
double elapsed_secsComm = 0;
double elapsed_secsIO = 0;
double elapsed_secsIntegInit =0;
double elapsed_secsIntegSteup =0;
double elapsed_secsPush =0;
double elapsed_secsResult =0;
double elapsed_secsVTKRead =0;
double elapsed_secsExtrVec =0;
float dt;

static int
randomRank(int r, int nR)
{
    while (1)
    {
        int rr = rand() % nR;
        if (rr != r)
            return rr;
    }
}

//----------------------------------------------------
bool value_sort(const Particle& pa, const Particle& pb)
{
   // cerr<<"in sort funciton \n";
    //return pa.coords[2] < pb.coords[2];
    return pa.blockID < pb.blockID;
}

bool frequent_sort(const Particle& pa, const Particle& pb)
{

          return pa.blockCount > pb.blockCount;
}
deque<Particle> SortParticels(deque<Particle> particles, int nblocks, deque<DataBlock> LoadedBlocks)// vector<blockTable> LoadedBlocks)
{

   std::map<int, int> myMap;
   for(int i=0; i< nblocks; i++)
      myMap[i] = 0;

   particles[0].blockCount = 1;
   for(int i =0; i< particles.size(); i++)
    {
        myMap[particles[i].blockID] = myMap[particles[i].blockID] + 1;


     }//for

   for(int i=0; i<LoadedBlocks.size(); i++)
      myMap[LoadedBlocks[i].blockID] = INT_MAX;

   for(int i=0; i< particles.size(); i++)
      particles[i].blockCount = myMap[particles[i].blockID];

//for(int i=0; i< nblocks; i++)
// cerr<<"Map "<<i<<" has "<<myMap[i]<<"\n";

   std::sort(particles.begin() , particles.end() ,frequent_sort);
   return particles;

}
//----------------------------------------------------
int GetBlockID(Particle p)
{

   int id =-1;

    for(int i=0; i<blockBooking.size();i++)
      {
           // if(blockBooking[i].BlockBox[0] <= p.coords[0] && p.coords[0] <= blockBooking[i].BlockBox[3] &&
           //    blockBooking[i].BlockBox[1] <= p.coords[1] && p.coords[1] <= blockBooking[i].BlockBox[4] &&
           //    blockBooking[i].BlockBox[2] <= p.coords[2] && p.coords[2] <= blockBooking[i].BlockBox[5]) 
           if(blockBooking[i].BlockBox.min[0] <= p.coords[0] && p.coords[0] <= blockBooking[i].BlockBox.max[0] &&
              blockBooking[i].BlockBox.min[1] <= p.coords[1] && p.coords[1] <= blockBooking[i].BlockBox.max[1] &&
              blockBooking[i].BlockBox.min[2] <= p.coords[2] && p.coords[2] <= blockBooking[i].BlockBox.max[2])
            id = blockBooking[i].blockID;

      }

     return id;
}
//----------------------------------------------------
void printBlock(int id)
{

   cerr<<"Block "<<blockBooking[id].blockID<<" from "<<blockBooking[id].BlockBox.min[0]<<" "<<blockBooking[id].BlockBox.min[1]<<" "<<blockBooking[id].BlockBox.min[2];
   cerr<<" To "<<blockBooking[id].BlockBox.max[0]<<" "<<blockBooking[id].BlockBox.max[1]<<" "<<blockBooking[id].BlockBox.max[2]<<"\n";

}
//----------------------------------------------------
void DivideDataToBlocks(int numBlocks, Bounds DataBox,int *divisonDim)
{


  cerr<<"numBlocks "<<numBlocks<<"\n";

  BlocksTracking.resize(numBlocks);
  int XRatio = (DataBox.max[0] - DataBox.min[0])/divisonDim[0];
  int YRatio = (DataBox.max[1] - DataBox.min[1])/divisonDim[1];
  int ZRatio = (DataBox.max[2] - DataBox.min[2])/divisonDim[2];
 
  int x0,x1,y0,y1,z0,z1, id;
    for(int i=0; i< divisonDim[0]; i++)
       {
         x0 = DataBox.min[0] + i * XRatio;
         x1 = x0 + (XRatio - 1);

         for(int j=0; j< divisonDim[1]; j++)
           {
             y0 = DataBox.min[1] + j * YRatio;
             y1 = y0 + (YRatio - 1);

             for(int k=0; k< divisonDim[2]; k++)
               {
                  z0 = DataBox.min[2] + k * ZRatio;
                  z1 = z0 + (ZRatio - 1);

                  int id = i + j * divisonDim[0] + k * divisonDim[0] * divisonDim[1];

                  blockBooking[id].blockID = id;
                  blockBooking[id].BlockBox.min[0] = x0;
                  blockBooking[id].BlockBox.min[1] = y0;
                  blockBooking[id].BlockBox.min[2] = z0;
                  blockBooking[id].BlockBox.max[0] = x1 + 1;
                  blockBooking[id].BlockBox.max[1] = y1 + 1;
                  blockBooking[id].BlockBox.max[2] = z1 + 1;

                  #if 1
                  if(blockBooking[id].BlockBox.max[0] < DataBox.max[0])
                  blockBooking[id].BlockBox.max[0] +=1;
                  if(blockBooking[id].BlockBox.max[1] < DataBox.max[1])
                  blockBooking[id].BlockBox.max[1] +=1;
                  if(blockBooking[id].BlockBox.max[2] < DataBox.max[2] )
                  blockBooking[id].BlockBox.max[2] +=1;

                 if(blockBooking[id].BlockBox.min[0] > 0)
                 blockBooking[id].BlockBox.min[0] -=1;
                 if(blockBooking[id].BlockBox.min[1] > 0)
                 blockBooking[id].BlockBox.min[1] -=1;
                 if(blockBooking[id].BlockBox.min[2] > 0)
                 blockBooking[id].BlockBox.min[2] -=1;
                  #endif
                  if(id == (blockBooking.size() - 1) )
                  {
                     if(blockBooking[id].BlockBox.max[0] < DataBox.max[0])
                        blockBooking[id].BlockBox.max[0] = DataBox.max[0];
                     if(blockBooking[id].BlockBox.max[1] < DataBox.max[1])
                        blockBooking[id].BlockBox.max[1] = DataBox.max[1];
                     if(blockBooking[id].BlockBox.max[2] < DataBox.max[2])
                        blockBooking[id].BlockBox.max[2] = DataBox.max[2];

                  }

               }//k loop

           }//j loop

       }// i loop

}
//----------------------------------------------------
void AddSeeds(int myRank, int numRanks, const options::Options *opts, int numSeeds, int start, int end)
{


  Particle p;
  Bounds seedsBox;
  if (opts->seedMethod == options::RANDOM ||
                 opts->seedMethod == options::RANDOM_BOX ||
                 opts->seedMethod == options::RANDOM_BLOCK)
   {
      if (opts->seedMethod == options::RANDOM_BOX)
        {
            seedsBox.min[0] = opts->seedBBox[0];
            seedsBox.min[1] = opts->seedBBox[2];
            seedsBox.min[2] = opts->seedBBox[4];
            seedsBox.max[0] = opts->seedBBox[1];
            seedsBox.max[1] = opts->seedBBox[3];
            seedsBox.max[2] = opts->seedBBox[5];
        }
      else if (opts->seedMethod == options::RANDOM)
        seedsBox = box;
 
      srand(opts->randSeed);
      for (int i = 0; i < opts->numSeeds; i++)
       {
           int id = i;
          // if (id % numRanks == myRank)
            //if(myRank == 0)
           
           {
             randPoint(p, seedsBox);
              if(blockBooking[start].BlockBox.min[0] <= p.coords[0] && p.coords[0] < blockBooking[end].BlockBox.max[0] &&
              blockBooking[start].BlockBox.min[1] <= p.coords[1] && p.coords[1] < blockBooking[end].BlockBox.max[1] &&
              blockBooking[start].BlockBox.min[2] <= p.coords[2] && p.coords[2] < blockBooking[end].BlockBox.max[2])
           {
             p.id = id;
             p.nSteps = 0;
             p.blockID = GetBlockID(p);//myRank;
             p.blockCount = 0;
             //cerr<<"P "<<p.coords[0]<<", "<<p.coords[1]<<", "<<p.coords[2]<<" blockID "<<p.blockID<<"\n";
             active.push_back(p);
             //active.push_back(Particle(Vec(rand01(), rand01(), rand01()), id));
             //cerr<<"myRank "<<myRank<<" added "<<p<<"\n";
            }
        }
       }

   }// if random

  else if(opts->seedMethod == options::TEST_SEED)
   {
       p.coords[0] = opts->testSeed[0];
       p.coords[1] = opts->testSeed[1];
       p.coords[2] = opts->testSeed[2];
       p.id = 1;
       p.nSteps = 0;
       p.blockID = myRank;
       p.blockCount = 0;
       active.push_back(p);



    }//if testseed 


        //std::sort(active.begin() , active.end() ,value_sort);
        //active = SortParticels(active, opts->numDataBlocks, LoadedBounds);

       #if 1
        if(myRank == 0)
       for(int i=0; i<active.size(); i++)
          cerr<<"Rank "<<myRank<<" P id "<<active[i].id<<" coord "<<active[i].coords[0]<<", "<<active[i].coords[1]<<", "<<active[i].coords[2]<<" BlockID "<<active[i].blockID<<"\n";
       #endif
  
}
//------------------------------------------------------
static void
ResampleGrid(vtkRectilinearGrid *rgrid, float *ptr, float *samples, int numComponents,
             float *brick_bounds, int *brick_dims, bool outputZonal)
{
    int  i, j, k;

    float x_step = (brick_bounds[1] - brick_bounds[0]) / (brick_dims[0]-1);
    float y_step = (brick_bounds[3] - brick_bounds[2]) / (brick_dims[1]-1);
    float z_step = (brick_bounds[5] - brick_bounds[4]) / (brick_dims[2]-1);

    cerr<<"In ResampleGrid \n";
    if (outputZonal)
    {
        brick_bounds[0] += x_step/2.0;
        brick_dims[0]   -= 1;
        brick_bounds[2] += y_step/2.0;
        brick_dims[1]   -= 1;
        brick_bounds[4] += z_step/2.0;
        brick_dims[2]   -= 1;
    }
    int grid_dims[3];
    rgrid->GetDimensions(grid_dims);

    float *x_prop = new float[brick_dims[0]];
    int   *x_ind  = new int[brick_dims[0]];
    for (i = 0 ; i < brick_dims[0] ; i++)
    {
        float x = brick_bounds[0] + x_step*i;
        if (i == brick_dims[0]-1)
            x = brick_bounds[1];  // floating-point roundoff screws up <,>
        for (j = 0 ; j < grid_dims[0]-1 ; j++)
        {
            float x1 = rgrid->GetXCoordinates()->GetTuple1(j);
            float x2 = rgrid->GetXCoordinates()->GetTuple1(j+1);
            if ((x1 <= x && x <= x2) || (j == grid_dims[0]-2))
            {
                x_ind[i] = j;
                float dist = x2-x1;
                if (dist == 0)
                    x_prop[i] = 0.;
                else
                {
                    float offset = x-x1;
                    x_prop[i] = offset / dist;
                }
                break;
            }
        }
    }
 float *y_prop = new float[brick_dims[1]];
    int   *y_ind  = new int[brick_dims[1]];
    for (i = 0 ; i < brick_dims[1] ; i++)
    {
        float y = brick_bounds[2] + y_step*i;
        if (i == brick_dims[1]-1)
            y = brick_bounds[3];  // floating-point roundoff screws up <,>
        for (j = 0 ; j < grid_dims[1]-1 ; j++)
        {
            float y1 = rgrid->GetYCoordinates()->GetTuple1(j);
            float y2 = rgrid->GetYCoordinates()->GetTuple1(j+1);
            if ((y1 <= y && y <= y2) || (j == grid_dims[1]-2))
            {
                y_ind[i] = j;
                float dist = y2-y1;
                if (dist == 0)
                    y_prop[i] = 0.;
                else
                {
                    float offset = y-y1;
                    y_prop[i] = offset / dist;
                }
                break;
            }
        }
    }

    float *z_prop = new float[brick_dims[2]];
    int   *z_ind  = new int[brick_dims[2]];
    for (i = 0 ; i < brick_dims[2] ; i++)
    {
        float z = brick_bounds[4] + z_step*i;
        if (i == brick_dims[2]-1)
            z = brick_bounds[5];  // floating-point roundoff screws up <,>
        for (j = 0 ; j < grid_dims[2]-1 ; j++)
        {
            float z1 = rgrid->GetZCoordinates()->GetTuple1(j);
            float z2 = rgrid->GetZCoordinates()->GetTuple1(j+1);
            if ((z1 <= z && z <= z2) || (j == grid_dims[2]-2))
            {
                z_ind[i] = j;
                float dist = z2-z1;
                if (dist == 0)
                    z_prop[i] = 0.;
                else
                {
                    float offset = z-z1;
                    z_prop[i] = offset / dist;
                }
                break;
            }
        }
    }

    for (k = 0 ; k < brick_dims[2] ; k++)
    {
        for (j = 0 ; j < brick_dims[1] ; j++)
        {
            for (i = 0 ; i < brick_dims[0] ; i++)
            {
                for ( int c = 0; c < numComponents; c++)
                {
                    // Tri-linear interpolation.
                    float val = 0.;
                    for (int l = 0 ; l < 8 ; l++)
                    {
                        int i_ = x_ind[i] + (l & 1 ? 1 : 0);
                        int j_ = y_ind[j] + (l & 2 ? 1 : 0);
                        int k_ = z_ind[k] + (l & 4 ? 1 : 0);
                        float x_prop_ = (l & 1 ? x_prop[i] : 1.-x_prop[i]);
                        float y_prop_ = (l & 2 ? y_prop[j] : 1.-y_prop[j]);
                        float z_prop_ = (l & 4 ? z_prop[k] : 1.-z_prop[k]);
                        int pt = k_*grid_dims[1]*grid_dims[0]*numComponents +
                                 j_*grid_dims[0]*numComponents +
                                 i_*numComponents +
                                 c;
                        val += x_prop_*y_prop_*z_prop_*ptr[pt];
                    }
                    int pt = k*brick_dims[1]*brick_dims[0]*numComponents + j*brick_dims[0]*numComponents + i*numComponents +c;
                    samples[pt] = val;
                }
            }
        }
    }

    delete [] x_prop;
    delete [] x_ind;
    delete [] y_prop;
    delete [] y_ind;
    delete [] z_prop;
    delete [] z_ind;
}

//------------------------------------------------------
// Create a Mesh and get data info
    vtkRectilinearGrid *
    MakeMesh(int nx, int ny, int nz,
             float x0, float x1,
             float y0, float y1,
             float z0, float z1) //const
    {
        int res[3] = {nx, ny, nz};

        //cerr<<"MakeMesh dims "<<nx<<", "<<ny<<", "<<nz<<"\n";
        vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
        rg->SetDimensions(res);

        vtkFloatArray *x = vtkFloatArray::New();
        vtkFloatArray *y = vtkFloatArray::New();
        vtkFloatArray *z = vtkFloatArray::New();

        x->SetNumberOfTuples(nx);
        y->SetNumberOfTuples(ny);
        z->SetNumberOfTuples(nz);
        rg->SetXCoordinates(x);
        rg->SetYCoordinates(y);
        rg->SetZCoordinates(z);
        float dx = (x1-x0) / (nx-1);
        float dy = (y1-y0) / (ny-1);
        float dz = (z1-z0) / (nz-1);

        for (int i = 0; i < nx; i++)
            x->SetTuple1(i, x0+i*dx);
        x->SetTuple1(nx-1, x1);

        for (int i = 0; i < ny; i++)
            y->SetTuple1(i, y0+i*dy);
        y->SetTuple1(ny-1, y1);
for (int i = 0; i < nz; i++)
            z->SetTuple1(i, z0+i*dz);
        z->SetTuple1(nz-1, z1);

        x->Delete();
        y->Delete();
        z->Delete();

        return rg;
    }

//------------------------------------------------------
vtkRectilinearGrid *
     ReadFile2(const options::Options &opts, const Bounds &domain, const Bounds &currentBounds, int currentdomainID) //const
{
   vtkRectilinearGrid *rg = NULL;
        //------- Copy these lines of code to ReadFile 
        int x0 = currentBounds.min[0], x1 = currentBounds.max[0];
        int y0 = currentBounds.min[1], y1 = currentBounds.max[1];
        int z0 = currentBounds.min[2], z1 = currentBounds.max[2];
        int dx0 = domain.min[0], dx1 = domain.max[0];
        int dy0 = domain.min[1], dy1 = domain.max[1];
        int dz0 = domain.min[2], dz1 = domain.max[2];
         int dx = currentBounds.max[0]-currentBounds.min[0];
        int dy = currentBounds.max[1]-currentBounds.min[1];
        int dz = currentBounds.max[2]-currentBounds.min[2];
        //vtkRectilinearGrid *ds = MakeMesh(dx,dy,dz,x0,x1, y0,y1, z0,z1);

         char str[1024];
        sprintf(str,"%s_%d.vtk",opts.filename.c_str(),currentdomainID);
        cerr<<" filename "<<str<<"\n";

        FILE *fp = fopen(str, "rb");
        if (fp == NULL)
            throw "Error: input file not found.";

        gettimeofday(&startTime, 0);
	// -----------------------
        //cerr<<" MakseMesh done \n";
        int np = dx*dy*dz * 3;
        vtkDataSetReader *rdr = vtkDataSetReader::New();
        rdr->SetFileName(str);
        rdr->Update();
        cerr<<"vtkDataSetReader \n";
        vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
        gettimeofday(&endTime, 0);
        dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
        elapsed_secsVTKRead += dt;
        elapsed_secsIO += dt;
        cerr<<"np "<<np<<"\n";
        //float *p = new float[np];
        //size_t nr = fread(p, sizeof(float), np, fp);
        //fclose(fp);

        cerr<<"Fread done \n";

        #if 0
        float *p = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
        int nPts = dx*dy*dz ;
        vtkFloatArray *var = vtkFloatArray::New();
        var->SetNumberOfComponents(3);
        var->SetNumberOfTuples(nPts);
        float *ptr = (float*)var->GetVoidPointer(0);
        cerr<<"Before loop \n";
        for (int i = 0; i < nPts*3; i++)
        {

          ptr[i] = p[i];

           //cerr<<"ptr "<<ptr[i]<<" ";     
         // if(currentdomainID == 0)
           {
              //fprintf(fp1,"%f ",ptr[i]);
              //fflush(fp1);
            }
        }
        //fclose(fp1);
        cerr<<"After loop \n";
        var->SetName("vec");
        ds->GetPointData()->SetVectors(var);

        cerr<<"SetVectors done \n";
        delete [] p;
        var->Delete();
        #endif
        //return ds;
        return rgrid;


}

//------------------------------------------------------
 vtkRectilinearGrid *
     ReadFile(const options::Options &opts, const Bounds &domain) //const
    {
        vtkRectilinearGrid *rg = NULL;
        rg = MakeMesh(domain.max[0], domain.max[1], domain.max[2],
                      0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

        cerr<<" ReadFile "<<domain.max[0]<<" "<<domain.max[1]<<" "<<domain.max[2]<<"\n";

        FILE *fp = fopen(opts.filename.c_str(), "rb");
        if (fp == NULL)
            throw "Error: input file not found.";

        int n = domain.max[0] * domain.max[1] * domain.max[2] *3;
        float *p = new float[n];
        size_t nr = fread(p, sizeof(float), n, fp);
        fclose(fp);

        int nPts = domain.max[0] * domain.max[1] * domain.max[2];
        vtkFloatArray *var = vtkFloatArray::New();
        var->SetNumberOfComponents(3);
        var->SetNumberOfTuples(nPts);
        float *ptr = (float*)var->GetVoidPointer(0);
        for (int i = 0; i < nPts*3; i++)
        //for (int i = 0; i < nPts*3; )
          {
             ptr[i] = p[i];
          //     ptr[i] = 0;
          // if (i%3 ==0) ptr[i] = 1.0;
        #if 0 
            ptr[i]   = 0;
             ptr[i+1] = 0;
             ptr[i+2] = 1;
             i = i+3;
        #endif
   }

        var->SetName("vec");
        rg->GetPointData()->SetVectors(var);
        delete [] p;
        var->Delete();
        //if (gid == 0)
        // dumpDS(rg, "output/data.vtk");

        return rg;
    }


//------------------------------------------------------
vtkRectilinearGrid *
    Resample(int gid,
             vtkRectilinearGrid *fileDS,
             const Bounds &currentBounds,
             const Bounds &domain) 
    {
        char tmp[32];

        //NOTE: in this function we calculate the bounds for this block
        //    
        int x0 = currentBounds.min[0], x1 = currentBounds.max[0];
        int y0 = currentBounds.min[1], y1 = currentBounds.max[1];
        int z0 = currentBounds.min[2], z1 = currentBounds.max[2];
        int dx0 = domain.min[0], dx1 = domain.max[0];
        int dy0 = domain.min[1], dy1 = domain.max[1];
        int dz0 = domain.min[2], dz1 = domain.max[2];
    //   cerr<<"bounds "<<currentBounds.min[0]<<", "<<currentBounds.min[1]<<", "<<currentBounds.min[2]<<", "<<currentBounds.max[0]<<", "<<currentBounds.max[1]<<", "<<currentBounds.max[2]<<"\n";

        int dx = currentBounds.max[0]-currentBounds.min[0];
        int dy = currentBounds.max[1]-currentBounds.min[1];
        int dz = currentBounds.max[2]-currentBounds.min[2];
        vtkRectilinearGrid *ds = MakeMesh(dx,dy,dz,
                                          x0,x1, y0,y1, z0,z1);

        cerr<<" MakseMesh done \n";
        int np = dx*dy*dz;
        vtkFloatArray *var = vtkFloatArray::New();
        var->SetNumberOfComponents(3);
        var->SetNumberOfTuples(np);
        var->SetName("vec");
        ds->GetPointData()->SetVectors(var);
 float b_bounds[6];
        int b_dims[3] = {x1-x0,y1-y0,z1-z0};
        b_bounds[0] = (float)x0 / (float)(dx1-dx0);
        b_bounds[1] = (float)x1 / (float)(dx1-dx0);
        b_bounds[2] = (float)y0 / (float)(dy1-dy0);
        b_bounds[3] = (float)y1 / (float)(dy1-dy0);
        b_bounds[4] = (float)z0 / (float)(dz1-dz0);
        b_bounds[5] = (float)z1 / (float)(dz1-dz0);

       cerr<<"b_dims "<<b_dims[0]<<" "<<b_dims[1]<<" "<<b_dims[2]<<"\n b_bounds "<<b_bounds[1]<<" "<<b_bounds[3]<<" "<<b_bounds[5]<<"\n";

 ResampleGrid(fileDS, (float*)fileDS->GetPointData()->GetVectors()->GetVoidPointer(0),
                     (float *)var->GetVoidPointer(0),
                     3, b_bounds, b_dims, false);
        cerr<<" ResampleGrid done \n";

   #if 0
        for (int i = 0; i < np; i++)
        {
           
            var->SetComponent(i, 0, 0.0);
            var->SetComponent(i, 1, 0.0);
            var->SetComponent(i, 2, 1);
            
            /*
            var->SetComponent(i, 0, rand01()*rand01());
            var->SetComponent(i, 1, rand01()*rand01());
            var->SetComponent(i, 2, rand01());
            */
        }
       #endif

        //dumpData initallty is false  
       // if (dumpData)
      //  {
            sprintf(tmp, "output/output_%02d.vtk", gid);
            dumpDS(ds, tmp, true);
       // }

        var->Delete();
        return ds;
    }

//------------------------------------------------------
void GetData(int myRank, const options::Options &opts, const Bounds &currentdomain,int currentdomainID)
{

        Bounds fileBounds;
        fileBounds.min[0] = fileBounds.min[1] = fileBounds.min[2] = 0;
        fileBounds.max[0] = opts.fileSize[0]; 
        fileBounds.max[1] = opts.fileSize[1]; 
        fileBounds.max[2] = opts.fileSize[2]; 
   	//vtkRectilinearGrid *fileDS = ReadFile(opts, fileBounds);
   	//vtkRectilinearGrid *fileDS = ReadFile(opts, currentdomain);
   	//vtkRectilinearGrid *resampDS = Resample(myRank, fileDS, domain, box);// domain);
        cerr<<" box "<<box.max[0]<<" "<<box.max[1]<<" "<<box.max[2]<<"\n";
        cerr<<" domain "<<currentdomain.max[0]<<" "<<currentdomain.max[1]<<" "<<currentdomain.max[2]<<"\n";
   	//vtkRectilinearGrid *resampDS = Resample(myRank, fileDS, currentdomain, box);
        vtkRectilinearGrid *resampDS = ReadFile2(opts, box, currentdomain, currentdomainID);
        cerr<<" Resample done \n";
   	vtkDataArray *var = resampDS->GetPointData()->GetVectors();
   	nVecs = var->GetNumberOfTuples();
        gettimeofday(&startTime, 0);
        vel.resize(3);
        vel[0].resize(nVecs);
        vel[1].resize(nVecs);
        vel[2].resize(nVecs);
        float *ptr = (float *) var->GetVoidPointer(0);
        float *vel_x = &(vel[0][0]);
        float *vel_y = &(vel[1][0]);
        float *vel_z = &(vel[2][0]);
        vtkm::Float32 x, y, z;
        //Abhishek
        vtkmVec2 = new vectorFeild[nVecs];
        for (int i = 0; i < nVecs; i++)
        {
/*
          vel[0][i] = var->GetComponent(i, 0);
          vel[1][i] = var->GetComponent(i, 1);
          vel[2][i] = var->GetComponent(i, 2);
 */
          *vel_x = *ptr;
          vtkmVec2[i].vec[0] = *ptr;
          x = *ptr;
          vel_x++;
          ptr++;
         *vel_y = *ptr;
          vtkmVec2[i].vec[1] = *ptr;
          y = *ptr; 
          vel_y++;
          ptr++;
          *vel_z = *ptr;
           vtkmVec2[i].vec[2]= *ptr;
          z = *ptr;
          vel_z++;
          ptr++;

          vtkmVec.push_back(vtkm::Vec<vtkm::Float32, 3>(x,y,z));
          
        }
        //dumpData = dumpData;
        //fileDS->Delete();
        cerr<<"Done reading \n";
        gettimeofday(&endTime, 0);
        dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
        elapsed_secsExtrVec += dt;
        //cerr<<"elapsed_secsExtrVec "<<dt<<"\n";
        elapsed_secsIO += dt;
        resampDS->Delete();

      cerr<<"End of get data \n";
}
//------------------------------------------------------
int  isLoaded(int id)
{
  for(int i=0; i< LoadedBounds.size(); i++)
      if(id == LoadedBounds[i].blockID)
             return i;


    return -1;
}
//------------------------------------------------------
int ownerRank(int block, int numRanks)
{
   for(int i=0; i <numRanks; i++)
     {
        if(startInd[i] <= block && block <= endInd[i])
           return i;
      }

   return -1;
} 
//------------------------------------------------------
int main(int argc, char **argv)
{

    int numRanks, myRank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

     MPI_Barrier(MPI_COMM_WORLD); //timing
    ofstream dbg;
    ofstream timDet;

    char dbgFile[1024];
    sprintf(dbgFile, "rank%d.txt",myRank);
    //if (myRank == 0)
        dbg.open(dbgFile);
    //else
      //  dbg.open("rank1.txt");
       

    options::Options opts;
    //opts.numDIYBlocks = 1;
    opts.Parse(argc, argv);
    string fold = "";
    if(opts.foldername.size () != 0)
    fold = opts.foldername + "/";

    char timingDet[1024];
    sprintf(timingDet, "%stimingDet%d.txt",fold.c_str(),myRank);
    timDet.open(timingDet);

    char str[1024];
    sprintf(str,"%sPOD%d_%d_MaxStep%d_nSeeds%d_vtkm:%d",fold.c_str(),myRank,numRanks,opts.maxAdvectSteps,opts.numSeeds,opts.vtkmON);
    if(myRank ==0) cerr<<" filePath "<<str<<" fold "<<fold<<"\n";

    TimingsManager::Initialize(str);
    visitTimer->WithholdOutput(true);
    visitTimer->Enable();

     MPI_Barrier(MPI_COMM_WORLD); //timing
    int t1 = visitTimer->StartTimer();
    int t2 = visitTimer->StartTimer();
    box.min[0] = box.min[1] = box.min[2] = 0;
    box.max[0] = opts.fileSize[0];
    box.max[1] = opts.fileSize[1];
    box.max[2] = opts.fileSize[2];
    if(opts.resample && 1) 
    {
       box.max[0] = opts.resampleSize[0];
       box.max[1] = opts.resampleSize[1];
       box.max[2] = opts.resampleSize[2];


    }
    blockBooking.resize(opts.numDataBlocks);
    blockTable currentBox;
    currentBox.blockID = -1;
     MPI_Barrier(MPI_COMM_WORLD); //timing

    visitTimer->StopTimer(t2, "Init");
   
    int t3 = visitTimer->StartTimer();
    DivideDataToBlocks(opts.numDataBlocks,box,opts.blockDivisions);
    MPI_Barrier(MPI_COMM_WORLD); //timing
    visitTimer->StopTimer(t3, "DivideDataToBlocks");

    int t4 = visitTimer->StartTimer();
    //avtParICAlgorithm a(MPI_COMM_WORLD);
    //int msgSize = 2; //(sender, msg)
    //int numParticleRecvs = std::min(32, numRanks-1);
    //int numMsgRecvs = std::min(32, numRanks-1);
    //a.InitializeBuffers(msgSize, numMsgRecvs, numParticleRecvs);
 
    int numParticles = opts.numSeeds;
    int numSeedsPerProc = opts.numSeeds/numRanks;
    int numOfReminParticles = opts.numSeeds - (numSeedsPerProc*numRanks);
    if(myRank < numOfReminParticles) numSeedsPerProc = numSeedsPerProc + 1;
    if(opts.numSeeds < numRanks && myRank < opts.numSeeds) numSeedsPerProc = 1;
    
     visitTimer->StopTimer(t4, "Setup");

  
    int numBlocksPerMPI = opts.numDataBlocks/numRanks;
    int numOfReminBlocks = opts.numDataBlocks - (numBlocksPerMPI*numRanks);
    int* numBlocksForMPI = new int[numRanks];
    for(int i=0; i < numRanks; i++)
    {

       numBlocksForMPI[i] = numBlocksPerMPI;
       if(i < numOfReminBlocks) numBlocksForMPI[i] = numBlocksForMPI[i] + 1;
       if(opts.numDataBlocks < numRanks && i < opts.numDataBlocks) numBlocksForMPI[i] = 1;

    }
    //vector<int> ownBlocks; 
    startInd = new int[numRanks];
    endInd = new int[numRanks];
    startInd[0]=0;
    endInd[0]= startInd[0] + numBlocksForMPI[0] - 1;
    for(int i=1; i < numRanks; i++)
    {
        startInd[i] = startInd[i-1] + numBlocksForMPI[i-1];
        endInd[i] = startInd[i] + numBlocksForMPI[i] - 1;

    }
   
    //for(int i=0; i < numRanks; i++)
      //  cerr<<"Rank "<<i<<" start "<<start[i]<<" end "<<end[i]<<"\n";
     ////////// 

    #if 0
    int* start = new int[numRanks];
    int* end = new int[numRanks];
    start[0]=0;
    for(int i=0; i < numRanks; i++)
      {
        
      }
    #endif
     ////////// 
     if(myRank == 0)
     {

        cerr<<"****** Data Divison: \n";
        for(int i=0; i< blockBooking.size(); i++)
          printBlock(blockBooking[i].blockID);

        cerr<<" *********************** \n\n";

     }
 
    int t5 = visitTimer->StartTimer();
    AddSeeds(myRank,numRanks,&opts,numSeedsPerProc,startInd[myRank],endInd[myRank]);
    MPI_Barrier(MPI_COMM_WORLD); //timing
    visitTimer->StopTimer(t5, "AddSeeds");

   //int tData = visitTimer->StartTimer();
   //GetData(myRank,opts,blockBooking[myRank].BlockBox,blockBooking[myRank].blockID,start[myRank],end[myRank]);
    //MPI_Barrier(MPI_COMM_WORLD); //timing
    //visitTimer->StopTimer(tData, "GetData");

    int t6 = visitTimer->StartTimer();
    bool done = false;
    int numTerminated = 0;
    numSeedsPerProc = active.size(); //TEMP
    int numActive = numParticles;
 
    int numOut = 0;
    int numTerm = 0;
     int numComm = 0; 
    //GetData(myRank,opts,box);
    
    //cerr<<" Rank "<<myRank<<" numParticles = "<<numParticles<<" numSeedsPerProc "<<numSeedsPerProc<<" active size "<<active.size()<<"\n";
    Integrator *integrator; 


/*

    gettimeofday(&startTime, 0);
    //int tInt = visitTimer->StartTimer();
    if (integrator == NULL)
    {

      if(!opts.vtkmON)
            {
                integrator = new RK4Integrator(opts.advectionStep, opts.recordParticlePaths, blockBooking[myRank].BlockBox, &vel,myRank);
            //    integrator = new RK4Integrator(opts.advectionStep, opts.recordParticlePaths, box, &vel,myRank);
               //cerr<<"RK4 \n";

            }
             else
                 integrator = new VTKmIntegrator(opts.advectionStep,numSeedsPerProc,opts.recordParticlePaths, blockBooking[myRank].BlockBox, &vel, active, myRank);
                 //integrator = new VTKmIntegrator(opts.advectionStep,numSeedsPerProc,opts.recordParticlePaths, box, &vel, active, myRank);

           integrator->Init();



    }
    MPI_Barrier(MPI_COMM_WORLD); //timing
    //visitTimer->StopTimer(tInt, "Init Integrator");
    gettimeofday(&endTime, 0);
        dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
        elapsed_secsIntegInit += dt;
*/
int tread = visitTimer->StartTimer();
Particle p;
Particle tp;

//while (!done )
//{

  tp = active.front();

  int isLoad = isLoaded(tp.blockID);   
  if(currentBox.blockID != tp.blockID && isLoad == -1 && startInd[myRank] <= tp.blockID && tp.blockID <= endInd[myRank])
    {

      currentBox.blockID = tp.blockID;
      currentBox.BlockBox = blockBooking[tp.blockID].BlockBox;
      GetData(myRank,opts,currentBox.BlockBox,currentBox.blockID);
      if( opts.numLoadedBlocks <= LoadedBounds.size())
                    LoadedBounds.pop_back();
           DataBlock d;
           d.blockID = tp.blockID;
           d.bounds = currentBox.BlockBox;// bounds;
           d.vel = vel;
           LoadedBounds.push_front(d);

           cerr<<"Loaded data currentBox.blockID = "<<currentBox.blockID<<"\n";

    }
    else if(currentBox.blockID != tp.blockID && isLoad != -1 )
        {
          currentBox.blockID = tp.blockID;
          currentBox.BlockBox = blockBooking[tp.blockID].BlockBox;
          vel = LoadedBounds[isLoad].vel;
          DataBlock d;
          d.blockID = tp.blockID;
          d.bounds = currentBox.BlockBox; //bounds;
          d.vel = vel;
         }
        else
          {
            currentBox.BlockBox = blockBooking[tp.blockID].BlockBox;
            //currentBox.BlockBox = blockBooking[isLoad].BlockBox;
            currentBox.blockID = tp.blockID;
            vel = LoadedBounds[isLoad].vel;


           }

 visitTimer->StopTimer(tread, "Read Time");    
       if(currentBox.blockID == tp.blockID)
         {
            if(!opts.vtkmON)
            {
                integrator = new RK4Integrator(opts.advectionStep, opts.recordParticlePaths, currentBox.BlockBox, &vel,myRank);
            //    integrator = new RK4Integrator(opts.advectionStep, opts.recordParticlePaths, box, &vel,myRank);
               //cerr<<"RK4 \n";

            }
             else
//                 integrator = new VTKmIntegrator(opts.advectionStep,numSeedsPerProc,opts.recordParticlePaths, currentBox.BlockBox, &vel, vtkmVec, active, myRank);
                    //Abhishek            
                    integrator = new VTKmIntegrator(opts.advectionStep,numSeedsPerProc,opts.recordParticlePaths, blockBooking[myRank].BlockBox, &vel, vtkmVec, vtkmVec2, active, myRank,nVecs);

                 //integrator = new VTKmIntegrator(opts.advectionStep,numSeedsPerProc,opts.recordParticlePaths, box, &vel, active, myRank);

           integrator->Init();

          //cerr<<"Created integrator\n";

         }
int tLoop = visitTimer->StartTimer();
   if(!opts.vtkmON)
       {

          while(!active.empty())
          {
             p = active.front();
             active.pop_front();
             int stepsTaken = integrator->Push(p, opts);
          
              if (p.nSteps >= opts.maxAdvectSteps ) //|| !p.inside(domain) )
              {
                numTerminated++;
                terminated.push_back(p);             
               }
               else if( p.coords[0] < box.min[0] || box.max[0]-1 <= p.coords[0] ||
                     p.coords[1] < box.min[1] || box.max[1]-1 <= p.coords[1] ||
                     p.coords[2] < box.min[2] || box.max[2]-1 <= p.coords[2])
                 {
                     numTerminated++;
                     terminated.push_back(p);
                 }
               else if(integrator->status == Integrator::OK)
            {
                if(blockBooking[myRank].BlockBox.min[0] <= p.coords[0] && p.coords[0] < blockBooking[myRank].BlockBox.max[0] &&
              blockBooking[myRank].BlockBox.min[1] <= p.coords[1] && p.coords[1] < blockBooking[myRank].BlockBox.max[1] &&
              blockBooking[myRank].BlockBox.min[2] <= p.coords[2] && p.coords[2] < blockBooking[myRank].BlockBox.max[2])
                active.push_back(p);
               else
                communicate.push_back(p); 
            }
           
              

          }//active ! empty
       }
    else if(opts.vtkmON && !active.empty())
    {

       cerr<<"before push \n";
       Particle tempP;
       gettimeofday(&startTime, 0);
       int stepsTaken = integrator->Push(tempP, opts);
       gettimeofday(&endTime, 0);
       dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
       elapsed_secsPush += dt;
       elapsed_secsAdvection += dt;

       active.clear();
       gettimeofday(&startTime, 0);
       cerr<<"after push\n";
       while (!integrator->outParticles.empty())
       {
 
          p = integrator->outParticles.front();
	  integrator->outParticles.pop_front();
                     dbg<<"Steps "<<stepsTaken<<" P "<<p.coords[0]<<", "<<p.coords[1]<<", "<<p.coords[2]<<" integrator->status "<<integrator->status<<" p.nSteps "<<p.nSteps<<" ";

         if(p.nSteps >= opts.maxAdvectSteps )
           {
             numTerminated++;
             terminated.push_back(p);
              numTerm++;
              dbg<<"TERM nSteps\n";


           }
          else if(CheckBit(p.status,STATUS_TERMINATED))
                   { numTerminated++;
                        terminated.push_back(p);
                        numTerm++;
                        dbg<<"TERM Status\n";
                      }
          else if(CheckBit(p.status,EXITED_SPATIAL_BOUNDARY)  || CheckBit(p.status,EXITED_TEMPORAL_BOUNDARY))
          {
           if( p.coords[0] < box.min[0] || box.max[0]-1 <= p.coords[0] ||
                           p.coords[1] < box.min[1] || box.max[1]-1 <= p.coords[1] ||
                           p.coords[2] < box.min[2] || box.max[2]-1 <= p.coords[2])
            {

               numTerminated++;
                              terminated.push_back(p);
                              numTerm++;
                        dbg<<"TERM Out of Data\n";

            }
             else
              {
                   int b = GetBlockID(p);
                        // To solve 
                        if(b == currentBox.blockID)
                         {
                            if(p.coords[0] < box.max[0])
                            p.coords[0] = p.coords[0] +1;
                            if(p.coords[1] < box.max[1])
                            p.coords[1] = p.coords[1] +1;
                            if(p.coords[2] < box.max[2])
                            p.coords[2] = p.coords[2] +1;

                            b = GetBlockID(p);
                         }
                    p.blockID = b;
                    //dbg<<"startInd "<<startInd[myRank]<<" endInd "<<endInd[myRank]<<" p.blockID "<<p.blockID<<"\n";
                   if(startInd[myRank] <= p.blockID && p.blockID <= endInd[myRank])
                    {
                      active.push_back(p);
                      dbg<<"active.push_back\n";
                    }
                   else   
                   {
                    communicate.push_back(p);
                    dbg<<"COMM to "<<p.blockID<<"\n";
                   }

               }
          }
            else if(CheckBit(p.status,STATUS_OK))
            {
if(blockBooking[startInd[myRank]].BlockBox.min[0] <= p.coords[0] && p.coords[0] < blockBooking[endInd[myRank]].BlockBox.max[0] &&
              blockBooking[startInd[myRank]].BlockBox.min[1] <= p.coords[1] && p.coords[1] < blockBooking[endInd[myRank]].BlockBox.max[1] &&
              blockBooking[startInd[myRank]].BlockBox.min[2] <= p.coords[2] && p.coords[2] < blockBooking[endInd[myRank]].BlockBox.max[2])
                {
                 active.push_back(p);
                   dbg<<"OK \n";
                 }
	      /* else
               {
                if(startInd[myRank] <= b && endInd[myRank] <= b)
                    active.push_back(p);
                else
                communicate.push_back(p);
                   dbg<<"COMM \n";
               }*/

             }


       }//while outParticle 

    gettimeofday(&endTime, 0);
          dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
          elapsed_secsResult += dt;
          elapsed_secsAdvection += dt; 
    }//vtkmon

   if(numActive < 0 ) numActive = 0;
   cerr<<"numActive "<<numActive<<"\n\n";
  
/* 
  gettimeofday(&startTime, 0); 
   if (numTerminated > 0)
        {
           // cerr<<"numTerminated "<<numTerminated<<"\n";
            vector<int> msg(2);
            msg[0] = MSG_TERMINATED;
            msg[1] = numTerminated;
            a.SendAllMsg(msg);
            dbg<<" Msg TERMINATE: "<<numTerminated<<endl;

            numActive -= numTerminated;
            if(numActive < 0 ) numActive = 0;
            numTerminated = 0;

            //cerr<<"numActive "<<numActive<<"\n";
        }
   // Send Particles in communicate 
    //numComm += communicate.size();

    int target;
     list <Particle> sendingParticles;
     for(int i=0; i< communicate.size(); i++)
     {
       target = ownerRank(communicate[i].blockID,numRanks); //GetTarget(communicate[i].blockID, opts.numDataBlocks, numRanks);
       cerr<<"target "<<target<<"\n";
       if(target < 0 || numRanks <= target)
       {
        cerr<<"**** Warning Particle "<<communicate[i]<<" target is out of scope \n";
        numTerminated++;
        terminated.push_back(communicate[i]);
       }
      else
       {   numComm++;
          sendingParticles.push_back(communicate[i]);
          a.SendICs(target,sendingParticles);
           sendingParticles.clear();
       }

     }
      communicate.clear();

a.CheckPendingSendRequests();


      bool block = active.empty() && numActive > 0;
      dbg<<"***** active.empty() "<<active.empty()<<" numActive "<<numActive<<" communicate "<<sendingParticles.size()<<endl;
      vector<MsgCommData> M;
      list<ParticleCommData<Particle>> P;
        dbg<<"Check messages....("<<block<<") ";
	 if (a.RecvAny(&M, &P, NULL, block))
        {
         for (int i = 0; i < M.size(); i++)
            {
               int sender = M[i].rank;
                vector<int> &msg = M[i].message;
                if (msg[0] == MSG_TERMINATED)
                {
                    numActive -= msg[1];
                    dbg<<"   "<<i<<" rcv term "<<msg[1]<<" from "<<sender<<" numActive "<<numActive<<endl;
                }
              else
                    cout<<"Unknown message type: "<<msg[0]<<endl;
            }// loop over MSGs
           if (M.size() > 0)
            dbg<<myRank<<": MSG Size: "<<M.size()<<endl;
	
	  M.resize(0);
          //if(myRank == 0)
          dbg<<" Rank "<<myRank<<" numActive "<<numActive<<" active.size "<<active.size()<<"\n";
          // If I got any particles then loop over the particles
            // and push them to my active list
          dbg<<" P size "<<P.size()<<"\n";
          for (auto pi = P.begin(); pi != P.end(); pi++)
            {
                cout<<"myRank "<<myRank<<" ---> Recv: "<<pi->p.id<<endl;
                cerr<<"    recv particle: "<<pi->p.id<<endl;
                active.push_back(pi->p);
            }
            if (P.size() > 0)
                cout<<myRank<<": PAR*******: "<<P.size()<<endl;
            P.resize(0);
        }// Recv 


    gettimeofday(&endTime, 0);
          dt = float(endTime.tv_sec - startTime.tv_sec) + float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
          elapsed_secsComm += dt;
 dbg<<"\n numActive --------- "<<numActive<<"\n";
   if (numActive == 0)
       done = true;
*/
//}//while !done


 MPI_Barrier(MPI_COMM_WORLD); //timing
    visitTimer->StopTimer(tLoop, "Adection & Communication");    
    cerr<<"numComm "<<numComm<<"\n";
    
    timDet<<"VTK reader time: "<<elapsed_secsVTKRead<<"\n";
    timDet<<"Extracting vector feild time: "<<elapsed_secsExtrVec<<"\n";
    timDet<<"IO time: "<<elapsed_secsIO<<"\n";
    timDet<<"Integrator Init time: "<<elapsed_secsIntegInit<<"\n";
    timDet<<"Integrator Setup time: "<<elapsed_secsIntegSteup<<"\n";
    timDet<<"Calling Advection (Push) time: "<<elapsed_secsPush<<"\n";
    timDet<<"Processing Results from Advection "<<elapsed_secsResult<<"\n";
    timDet<<"Advection time: "<<elapsed_secsAdvection<<"\n";
    timDet<<"Comm time: "<<elapsed_secsComm<<"\n";

     MPI_Barrier(MPI_COMM_WORLD); //timing
    //cout<<myRank<<": I have "<<active.size()<<endl;    

    //if(myRank == 0) cerr<<"numRanks "<<numRanks<<"\n";
    visitTimer->StopTimer(t1, "Total");
   TimingsManager::Finalize();

  MPI_Finalize();
};
