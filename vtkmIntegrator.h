#ifndef __VTKM_INTEGRATOR_H
#define __VTKM_INTEGRATOR_H

#ifndef VTKM_DEVICE_ADAPTER
//#define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_SERIAL
#define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_TBB
#endif

#define DIY_USE_SPDLOG 1

#include "integrator.h"

#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>
#include <vtkm/worklet/particleadvection/Particles.h>

#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/ParticleAdvectionWorklets.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/worklet/StreamLineUniformGrid.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

//#ifdef __BUILDING_TBB_VERSION__
#include <tbb/task_scheduler_init.h>
//#endif

#include <time.h>
#include <sys/time.h>

#define TIMEINFO timeval


//using namespace vtkm;

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;
int timingPrint = 0;

//Abhishek
struct vectorFeild
{
  //vtkm::Float32 x;
  //vtkm::Float32 y;
  //vtkm::Float32 z;
  vtkm::Vec<vtkm::Float32, 3> vec;
};

double
DiffTime(const struct TIMEINFO &startTime,
                         const struct TIMEINFO &endTime)
{
    double seconds = double(endTime.tv_sec - startTime.tv_sec) +
                     double(endTime.tv_usec - startTime.tv_usec) / 1000000.;

    return seconds;
}



//VTKm Integrator
class VTKmIntegrator : public Integrator
{
public:

    VTKmIntegrator(double _h,
                   int _ns,
                   bool _recordSteps,
                   const Bounds &_bounds,
                   vector<vector<float> > *_vec,
                   std::vector<vtkm::Vec<vtkm::Float32, 3> > _vtkmVec,
 		   vectorFeild* _vtkmVec2,
                   //list<Particle> _seeds,
                   deque<Particle> _seeds,
		   int _blockID,
		   int _numPoints) :
        Integrator(_h, _recordSteps, _vec), stepSize(_h), numSeeds(_ns), bounds(_bounds), particles(_seeds), vtkmVec (_vtkmVec), vtkmVec2(_vtkmVec2), blockID(_blockID), numPoints (_numPoints), tracer(NULL) {}
    virtual ~VTKmIntegrator()
    {
        if (tracer)
            delete tracer;
    }

    virtual void Init()
    {
        tracer = new vtkm::worklet::StreamLineFilterUniformGrid<vtkm::Float32, DeviceAdapter>();
        vtkm::cont::DataSetBuilderUniform dataSetBuilder;
        vtkm::Vec<vtkm::Float32, 3> origin(bounds.min[0],bounds.min[1],bounds.min[2]);
        vtkm::Vec<vtkm::Float32, 3> spacing(1,1,1);
        vtkm::Id3 dims(bounds.max[0]-bounds.min[0],
                       bounds.max[1]-bounds.min[1],
                       bounds.max[2]-bounds.min[2]);
                       
        ds = dataSetBuilder.Create(dims, origin, spacing);
        //cerr<<"bounds.max[0]-bounds.min[0] "<<bounds.max[0]-bounds.min[0]<<"\n"; 
       // cerr<<"vtkm particleAdvection is done \n";
//        cerr<<"bounds "<<bounds.min[0]<<" "<<bounds.min[1]<<" "<<bounds.min[2]<<" "<<bounds.max[0]<<" "<<bounds.max[1]<<" "<<bounds.max[2]<<"\n"; 
    }
   
  //---------------------------------------------------------------
    virtual int Push(Particle& p, int nSteps) 
    //virtual int PushParticles(list<Particle> particles, int nSteps) 
    {
   

        //#ifdef __BUILDING_TBB_VERSION__
        int nT = tbb::task_scheduler_init::default_num_threads(); 
        tbb::task_scheduler_init init(8);
        cerr<<"Running with tbb \n";
        //#endif
        struct TIMEINFO startTime5;
        gettimeofday(&startTime5, 0);
 
        //std::vector<vtkm::Vec<vtkm::Float32, 3> > arr;
        vtkm::Id nP = (*vec)[0].size();
        /*
        FILE * fp1;
        char vtkmFile[1024];
        sprintf(vtkmFile, "output/out_file_in_vtkm_%d.txt",blockID);
        fp1 = fopen( vtkmFile, "w" );
        */
         //cerr<<" ******* Data ******"<<"\n";
        #if 1
        
        std::vector<vtkm::Vec<vtkm::Float32, 3> > arr;
	arr.reserve(nP);
        for (int i = 0; i < nP; i++)
        {   
            vtkm::Float32 x = (*vec)[0][i];
            vtkm::Float32 y = (*vec)[1][i];
            vtkm::Float32 z = (*vec)[2][i];
            arr[i] = vtkm::Vec<vtkm::Float32, 3>(x,y,z);
            //arr.push_back(vtkm::Vec<vtkm::Float32, 3>(x,y,z));
             /*
               char printsFile[1024];
               sprintf(printsFile, "%f %f %f\n",x,y,z);
               fprintf(fp1, printsFile);
               fflush(fp1);
              */
            //if(blockID ==0)
            //cerr<<x<<" "<<y<<" "<<z<<"\n"; 
          //cerr<<"z "<<z<<"\n";          
        }
        #endif
          //fflush(fp1);
          //fclose(fp1);
        //field = vtkm::cont::make_ArrayHandle(arr);
        field = vtkm::cont::make_ArrayHandle(vtkmVec);
        vtkm::Id np = numPoints;
        //Abhishek
        // field2 = vtkm::cont::make_ArrayHandle(vtkmVec2,np);
        //field2 = vtkm::cont::make_ArrayHandle(vtkmVec2,np, vtkm::cont::CopyFlag::On);
        ds.AddField(vtkm::cont::Field("vecData", vtkm::cont::Field::ASSOC_POINTS, field));

        struct TIMEINFO endTime5;
        gettimeofday(&endTime5, 0);
        double t5 = DiffTime(startTime5, endTime5);
       //if(timingPrint)
        cerr << "The time to Build vtkm data is " << t5 << endl;

       // struct TIMEINFO startTime4;
       // gettimeofday(&startTime4, 0);

        typedef vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>> FieldHandle;
        typedef typename FieldHandle::template ExecutionTypes<DeviceAdapter>::PortalConst FieldPortalConstType;
        typedef vtkm::worklet::particleadvection::UniformGridEvaluate<FieldPortalConstType, vtkm::Float32,DeviceAdapter>  RGEvalType;
        //typedef vtkm::worklet::particleadvection::UniformGridEvaluate<FieldPortalConstType, vtkm::Float32>  RGEvalType;
        typedef vtkm::worklet::particleadvection::RK4Integrator<RGEvalType, vtkm::Float32> RK4RGType; 

        //RGEvalType eval(ds);
        //RK4RGType rk4(eval, stepSize, ds);
        RGEvalType eval(ds.GetCoordinateSystem(), ds.GetCellSet(0), field);
        RK4RGType rk4(eval, stepSize);

       //cerr<<"VTK-M stepSize "<<stepSize<<"\n";
      // struct TIMEINFO endTime4;
       // gettimeofday(&endTime4, 0);
        //double t4 = DiffTime(startTime4, endTime4);
      

       //if(timingPrint)
      // cerr << "The time of vtkm definitions is " << t4 << endl;
       // cerr<<"vtkm stepSize "<<stepSize<<"\n";
       
      // cerr<<"VTKM particles size "<<particles.size()<<"\n";  
        //struct TIMEINFO startTime3;
        //gettimeofday(&startTime3, 0);
        std::vector<vtkm::Vec<vtkm::Float32, 3>> seeds;
        std::vector<vtkm::Id> particleSteps;
        vtkm::Id pSteps;
         int count = 0;
         for (auto lit = particles.begin(); lit != particles.end() ; lit++)
	  { 
	        Particle mypt = (*lit); 
                vtkm::Vec<vtkm::Float32, 3> p;
                p[0] = static_cast<vtkm::Float32>(mypt.coords[0]); 
                p[1] = static_cast<vtkm::Float32>(mypt.coords[1]); 
                p[2] = static_cast<vtkm::Float32>(mypt.coords[2] );
                pSteps = static_cast<vtkm::Id>(mypt.nSteps);
                //p.nSteps = static_cast<vtkm::Int32>(mypt.nSteps);  

                

                //cerr<<"VTK-M Particle before "<<p[0]<<", "<<p[1]<<", "<<p[2]<<"\n";              
                seeds.push_back(p);
                particleSteps.push_back(pSteps);
                //count++;
         }
        // cerr<<"seeds array done size "<<seeds.size()<<" \n";

        vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>> seedArray;
         seedArray = vtkm::cont::make_ArrayHandle(seeds);
      
        vtkm::cont::ArrayHandle< vtkm::Id > particleStepsArray; 
        particleStepsArray = vtkm::cont::make_ArrayHandle(particleSteps);

        //struct TIMEINFO endTime3;
        //gettimeofday(&endTime3, 0);
        //double t3 = DiffTime(startTime3, endTime3);
       //if(timingPrint)
       // cerr << "The time of creating seeds vec is " << t3 << endl;

       struct TIMEINFO startTime2;
        gettimeofday(&startTime2, 0);
       vtkm::worklet::ParticleAdvectionResult<vtkm::Float32> res;
       vtkm::worklet::ParticleAdvection particleAdvection;

       //res = particleAdvection.Run(rk4, seedArray, nSteps, DeviceAdapter());
       res = particleAdvection.Run(rk4, seedArray, particleStepsArray, nSteps, DeviceAdapter());
       //res = particleAdvection.Run(rk4, seedArray, field, nSteps, DeviceAdapter());
      // cerr<<"particleAdvection.Run \n";
        struct TIMEINFO endTime2;
        gettimeofday(&endTime2, 0);
        double t2 = DiffTime(startTime2, endTime2);
       //if(timingPrint)
        cerr << "The time of Actual particle vtkm call is " << t2 << endl;
        //cerr<<"VTK-M nSteps "<<nSteps<<"\n"; 
        

        //struct TIMEINFO startTime1;
        //gettimeofday(&startTime1, 0);
        typedef typename vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3> >::PortalConstControl PortalPosition;
        PortalPosition readPortal = res.positions.GetPortalConstControl();

        typedef typename vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl PortalStatus;
        PortalStatus readPortal1 = res.status.GetPortalConstControl();

        typedef typename vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl PortalSteps;
        PortalSteps readPortal2 = res.stepsTaken.GetPortalConstControl();
        
          vtkm::Id index = 0;
      //for (vtkm::Id index = 0; index < readPortal.GetNumberOfValues(); index++)
         for (auto lit = particles.begin(); lit != particles.end() ; lit++)
        {
             Particle mypt = (*lit);
             Particle outpt;
             outpt.id = mypt.id;
             outpt.seq = mypt.seq; 
             outpt.status = readPortal1.Get(index);
             outpt.nSteps = readPortal2.Get(index);
             outpt.coords[0] = readPortal.Get(index)[0];// + 2; 
             outpt.coords[1] = readPortal.Get(index)[1];// + 2; 
             outpt.coords[2] = readPortal.Get(index)[2];// + 2;
             index++; 
             //cerr<<"VTK-M Particle "<<outpt.coords[0]<<", "<<outpt.coords[1]<<", "<<outpt.coords[2]<<" outpt.status "<<outpt.status<<" Block "<<blockID<<"\n"; 
             outParticles.push_back(outpt);

             
       }

/*
        while(!outParticles.empty())
       { 
          Particle p = outParticles.front();
          outParticles.pop_front(); 
          cerr<<"vtkm Particle outParticles "<<p.status<"\n";

       }

*/
    //  struct TIMEINFO endTime1;
      //  gettimeofday(&endTime1, 0);
       // double t1 = DiffTime(startTime1, endTime1);
       //if(timingPrint)
        //cerr << "The time of getting output is " << t1 << endl;
       // cerr<<"outParticles size "<<outParticles.size()<<"\n"; 
        status = OUT;
        return nSteps;
    }

  //---------------------------------------------------------------
   //virtual void setParticles(list<Particle> newParticles)
   virtual void setParticles(deque<Particle> newParticles)
   {
    //  cerr<<"size of new particles "<<newParticles.size()<<"\n";
      particles = newParticles;

   }
 //---------------------------------------------------------------
private:
    vtkm::worklet::StreamLineFilterUniformGrid<vtkm::Float32, DeviceAdapter> *tracer;
    vtkm::cont::DataSet ds;
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3> > field;
    vtkm::cont::ArrayHandle<vectorFeild > field2;
    Bounds bounds;
    double stepSize;
    int numSeeds;
    std::vector<vtkm::Vec<vtkm::Float32, 3> > vtkmVec;
    vectorFeild* vtkmVec2;
    int numPoints;
    //vector<vector<float> > *vec;
    //list<Particle> particles;
    deque<Particle> particles;
    int blockID;
};

#endif //__VTKM_INTEGRATOR_H
