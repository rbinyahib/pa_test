#ifndef __VTK_INTEGRATOR_H
#define __VTK_INTEGRATOR_H

#define DIY_USE_SPDLOG 1

#include "integrator.h"
#include <vtkImageData.h>
#include <vtkStreamTracer.h>

using namespace std;

class VTKIntegrator : public Integrator
{
public:
    VTKIntegrator(double _h,
                  bool _recordSteps,
                  const Bounds &_bounds,
                  vector<vector<float> > *_vec) :
        Integrator(_h, _recordSteps, _vec), bounds(_bounds), tracer(NULL) {}
    virtual ~VTKIntegrator()
    {
        if (tracer) tracer->Delete();
        if (ds) ds->Delete();
    }

    virtual void Init()
    {
        ds = vtkImageData::New();
        ds->SetExtent(bounds.min[0], bounds.max[0]-1,
                      bounds.min[1], bounds.max[1]-1,
                      bounds.min[2], bounds.max[2]-1);
        int nP = (*vec)[0].size();
        //cout<<(*vec)[0].size()<<" "<<(*vec)[1].size()<<" "<<(*vec)[2].size()<<endl;
        
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3);
        arr->SetNumberOfTuples(nP);
        for (int i = 0; i < nP; i++)
        {
            arr->SetComponent(i, 0, (*vec)[0][i]);
            arr->SetComponent(i, 1, (*vec)[1][i]);
            arr->SetComponent(i, 2, (*vec)[2][i]);
        }
        
        ds->GetPointData()->SetVectors(arr);
        arr->Delete();
        //ds->PrintSelf(cout, (vtkIndent)1);

        
        //cout<<"INIT: "<<bounds<<endl;
        tracer = vtkStreamTracer::New();
        tracer->SetMaximumPropagation(std::numeric_limits<double>::max());
        tracer->SetIntegratorTypeToRungeKutta4();
        tracer->SetIntegrationStepUnit(vtkStreamTracer::LENGTH_UNIT);
        tracer->SetIntegrationDirectionToForward();
        tracer->SetInitialIntegrationStep(h);
        tracer->SetComputeVorticity(false);
        tracer->SetInputData(ds);
    }
    
    virtual int Push(Particle &p, int nSteps)
    {
        tracer->SetMaximumNumberOfSteps(nSteps);
        //cout<<"PUSH: "<<p<<endl;
        tracer->SetStartPosition(p.coords);
        tracer->Update();
        vtkPolyData *output = tracer->GetOutput();
        int nPts = output->GetNumberOfPoints();
        
        int stepsTaken = 0;
        status = OK;
        if (nPts > 2)
        {
            stepsTaken = nPts-2;
            output->GetPoint(nPts-1, p.coords);
        
            int stat = output->GetCellData()->GetArray("ReasonForTermination")->GetTuple1(0);
            if (stat == vtkStreamTracer::OUT_OF_STEPS)
                status = OK;
            else if (stat == vtkStreamTracer::OUT_OF_DOMAIN)
                status = OUT;
            else
                status = ERROR;
            
            if (recordSteps)
                RecordSamples(p, output);
        }
        
        p.nSteps += stepsTaken;
        return stepsTaken;
    }

protected:
    void RecordSamples(const Particle &p, vtkPolyData *pd)
    {
        int nPts = pd->GetNumberOfPoints();

        pair<int,int> key(p.id, p.seq);
        auto it = samples.find(key);
        if (it == samples.end())
            samples[key] = vector<float>();

        //Now that we have samples array, allocate the memory.
        it = samples.find(key);
        int n = it->second.size();
        it->second.reserve(n + nPts*3);

        double pt[3];
        for (int i = 0; i < nPts; i++)
        {
            pd->GetPoint(i, pt);
            it->second.push_back(pt[0]);
            it->second.push_back(pt[1]);
            it->second.push_back(pt[2]);
        }
    }

    Bounds bounds;
    //vector<vector<float> > *vec;
    vtkStreamTracer *tracer;
    vtkImageData *ds;
};

#endif //__VTK_INTEGRATOR_H
