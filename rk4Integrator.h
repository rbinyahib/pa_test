#ifndef __RK4_INTEGRATOR_H
#define __RK4_INTEGRATOR_H

#define DIY_USE_SPDLOG 1

#include "integrator.h"
#include <vtkImageData.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>

//This one doesn't really work...
class RK4Integrator : public Integrator
{
public:
    RK4Integrator(double _h,
                  bool _recordSteps,
                  const Bounds &_bounds,
                  vector<vector<float> > *_vec,
                  int _blockID) :
        Integrator(_h, _recordSteps, _vec), bounds(_bounds),blockID(_blockID),
        locator(NULL), cell(NULL), ds(NULL) {}
    
    virtual ~RK4Integrator()
    {
        if (locator) locator->Delete();
        if (cell) cell->Delete();
        if (ds) ds->Delete();        
    }

    virtual void Init()
    {

       //cerr<<"Init vtk \n";
        ds = vtkImageData::New();
        ds->SetExtent(bounds.min[0], bounds.max[0]-1,
                      bounds.min[1], bounds.max[1]-1,
                      bounds.min[2], bounds.max[2]-1);
        int nP = (*vec)[0].size();
        
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3);
        arr->SetNumberOfTuples(nP);
       // cerr<<"********* np "<<nP<<"\n";
        /*
        FILE * fp1;
        char vtkFile[1024];
        sprintf(vtkFile, "output/out_file_in_rk4_%d.txt",blockID);
        fp1 = fopen( vtkFile, "w" );
        */
        for (int i = 0; i < nP; i++)
        {
            arr->SetComponent(i, 0, (*vec)[0][i]);
            arr->SetComponent(i, 1, (*vec)[1][i]);
            arr->SetComponent(i, 2, (*vec)[2][i]);
           
               //char printsFile[1024];
               //sprintf(printsFile, "%f %f %f\n",(*vec)[0][i],(*vec)[1][i],(*vec)[2][i]);
               //fprintf(fp1, printsFile);
               //fflush(fp1);
        }
        //fflush(fp1);
        //fclose(fp1);
        ds->GetPointData()->SetVectors(arr);
        arr->Delete();
        
        locator = vtkCellLocator::New();
        locator->SetDataSet(ds);
        locator->BuildLocator();
        cell = vtkGenericCell::New();
    }

    // Increases the number of steps and store the location at this step
    // it is like outputlocation array in my code
    virtual int Push(Particle &p, int nSteps)
    {
        int stepsTaken = 0;
        vector<float> steps;
        status = ERROR;

        while (true)
        {
        //cerr<<"In Push \n";
            if (RK4Step(p, h))
            {
                stepsTaken++;
                if (recordSteps)
                {
                    
                    steps.push_back(p.coords[0]);
                    steps.push_back(p.coords[1]);
                    steps.push_back(p.coords[2]);
                    
                    #if 0
                     steps.push_back(p.x);
                     steps.push_back(p.y);
                     steps.push_back(p.z);
                    #endif
                }
                
                if (stepsTaken == nSteps)
                {
                    status = OK;
                    break;
                }
            }
            else
            {
                status = OUT;
                break;
            }
        }

        if (recordSteps)
            RecordSamples(p, steps);
        p.nSteps += stepsTaken;
        return stepsTaken;
    }

private:

    void RecordSamples(Particle &p, vector<float> &steps)
    {
        pair<int,int> key(p.id, p.seq);
        auto it = samples.find(key);
        if (it == samples.end())
            samples[key] = steps;
        else
            it->second.insert(it->second.end(), steps.begin(), steps.end());
    }

    // Calculate next step
    inline void INC(Vec &p, const Vec &v, const double &s) const
    {
        p.coords[0] += v.coords[0] * s;
        p.coords[1] += v.coords[1] * s;
        p.coords[2] += v.coords[2] * s;
   
      //cerr<<"S "<<s<<"\n";
   }    

    // Calculate next step with RK4
    bool RK4Step(Vec &inP, double h)
    {
        Vec k1, k2, k3, k4;
        Vec p = inP;
    
           

        // Pass the point p and get the value of that point at k1

        if (!GetVel(p, k1))
          { for (int i = 0; i < 3; i++)
            inP.coords[i] = p.coords[i];
          //  cerr<<"False 1\n";
            return false;
          }
        INC(p,k1,h/2);
        if (!GetVel(p, k2))
          { for (int i = 0; i < 3; i++)
            inP.coords[i] = p.coords[i];
            //cerr<<"False 2\n";
            return false;
          }
        INC(p,k2,h/2);
        if (!GetVel(p, k3))
          {for (int i = 0; i < 3; i++)
           inP.coords[i] = p.coords[i]; 
            //cerr<<"False 3\n";
           return false;
          }
        INC(p,k3,h);
        if (!GetVel(p, k4))
           {for (int i = 0; i < 3; i++)
             inP.coords[i] = p.coords[i];
            //cerr<<"False 4\n";
             return false;
           }

        //cerr<<" h "<<h<<"\n"; 
        for (int i = 0; i < 3; i++)
         {   inP.coords[i] += h*(1/6.0)*(k1.coords[i] + 2*k2.coords[i] + 2*k3.coords[i] + k4.coords[i]);
     
            //   if(i ==2)
              // cerr<<"k1 "<<k1.coords[i]<<" k2 "<<k2.coords[i]<<" k3 "<<k3.coords[i]<<" k4 "<<k4.coords[i]<<"   Val "<<inP.coords[i]<<"\n";
         }

         return true;
    }
    
    bool GetVel(Vec &p, Vec &vec)
    {
        vtkIdType cellID;
        double pcoords[3], weights[32], tol=1e-4;

        vec.coords[0] = 0;
        vec.coords[1] = 0;
        vec.coords[2] = 0;
       
        //cerr<<"in GetVel \n"; 
        cellID = locator->FindCell(p.coords, tol, cell, pcoords, weights);
        // cerr<<"cellID "<<cellID<<"\n";

        if (cellID < 0)
          {
            vec.coords[0] = p.coords[0];         
            vec.coords[1] = p.coords[1];         
            vec.coords[2] = p.coords[2];         
          //cerr<<" cellID "<<cellID<<" p.coords "<<p.coords[0]<<" "<<p.coords[1]<<" "<<p.coords[2]<<"\n";
           return false;
          }
        //Evaluate the feild at location passed as P and returns it in vec         
        for (int i = 0; i < cell->GetNumberOfPoints(); i++)
        {
            double *v = ds->GetPointData()->GetVectors()->GetTuple(cell->GetPointId(i));
            vec.coords[0] += v[0]*weights[i];
            vec.coords[1] += v[1]*weights[i];
            vec.coords[2] += v[2]*weights[i];
   
         //cerr<<"v[2] "<<v[2]<<" weights[i] "<<weights[i]<<"\n";
      }

        return true;
    }
    
    Bounds bounds;
    //vector<vector<float> > *vec;
    vtkCellLocator *locator;
    vtkGenericCell *cell;
    vtkImageData *ds;
    int blockID;
};

#endif //__RK4_INTEGRATOR_H:q

