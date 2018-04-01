#ifndef __STATS_H
#define __STATS_H

//#define DIY_USE_SPDLOG 1

#include <mpi.h>
#include <math.h>

#include <time.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <unistd.h>
#define TIMEINFO timeval

#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

#include "util.h"

using namespace std;

namespace stats {

class stopWatch
{
public:
    stopWatch(bool _keepHistory=false) : t(0.0f), isRunning(false), keepHistory(_keepHistory) {}
    ~stopWatch() {}

    void start()
    {
        gettimeofday(&startTime, 0);
        isRunning = true;
    }

    float stop()
    {
        gettimeofday(&endTime, 0);
        float dt = float(endTime.tv_sec - startTime.tv_sec) +
            float(endTime.tv_usec - startTime.tv_usec) / 1000000.0f;
        t += dt;
        if (keepHistory)
            history.push_back(dt);
        isRunning = false;
        return t;
    }
    void reset() {t=0.0f; history.resize(0);}
    bool IsRunning() const {return isRunning;}

    bool isRunning, keepHistory;
    float t;
    struct TIMEINFO startTime, endTime;
    vector<float> history;
};

class statisticsDB
{
public:

    template <typename T>
    struct statValue
    {
        statValue() {}
        statValue(vector<T> vals)
        {
            sum = accumulate(vals.begin(), vals.end(), 0);

            //Get min/max info.
            auto res = minmax_element(vals.begin(), vals.end());
            minI = res.first - vals.begin();
            maxI = res.second - vals.end();
            min = (*res.first);
            max = (*res.second);

            //compute mean/median
            int n = vals.size();
            sort(vals.begin(), vals.end());
            mean = (float)sum / (float)n;
            med = vals[vals.size()/2];

            //compute standard deviation
            float x = 0;
            for (int i = 0; i < n; i++)
                x += (vals[i]-mean)*(vals[i]-mean);
            std_dev = sqrt(x/(float)n);
        }

        int minI, maxI;
        T min,max, med, sum;
        float mean, std_dev;
    };

    statisticsDB() : statsComputed(false)
    {
    }
    ~statisticsDB() {}

    //Timers.
    void addTimer(const string &nm, bool keepHistory=false)
    {
        if (timers.find(nm) != timers.end())
            throw nm + " timer already exists!";
        timers[nm] = stopWatch(keepHistory);
        timers[nm].reset();
    }
    void start(const string &nm) {vt(nm); timers[nm].start();}
    float stop(const string &nm) {vt(nm); return timers[nm].stop();}
    float time(const string &nm) {vt(nm); return timers[nm].t;}
    void reset(const string &nm) {vt(nm); timers[nm].reset();}

    //Counters.
    void addCounter(const string &nm)
    {
        if (counters.find(nm) != counters.end())
            throw nm + " counter already exists!";
        counters[nm] = 0;
    }
    void increment(const string &nm) {vc(nm); counters[nm]++;}
    void increment(const string &nm, unsigned long val) {vc(nm); counters[nm]+=val;}
    unsigned long val(const string &nm) {vc(nm); return counters[nm];}

    statValue<float> timerStat(const string &nm) {cs(); return timerStats[nm];}
    statValue<unsigned long> counterStat(const string &nm) {cs(); return counterStats[nm];}
    unsigned long totalVal(const string &nm) {cs(); return counterStats[nm].sum;}

    void cs() {calcStats();}
    void calcStats()
    {
        if (statsComputed)
            return;
        
        int rank, nProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int sz = 0;

        sz = timers.size();
        if (sz > 0)
        {
            vector<float> vals(sz);
            map<string,stopWatch>::iterator it = timers.begin();
            for (int i = 0; it != timers.end(); it++, i++)
                vals[i] = it->second.t;

            it = timers.begin();
            for (int i = 0; i < sz; i++, it++)
            {
                vector<float> res(nProcs, 0.0);
                if (nProcs == 1)
                    res[0] = vals[i];
                else
                {
                    vector<float>tmp(nProcs,0.0);
                    tmp[rank] = vals[i];
                    MPI_Reduce(&tmp[0], &res[0], nProcs, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
                }
                if (rank == 0)
                    timerStats[it->first] = statValue<float>(res);
            }
        }
        
        sz = counters.size();
        if (sz > 0)
        {
            vector<unsigned long> vals(sz);
            map<string,unsigned long>::iterator it = counters.begin();
            for (int i = 0; it != counters.end(); it++, i++)
                vals[i] = it->second;

            it = counters.begin();
            for (int i = 0; i < sz; i++, it++)
            {
                vector<unsigned long> res(nProcs,0);
                if (nProcs == 1)
                    res[0] = vals[i];
                else
                {
                    vector<unsigned long> tmp(nProcs,0);
                    tmp[rank] = vals[i];
                    MPI_Reduce(&tmp[0], &res[0], nProcs, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
                }
                if (rank == 0)
                    counterStats[it->first] = statValue<unsigned long>(res);
            }
        }
        
        statsComputed = true;
    }

private:
    void vt(const string &nm)
    {
        if (timers.find(nm) == timers.end())
            throw nm + " timer not found!";
    }
    void vc(const string &nm)
    {
        if (counters.find(nm) == counters.end())
            throw nm + " counter not found!";
    }
    map<string, stopWatch> timers;
    map<string, unsigned long> counters;

    bool statsComputed;
    map<string, statValue<float> > timerStats;
    map<string, statValue<unsigned long> > counterStats;
};

class funcTimer
{
public:
    funcTimer(const string &nm_, statisticsDB *sdb_) : nm(nm_), sdb(sdb_)
    {
        sdb->start(nm);
    }
    ~funcTimer()
    {
        sdb->stop(nm);
    }
private:
    string nm;
    statisticsDB *sdb;
};
    
}

template <typename T>
ostream &operator<<(ostream &os, const stats::statisticsDB::statValue<T> &s)
{
    return os<<"AVG: "<<s.mean<<" MED: "<<s.med<<" ("<<s.min<<","<<s.max<<":"<<s.std_dev<<")";
}




#endif //__STATS_H
