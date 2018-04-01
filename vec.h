#ifndef __VEC_H
#define __VEC_H

//#define DIY_USE_SPDLOG 1

//#include "diy_utils.h"
#include <iostream>
#include <cmath>

class Vec
{
public:
    Vec() {coords[0]=0;coords[1]=0;coords[2]=0;}
    Vec(const float *c) {coords[0]=c[0];coords[1]=c[1];coords[2]=c[2];}
    Vec(const double *c) {coords[0]=c[0];coords[1]=c[1];coords[2]=c[2];}
    Vec(float x, float y, float z) {coords[0]=x;coords[1]=y;coords[2]=z;}
    Vec(double x, double y, double z) {coords[0]=x;coords[1]=y;coords[2]=z;}
    Vec(const Vec &v) {coords[0]=v.coords[0];coords[1]=v.coords[1];coords[2]=v.coords[2];}

    inline Vec inc(const Vec &v, const double &s) const;
    inline Vec inc(const Vec &v, const float &s) const;
    inline void inc(const Vec &v, const double &s);
    inline void inc(const Vec &v, const float &s);

    inline Vec& operator=(const Vec &v);
    inline Vec operator+(const Vec &v) const;
    inline Vec operator-(const Vec &v) const;
    inline Vec operator*(const float &s) const;
    inline Vec operator*(const double &s) const;

    inline void operator+=(const Vec &v);
    inline void operator-=(const Vec &v);
    inline void operator*=(const float &s);
    inline void operator*=(const double &s);

    inline double dot(const Vec &v) const;
    inline Vec cross(const Vec &v) const;
    inline double length() const { return sqrt(length2()); }
    inline double length2() const;

    //inline bool inside(const Bounds &bounds) const;

    double coords[3];
};

/*
std::ostream &operator<<(std::ostream &os, Vec const &v)
{
    os<<"["<<v.coords[0]<<" "<<v.coords[1]<<" "<<v.coords[2]<<"]";
    return os;
}
*/

inline Vec Vec::inc(const Vec &v, const double &s) const
{
    return Vec(coords[0] + v.coords[0]*s,
               coords[1] + v.coords[1]*s,
               coords[2] + v.coords[2]*s);
}

inline void Vec::inc(const Vec &v, const double &s)
{
    coords[0] += v.coords[0]*s;
    coords[1] += v.coords[1]*s;
    coords[2] += v.coords[2]*s;
}
inline void Vec::inc(const Vec &v, const float &s)
{
    coords[0] += v.coords[0]*s;
    coords[1] += v.coords[1]*s;
    coords[2] += v.coords[2]*s;
}

inline Vec& Vec::operator=(const Vec &v)
{
    coords[0] = v.coords[0];
    coords[1] = v.coords[1];
    coords[2] = v.coords[2];
    return *this;
}
inline Vec Vec::operator+(const Vec &v) const
{
    return Vec(coords[0]+v.coords[0], coords[1]+v.coords[1], coords[2]+v.coords[2]);
}
inline Vec Vec::operator-(const Vec &v) const
{
    return Vec(coords[0]-v.coords[0], coords[1]-v.coords[1], coords[2]-v.coords[2]);
}
inline Vec Vec::operator*(const float &s) const
{
    return Vec(coords[0]*s, coords[1]*s, coords[2]*s);
}
inline Vec Vec::operator*(const double &s) const
{
    return Vec(coords[0]*s, coords[1]*s, coords[2]*s);
}

inline void Vec::operator+=(const Vec &v)
{
    coords[0] += v.coords[0];
    coords[1] += v.coords[1];
    coords[2] += v.coords[2];
}
inline void Vec::operator-=(const Vec &v)
{
    coords[0] -= v.coords[0];
    coords[1] -= v.coords[1];
    coords[2] -= v.coords[2];
}
inline void Vec::operator*=(const float &s)
{
    coords[0] *= s;
    coords[1] *= s;
    coords[2] *= s;
}
inline void Vec::operator*=(const double &s)
{
    coords[0] *= s;
    coords[1] *= s;
    coords[2] *= s;
}

inline double Vec::dot(const Vec &v) const
{
    return coords[0]*v.coords[0] + coords[1]*v.coords[1] + coords[2]*v.coords[2];
}
inline Vec Vec::cross(const Vec &v) const
{
    return Vec(coords[1]*v.coords[2] - coords[2]*v.coords[1],
               coords[2]*v.coords[0] - coords[0]*v.coords[2],
               coords[0]*v.coords[1] - coords[1]*v.coords[0]);
}
inline double Vec::length2() const
{
    return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
}

/*
inline bool Vec::inside(const Bounds &bounds) const
{
    if (coords[0] < bounds.min[0] || coords[0] >= bounds.max[0]) return false;
    if (coords[1] < bounds.min[1] || coords[1] >= bounds.max[1]) return false;
    if (coords[2] < bounds.min[2] || coords[2] >= bounds.max[2]) return false;
    return true;
}
*/
#endif //__VEC_H
