#ifndef TRIANGLE_SURFACE_H
#define TRIANGLE_SURFACE_H

#include "tracingRay.h"

typedef struct
{
  double v0[3], v1[3], v2[3]; //vertices of triangle in world coord
  double reflectionIndex; // 0 to 1 on how dispersive(0)/reflective(1) surface is
  double lightSource; //how much light the surface produces (0 to 1)
  double hue, sat; //hue and sat are from material

} triangleSurface;


bool intersect(triangleSurface *surf, tracingRay *ray, &intersectionResults res);


#endif