#include "triangleSurface.h"

bool intersect(triangleSurface *surf, tracingRay *ray, &intersectionResults res)
{
  //find if ray and surface intersect
   

  //find point at which intersection takes place

  //generate reflected ray based on material properties 
  //random direction from reflectionIndex

  //change hue and saturation if first intersection with surface? average across intersections
  //with materials?

  // get lightness from light sources if intersection?
  


}

void crossProduct(double a[3], double b[3], double result[3])
{
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}
