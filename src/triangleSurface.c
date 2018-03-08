#include "triangleSurface.h"
#include "tracingRay.h"
#include "Math.h"
#include "gsl_rng.h"

double dotProduct(double E1[3], double E2[3]){
  double result;
  for (int i=0; i<3; i++){
    result +=E1[i]*E2[i]
  }
  return result;
}
void crossProduct(double a[3], double b[3], double result[3])
{
  //find if ray and surface intersect
  //find point at which intersection takes place

  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

collisionState intersect(triangleSurface *surf, tracingRay *ray, &intersectionResults res)
{
  //find if the ray intersects the plane
  //find if ray and surface intersect
  double edge0[3],edge1[3],edge2[3],N[3],t;
  for(int i=0; i<3, i++){
    edge0[i]=(surf->v0[i])-(surf->v1[i]);
    edge1[i]=(surf->v1[i])-(surf->v2[i]);
    edge2[i]=(surf->v2[i])-(surf->v0[i]);
  }
  crossProduct(edge0,edge1,N);
  double denom = dotProduct(N, ray->direction);
  if (denom > 1e-6){
    double numer[3];
    for (int i=0;i<3;i++){
      p0l0[i]=surf->v0[i]-ray->origin[i];
    }
    t = dotProduct(p0l0,N)/denom;
  }
  else{
    return MISS;
  }
  double D = dotProduct(N,surf->v0);
  t=(dotProduct(N, ray->origin)+D)/(dotProduct(N,ray->direction));
  double P[3],C[3],vp0[3],vp1[3],vp2[3];
  //check edge0
  for (int i=0; i<3;i++){
    P[i] = ray->origin[i]+(t*ray->direction[i]);
    vp0[i]=(P[i]-surface->v0[i]);
  }
  crossProduct(edge0,vp0,C);
  if (dotProduct(N,C)<0){
    return MISS;
  }
  //check edge1
  for (int i=0; i<3;i++){
    vp1[i]=(P[i]-surface->v1[i]);
  }
  crossProduct(edge1,vp1,C);
  if (dotProduct(N,C)<0){
    return MISS;
  }
  //check edge2
  for (int i=0; i<3;i++){
    vp2[i]=(P[i]-surface->v2[i]);
  }
  crossProduct(edge2,vp2,C);
  if (dotProduct(N,C)<0){
    return MISS;
  }
  res->distance = (P[0])(P[0])+(P[1])*(P[1])+(P[2])*(P[2]);
  //generate reflected ray based on material properties 
  //random direction from reflectionIndex
  if (/*insert surface condition/random number thing here*/) {
    
  }
  else if (/*insert surface condition/random number thing here*/){
    return ABSORBED;
  }
  else {
      double DP,MagN,MagDir;
    DP = dotProduct(ray->direction,N);
    for (int i=0;i<3;i++){
      MagN+=N[i]*N[i];
    }
    MagN=sqrt(MagN);
    for (int i =0; i<3; i++){
      N[i]/=MagN;
    }
    for(int i=0;i<3;i++){
      res->ray.direction[i]=ray->direction[i]-2(DP)N[i];
    }
    for (int i=0;i<3;i++){
      res->ray.origin[i]=P[i];
    }
    res->ray.reflections=ray->reflections+1;
  }
  //change hue and saturation if first intersection with surface? average across intersections
  //with materials?

  // get lightness from light sources if intersection?
  


}

