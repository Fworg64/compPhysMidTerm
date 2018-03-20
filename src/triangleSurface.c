#include "triangleSurface.h"
#include "tracingRay.h"
#include <math.h>
#include <gsl/gsl_rng.h>

double dotProduct(double E1[3], double E2[3]){
  double result=0;
  for (int i=0; i<3; i++){
    result +=E1[i]*E2[i];
  }
  return result;
}
void crossProduct(double a[3], double b[3], double * result)
{
  //find if ray and surface intersect
  //find point at which intersection takes place

  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

collisionState intersect(triangleSurface *surf, tracingRay *ray, intersectionResults *res)
{
  //find if the ray intersects the plane
  //find if ray and surface intersect
  double edge0[3],edge1[3],edge2[3];
  double N[3];
  double t;
  double p0l0[3];
  //printf("Print from function \n");
  for(int i=0; i<3; i++){
    edge0[i]=(surf->v0[i])-(surf->v1[i]);
    edge1[i]=(surf->v1[i])-(surf->v2[i]);
    edge2[i]=(surf->v2[i])-(surf->v0[i]);
  }
  crossProduct(edge0,edge1,N);
  //printf("N %.4f, %.4f, %.4f\n",N[0],N[1],N[2]);
  double denom = dotProduct(N, ray->direction);
  if (denom > 1e-6){
    double numer[3];
    for (int i=0;i<3;i++){
      p0l0[i]=surf->v0[i] - ray->origin[i];
    }
    t = dotProduct(p0l0,N)/denom;
    //printf("t1 %4.f\n",t);
  }
  else{
    //printf("MISS1");
    return MISS;
  }
  double D = dotProduct(N,surf->v0);
  // printf("D %4.f\n",D);
  t=(dotProduct(N, ray->origin)+D)/(dotProduct(N,ray->direction));
  //printf("t2 %4.f\n",t);
  double P[3],C[3],vp0[3],vp1[3],vp2[3];
  //check edge0
  for (int i=0; i<3;i++){
    P[i] = ray->origin[i]+(t*ray->direction[i]);
    vp0[i]=(P[i] - surf->v0[i]);
  }
  crossProduct(vp0,edge0,C);
  //printf("C %.4f, %.4f, %.4f\n",C[0],C[1],C[2]);
  if (dotProduct(N,C)<0){
    //printf("MISS2\n");
    return MISS;
  }
  //check edge1
  for (int i=0; i<3;i++){
    vp1[i]=(P[i]-surf->v1[i]);
  }
  crossProduct(vp1,edge1,C);
  if (dotProduct(N,C)<0){
    //printf("MISS3\n");
    return MISS;
  }
  //check edge2
  for (int i=0; i<3;i++){
    vp2[i]=(P[i]-surf->v2[i]);
  }
  crossProduct(vp2,edge2,C);
  if (dotProduct(N,C)<0){
    //printf("MISS4\n");
    return MISS;
  }
  res->dist = (P[0])*(P[0])+(P[1])*(P[1])+(P[2])*(P[2]);
  //generate reflected ray based on material properties 
  //random direction from reflectionIndex
  if (0/*insert surface condition/random number thing here*/) {
    
  }
  else if ( 0/*insert surface condition/random number thing here*/){
    return ABSORBED;
  }
  else { //perfect reflection
      double DP,MagN,MagD;
    DP = dotProduct(ray->direction,N);
    for (int i=0;i<3;i++){
      MagN+=N[i]*N[i];
    }
    MagN=sqrt(MagN);
    for (int i =0; i<3; i++){
      N[i]/=MagN;
    }
    for(int i=0;i<3;i++){
      res->ray.direction[i]=(ray->direction[i]-2*(DP)*N[i]);
      MagD+=res->ray.direction[i]*res->ray.direction[i];
    }
    MagD=sqrt(MagD);
    for (int i=0; i<3; i++){
      res->ray.direction[i]/=MagD;
    }
    for (int i=0;i<3;i++){
      res->ray.origin[i]=P[i];
    }
    res->ray.reflections=ray->reflections+1;
  }
  //change hue and saturation if first intersection with surface? average across intersections
  //with materials?

  if (ray->reflections==0)
  {
    res->ray.hue = surf->hue;
    res->ray.sat = .5;
    res->ray.lightness = .5;
  }

  return INTERSECT;
  // get lightness from light sources if intersection?
  


}

