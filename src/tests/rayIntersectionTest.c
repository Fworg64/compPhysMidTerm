#include "tracingRay.h"
#include "triangleSurface.h"

#include <stdio.h>

tracingRay ray;
intersectionResults res;
triangleSurface surf;
collisionState col;

int main(){
  printf("starting test\n");

  surf.v0[0]=50;
  surf.v0[1]=50;
  surf.v0[2]=0;
  surf.v1[0]=50;
  surf.v1[1]=-25;
  surf.v1[2]=25;
  surf.v2[0]=50;
  surf.v2[1]=-25;
  surf.v2[2]=-25;

  ray.origin[0]=0;
  ray.origin[1]=0;
  ray.origin[2]=0;
  ray.direction[0]=1;
  ray.direction[1]=0;
  ray.direction[2]=0;
  ray.reflections=0;
  
  printf("going into routine\n");
  col = intersect(&surf, &ray, &res);
  printf("resulting ray: origin < %.4f, %.4f, %.4f> direction < %.4f, %.4f, %.4f>\n", res.ray.origin[0], res.ray.origin[1], res.ray.origin[2], res.ray.direction[0], res.ray.direction[1], res.ray.direction[2]);
  if (col== MISS){
    printf("Miss");
  }
  else if (col==ABSORBED){
    printf("Absorbed");
  }
  else{
    printf("Intersect");
  }
}
