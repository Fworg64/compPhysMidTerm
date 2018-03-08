#include "tracingRay.h"
#include "camera.h"

#include "renderer.h"
#include <stdio.h>

#define RAYCOL 120
#define RAYROW 120

#define LENGTH 50

tracingRay rays[120][120];

camera cam;

int main()
{

  cam.pose = gsl_vector_alloc(3);
  gsl_vector_set(cam.pose, 0, -.5*LENGTH);
  gsl_vector_set(cam.pose, 1, .3*LENGTH);
  gsl_vector_set(cam.pose, 2, .6 * LENGTH);
  
  cam.pan = 0;
  cam.tilt = .4;
  cam.zoom = 200;
  cam.roll =0;



  generateRaysFromCamera(cam, rays, 120, 120);
  
  for (int row =0; row < RAYROW; row++)
  {
    for (int col =0; col < RAYCOL; col++)
    {
      printf("ray: ");
    }
  }

}