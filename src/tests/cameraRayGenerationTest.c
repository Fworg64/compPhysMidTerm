#include "tracingRay.h"
#include "camera.h"

#include "renderer.h"
#include <stdio.h>

#define RAYCOL 120
#define RAYROW 120

#define LENGTH 50

tracingRay rays[120*120];

camera cam;

int main()
{
  printf("starting test\n");
  cam.pose = gsl_vector_alloc(3);
  gsl_vector_set(cam.pose, 0, -.5*LENGTH);
  gsl_vector_set(cam.pose, 1, .3*LENGTH);
  gsl_vector_set(cam.pose, 2, .6 * LENGTH);
  
  cam.pan = 0;
  cam.tilt = 0;
  cam.zoom = .2;
  cam.roll =0;


  printf("going into routine\n");
  generateRaysFromCamera(cam, rays, 120, 120);
  
  for (int row =0; row < RAYROW; row++)
  {
    for (int col =0; col < RAYCOL; col++)
    {
      printf("ray: origin < %.4f, %.4f, %.4f> direction < %.4f, %.4f, %.4f>\n",
                                      rays[row + RAYROW*col].origin[0], 
                                      rays[row + RAYROW*col].origin[1],
                                      rays[row + RAYROW*col].origin[2],
                                      rays[row + RAYROW*col].direction[0], 
                                      rays[row + RAYROW*col].direction[1],
                                      rays[row + RAYROW*col].direction[2]);
    }
  }

}