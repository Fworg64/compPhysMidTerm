#include <mygraph.h>
#include <unistd.h> //sleep
#include <time.h> //nanosleep
#include <string.h> //memcpy
#include "camera.h"
#include "renderer.h"
#include "tracingRay.h"
#include "triangleSurface.h"

#define LENGTH 50
#define RAYCOL 120
#define RAYROW 120

//pinhole projector example
tracingRay rays[RAYCOL*RAYROW];
intersectionResults res[RAYCOL*RAYROW];

camera cam;

triangleSurface mySurf = {.v0 = {05,0,-30},
                          .v1 = {30,0,-50},
                          .v2 = {10,30,-40},
                          .reflectionIndex = 0,
                          .lightSource =0,
                          .hue = 0, /*should be red*/
                          .sat = .5};

triangleSurface myAntiSurf = {.v0 = {10,30,-40},
                          .v1 = {30,0,-50},
                          .v2 = {05,0,-30},
                          .reflectionIndex = 0,
                          .lightSource =0,
                          .hue = 0, /*should be red*/
                          .sat = .5};

triangleSurface ground = {.v0 = {-50,0,0},
                          .v1 = {500,-500,0},
                          .v2 = {500,500,0},
                          .reflectionIndex = 0,
                          .lightSource =0,
                          .hue = 0, /*should be red*/
                          .sat = .5};

void init()
{
  cam.pose = gsl_vector_alloc(3);
  gsl_vector_set(cam.pose, 0, -50);
  gsl_vector_set(cam.pose, 1, 20);
  gsl_vector_set(cam.pose, 2, 40);
  cam.pan = -.5;
  cam.tilt = 0;
  cam.zoom = .125;
  cam.roll =0;
}


void draw3d(int xdim, int ydim)
{
  generateRaysFromCamera(cam, rays, RAYROW, RAYCOL);
  collisionState colState;
  for (int row =0; row<RAYROW; row++)
  {
    for (int col=0; col<RAYCOL; col++)
    {
      colState = intersect(&mySurf, &(rays[row + RAYROW*col]), &(res[row + RAYROW*col]));
      if (colState == MISS) 
      {
        // printf("MISS\n");
         for (int i=0; i<3;i++)
         {
           res[row + RAYROW*col].ray.origin[i] = rays[row + RAYROW*col].origin[i];
           res[row + RAYROW*col].ray.direction[i] = rays[row + RAYROW*col].direction[i];
         }
        res[row + RAYROW*col].ray.reflections =0;

        //bounce the ray off of the ground
        intersectionResults tempRes;
        colState = intersect(&ground, &(res[row + RAYROW*col].ray), &tempRes);
        if (colState == INTERSECT) //missed the triangle, hit the ground
        {
           //memcpy(&(res[row + RAYROW*col]), &tempRes, sizeof(tempRes));
           for (int i=0; i<3;i++)
           {
             res[row + RAYROW*col].ray.origin[i] = tempRes.ray.origin[i];
             res[row + RAYROW*col].ray.direction[i] =  tempRes.ray.direction[i];
           }
           res[row + RAYROW*col].ray.reflections =1;
           //check to see if it would hit the triangle after the ground
           colState = intersect(&myAntiSurf, &(res[row + RAYROW*col].ray), &tempRes);
           if (colState == INTERSECT) //this is a shadow
           {
             res[row + RAYROW*col].ray.hue = 120;
             res[row + RAYROW*col].ray.sat = .2;
             res[row + RAYROW*col].ray.lightness = .2;
             res[row + RAYROW*col].ray.reflections =2;
           }
           else //this is just the ground
           {
             res[row + RAYROW*col].ray.hue = 120;
             res[row + RAYROW*col].ray.sat = .5;
             res[row + RAYROW*col].ray.lightness = .5;
             res[row + RAYROW*col].ray.reflections =0;
           }
        }
        else //must have hit the sky
        {
           res[row + RAYROW*col].ray.hue = 240;
           res[row + RAYROW*col].ray.sat = .5;
           res[row + RAYROW*col].ray.lightness = .5;
           res[row + RAYROW*col].ray.reflections =0;
        }
      }
      else {} //hit the triangle
      //else if it intersected, color should have been set in intersect
    }
  }
  //match pixels to resulting array, draw to screen rays which have at least one refelction
  //to show that rays are being generated and the pinhole is working
  
  //iterate backwards? to rotate rays by 180 to make screen make sense


  for (int r =119; r>=0; r--)
  {
    for (int c=119; c>=0; c--)
    {
      if (res[r + RAYROW*c].ray.reflections ==2)
      {
         int xdraw,ydraw;
         xdraw = ((double)(119 - c) / 119.0) * xdim;
         ydraw = ((double)(119 - r) / 119.0) * ydim;
         myfilledcircle(5, xdraw, ydraw, 10);
      }
      if (res[r + RAYROW*c].ray.reflections ==1)
      {
         int xdraw,ydraw;
         xdraw = ((double)(119 - c) / 119.0) * xdim;
         ydraw = ((double)(119 - r) / 119.0) * ydim;
         myfilledcircle(4, xdraw, ydraw, 2);
      }
      else
      {
         int xdraw,ydraw;
         xdraw = ((double)(119 - c) / 119.0) * xdim;
         ydraw = ((double)(119 - r) / 119.0) * ydim;
         if (res[r + RAYROW*c].ray.direction[0] >0)
         {
           myfilledcircle(3, xdraw, ydraw, 2);
         }
         else
         {
           myfilledcircle(2, xdraw, ydraw, 2);
         }
      }
    }
  }

}

void doRender()
{
  render(cam, res, 120, 120, "output.png");
}

int main()
{
  struct timespec ts={0,100};
  int cont=0;
  int sstep=0;
  int done=0;
  int repeat=100;

  init();

  AddFreedraw("Particles",&draw3d);
  
  StartMenu("Newton",1);
  DefineDouble("zoom", &(cam.zoom));
  DefineDouble("pan", &(cam.pan));
  DefineDouble("tilt", &(cam.tilt));
  DefineDouble("roll", &(cam.roll));

  DefineDouble("camX", gsl_vector_ptr(cam.pose, 0));
  DefineDouble("camY", gsl_vector_ptr(cam.pose, 1));
  DefineDouble("camZ", gsl_vector_ptr(cam.pose, 2));
  DefineFunction("do Render", &doRender);

  DefineGraph(freedraw_,"graph2");
  DefineInt("num steps",&repeat);
  DefineBool("step",&sstep);
  DefineLong("NS slow",&ts.tv_nsec);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();

while (!done){
    Events(1);
    DrawGraphs();
    if (cont||sstep){
      sstep=0;
      for (int i=0; i<repeat; i++) 
      {
         //iterate(x,v,dt);
         //mytemp = findTemp();
         //tplot[tempindex] =mytemp;
         //if (++tempindex >= SAMPLELENGTH-1) {tempindex =0; cont =0;}
         //produceAverages(mytemp);
      }
      if (cont) nanosleep(&ts,NULL);
    }
    else sleep(1);
  }
}
