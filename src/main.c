#include <mygraph.h>
#include <unistd.h> //sleep
#include <time.h> //nanosleep
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

triangleSurface mySurf = {.v0 = {0,0,0},
                          .v1 = {50,0,0},
                          .v2 = {0,50,0},
                          .reflectionIndex = 0,
                          .lightSource =0,
                          .hue = 0, /*should be red*/
                          .sat = .5};

void init()
{
  cam.pose = gsl_vector_alloc(3);
  gsl_vector_set(cam.pose, 0, -.5*LENGTH);
  gsl_vector_set(cam.pose, 1, 0*.3*LENGTH);
  gsl_vector_set(cam.pose, 2, .6 * LENGTH);
  cam.pan = 0;
  cam.tilt = .2;
  cam.zoom = 2;
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
      if (colState == INTERSECT) printf("HIT\n");
    }
  }

  //match pixels to resulting array, draw to screen rays which have at least one refelction
  //to show that rays are being generated and the pinhole is working
  
  //iterate backwards? to rotate rays by 180 to make screen make sense

  for (int r =119; r>=0; r--)
  {
    for (int c=119; c>=0; c--)
    {
      if (res[r + RAYROW*c].ray.reflections >0)
      {
         int xdraw,ydraw;
         xdraw = ((double)(119 - c) / 119.0) * xdim;
         ydraw = ((double)(119 - r) / 119.0) * ydim;
         myfilledcircle(4, xdraw, ydraw, 2);
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