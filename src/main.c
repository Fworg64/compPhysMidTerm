#include <mygraph.h>
#include "camera.h"
#include "renderer.h"
#include "tracingRay.h"
#include "triangleSurface.h"

//pinhole projector example
tracingRay rays[120][120];
intersectionResults res[120][120];

camera cam = {.pose = {-20, 0, 10},
                .pan = 0,
                .tilt =0,
                .roll =0,
                .zoom = 1};

triangleSurface mySurf = {.v0 = {0,0,0},
                          .v1 = {5,0,0},
                          .v2 = {0,5,0},
                          .reflectionIndex = 0,
                          .lightSource =0,
                          .hue = .5,
                          .sat = .5};


void draw3d(int xdim, int ydim)
{
  generateRaysFromCamera(cam, rays, 120, 120);
  for (int r =0; r<120; r++)
  {
    for (int c=0; c<120; c++)
    {
      intersect(&mySurf, &(rays[r][c]), &(res[r][c]);
    }
  }
  //match pixels to resulting array, draw to screen rays which have at least one refelction
  //to show that rays are being generated and the pinhole is working
  
  //iterate backwords? to rotate rays by 180 to make screen make sense
  for (int r =119; r>=0; r--)
  {
    for (int c=119; c>=0; c--)
    {
      if (res[r][c].reflections >0)
      {
         int xdraw,ydraw;
         xdraw = ((double)(119 - c) / 119.0) * xdim;
         ydraw = ((double)(119 - r) / 119.0) * ydim;
         myfilledcircle(4, xdraw, ydraw, 2);
      }
    }
  }
}

int main()
{
  struct timespec ts={0,100};
  int cont=0;
  int sstep=0;
  int done=0;
  int repeat=100;

  AddFreedraw("Particles",&draw3d);
  
  StartMenu("Newton",1);
  DefineDouble("dt",&dt);
  DefineDouble("zoom", &(cam.zoom));
  DefineDouble("pan", &(cam.pan));
  DefineDouble("tilt", &(cam.tilt));
  DefineDouble("roll", &(cam.roll));

  DefineDouble("camX", cam.pose[0]);
  DefineDouble("camY", cam.pose[1]);
  DefineDouble("camZ", cam.pose[2]);
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