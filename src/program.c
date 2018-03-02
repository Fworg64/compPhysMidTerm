#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>


#define NUM_PARTICLES 240
#define DIM 3
#define SAMPLELENGTH 400000

#define RUNNINGAVGLENGTH 5000

double LENGTH=50; //side length of box

double masses[NUM_PARTICLES];
double x[NUM_PARTICLES][DIM];
double v[NUM_PARTICLES][DIM];
double tplot[SAMPLELENGTH] = {0};

double runningAvgPlot[SAMPLELENGTH] = {0};
int runningAvgPlotIndex =0;
double runningAvgBuff[RUNNINGAVGLENGTH] = {0};
double runningAvgSum =0;
int runningAvgIndex =0;

double dt =.001;
char buffer[500] = {0};
int firstrun =1;
double setTempVal = 1.5;
double setDensityVal = .05;

double nonIdealPressure;

int points = SAMPLELENGTH;

struct camera_S 
{
  gsl_vector * pose;
  double pan;
  double tilt;
  double roll;
  double zoom;
} mycam;

struct renderman_S
{
  gsl_vector * camU;
  gsl_vector * camV;
  gsl_vector * camW;
  gsl_matrix * rot;
  gsl_matrix * trans;
  gsl_matrix * renderCoord;
  gsl_matrix * worldCoordHolder;
  gsl_matrix * coordinateTransform;
  gsl_matrix * boxpoints;
  gsl_matrix * renderBoxpoints;
} renderman;

void calculateForce(double pos[NUM_PARTICLES][DIM], double vel[NUM_PARTICLES][DIM], 
                    double force[NUM_PARTICLES][DIM])
{
  nonIdealPressure =0;
  memset(force,0,NUM_PARTICLES*DIM*sizeof(double));
  for (int n=0; n<NUM_PARTICLES; n++)
  {
    for (int m=n+1; m<NUM_PARTICLES; m++)
    {
      double dr[DIM],dR=0;
      for (int d=0; d<DIM; d++)
      {
	    dr[d]=x[m][d]-(x[n][d]-LENGTH);
	    double ddr;
	    ddr=x[m][d]-(x[n][d]);
	    if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr;
	    ddr=x[m][d]-(x[n][d]+LENGTH);
	    if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr;
	    dR+=dr[d]*dr[d];
      }
      double dR6=dR*dR*dR;
      double dR12=dR6*dR6;
      double Fabs=12/dR12-6/dR6;
      nonIdealPressure += Fabs;
      Fabs/=dR;
      for (int d=0;d<DIM; d++){
	    force[n][d]-=Fabs*dr[d];
	    force[m][d]+=Fabs*dr[d];
      }
    }
  }

  for (int d=0; d<DIM; d++) nonIdealPressure/=LENGTH;

  return;
}

void iterate(double pos[NUM_PARTICLES][DIM], 
             double vel[NUM_PARTICLES][DIM], double dt)
{
  double force[NUM_PARTICLES][DIM];
  calculateForce(pos, vel, force);
  if (firstrun == 1)
  {
	for (int n=0; n<NUM_PARTICLES;n++)
    {
      for (int d=0; d<DIM; d++)
	  {
    	vel[n][d] += .5 * (force[n][d]) / masses[n] * dt;
      }
    }
    firstrun =0;
  }
  else
  {
    for (int n=0; n<NUM_PARTICLES;n++)
    {
  	  for (int d=0; d<DIM; d++)
  	  {
  	    vel[n][d] += (force[n][d]) / masses[n] * dt;
      }
    }
  }
  for (int n=0; n<NUM_PARTICLES;n++)
  {
    for (int d=0; d<DIM; d++)
    {
      pos[n][d] += vel[n][d] * dt;
      while (pos[n][d] < 0) pos[n][d] += LENGTH;
      while (pos[n][d] >= LENGTH) pos[n][d] -= LENGTH;
    }
  }
}

double findTemp()
{
    double temp =0;
    for (int n=0; n<NUM_PARTICLES;n++)
    {
      for (int d=0; d<DIM; d++)
	  {
    	temp += masses[n] * v[n][d]*v[n][d];
      }
    }
    return temp/NUM_PARTICLES/DIM;
}

double findDensity()
{
    double density =0;
    for (int n=0; n<NUM_PARTICLES;n++)
    {
       density += masses[n];
    }
    density = density / (double)(LENGTH*LENGTH);
    return density;
}

double findPressure(double temp, double density)
{
    return temp*density + nonIdealPressure;
}

void init()
{
  int M=pow(NUM_PARTICLES,1./DIM)+1;
  for (int n=0; n<NUM_PARTICLES; n++)
  {
    masses[n]=1;
    for (int d=0; d<DIM; d++)
    {
      int nn=n;
      for (int dd=0; dd<d; dd++) nn/=M;

      x[n][d]=(nn%M)*LENGTH/(double)M;
      if (d==1)
      {
	    if (x[n][0]<LENGTH/2.0) 
        {
          v[n][d]=1;
        }
	    else 
        {
          v[n][d]=-1;
        }
      }
      else 
      {
      v[n][d]=0;
      }
    }
  }
  v[0][0] =1;
  firstrun =1;
}

void initCollision()
{
  int M=pow(NUM_PARTICLES,1./DIM)+1;

  for (int n=0; n<NUM_PARTICLES/2; n++)
  {
    masses[n]=1;
    for (int d=0; d<DIM; d++)
    {
      int nn=n;
      for (int dd=0; dd<d; dd++) nn/=M;

      x[n][d]=(nn%M)*LENGTH/(double)M/4.0 + LENGTH/8.0;
      v[n][d] = LENGTH / 30.0;
    }
  }

  for (int n=NUM_PARTICLES/2; n<NUM_PARTICLES; n++)
  {
    masses[n]=1;
    for (int d=0; d<DIM; d++)
    {
      int nn=n;
      for (int dd=0; dd<d; dd++) nn/=M;

      x[n][d]=(nn%M)*LENGTH/(double)M/4.0 + LENGTH/2.0 + LENGTH/8.0;
      v[n][d] = -LENGTH / 30.0;

    }
  }
  firstrun =1;
}

void initDraw()
{
  mycam.pose = gsl_vector_alloc(3);
  gsl_vector_set(mycam.pose, 0, -.5*LENGTH);
  gsl_vector_set(mycam.pose, 1, .3*LENGTH);
  gsl_vector_set(mycam.pose, 2, .6 * LENGTH);
  
  mycam.pan = 0;
  mycam.tilt = .4;
  mycam.zoom = 200;
  mycam.roll =0;

  renderman.camU = gsl_vector_alloc(3);
  renderman.camV = gsl_vector_alloc(3);
  renderman.camW = gsl_vector_alloc(3);
  renderman.rot = gsl_matrix_alloc(4,4);
  renderman.trans = gsl_matrix_alloc(4,4);
  renderman.boxpoints = gsl_matrix_alloc(4,8);
  renderman.renderBoxpoints = gsl_matrix_alloc(4,8);

  gsl_matrix_set(renderman.rot, 3,3,1);
  for (unsigned int i=0; i<8; i++)
  {
    gsl_matrix_set(renderman.boxpoints, 3, i, 1);
    gsl_matrix_set(renderman.renderBoxpoints, 3, i, 1);
  }
  for (unsigned int i=0; i<4; i++)
       gsl_matrix_set(renderman.trans, i,i,1);

  renderman.renderCoord = gsl_matrix_alloc(NUM_PARTICLES, 4);
  renderman.worldCoordHolder = gsl_matrix_alloc(4, NUM_PARTICLES);
  renderman.coordinateTransform = gsl_matrix_alloc(4,4);
  gsl_matrix_set_all(renderman.worldCoordHolder, 1);
}

void setTemp()
{
   double currTemp = findTemp();
   double factor = sqrt(setTempVal/currTemp);

  for (int n=0; n<NUM_PARTICLES;n++)
  {
    for (int d=0; d<DIM; d++)
    {
       v[n][d] *= factor;
    }
  }
}

void setDensity()
{
  double factor = 1.0 /LENGTH;
  double currDensity = findDensity();
  
  LENGTH = pow(NUM_PARTICLES/setDensityVal, 1.0/DIM);
  factor *= LENGTH;
  for (int n=0; n<NUM_PARTICLES;n++)
  {
    for (int d=0; d<DIM; d++)
    {
       x[n][d] *= factor; //assumes life is in the 1st quadrant
    }
  } 
}


int compare(const void *x1,const void *x2){
  if (((double *) x1)[2]< ((double *)x2)[2]) return 1;
  else return -1;
}

void draw3d(int xdim, int ydim)
{
  //rotate by pan and tilt(rotation around Z/W)
  gsl_vector_set_basis(renderman.camU, 0);
  gsl_vector_set_basis(renderman.camV, 1);
  gsl_vector_set_basis(renderman.camW, 2);
  gsl_blas_drot(renderman.camU, renderman.camV, cos(mycam.roll), sin(mycam.roll));
  gsl_blas_drot(renderman.camU, renderman.camW, cos(mycam.pan), sin(mycam.pan));
  gsl_blas_drot(renderman.camV, renderman.camW, cos(mycam.tilt), sin(mycam.tilt));

  //form rotation and translation matrices
  for (unsigned int i=0; i<3; i++)
  {
    gsl_matrix_set(renderman.trans, i, 3, gsl_vector_get(mycam.pose, i));
    gsl_matrix_set(renderman.rot,   i, 0, gsl_vector_get(renderman.camU, i));
    gsl_matrix_set(renderman.rot,   i, 1, gsl_vector_get(renderman.camV, i));
    gsl_matrix_set(renderman.rot,   i, 2, gsl_vector_get(renderman.camW, i));
  }

  //get box vertices ready
  gsl_matrix_set(renderman.boxpoints, 0, 1, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 1, 2, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 2, 3, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 0, 4, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 1, 4, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 1, 5, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 2, 5, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 2, 6, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 0, 6, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 0, 7, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 1, 7, LENGTH);
  gsl_matrix_set(renderman.boxpoints, 2, 7, LENGTH);

  //printf("\nrotation matrix: \n");
  //for (unsigned int i=0; i<4; i++)
  //printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.rot, i,0),
  //                            gsl_matrix_get(renderman.rot, i,1),
  //                            gsl_matrix_get(renderman.rot, i,2), 
  //                            gsl_matrix_get(renderman.rot, i,3));

  //printf("\ntranslation matrix: \n");
  //for (unsigned int i=0; i<4; i++)
  //printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.trans, i,0),
  //                            gsl_matrix_get(renderman.trans, i,1),
  //                            gsl_matrix_get(renderman.trans, i,2), 
  //                            gsl_matrix_get(renderman.trans, i,3));


  //cast pose of each particle to our matrix type via nasty copy
  for (unsigned int i=0; i<NUM_PARTICLES; i++)
  {
    for (unsigned int d=0; d<DIM;d++)
    {
      gsl_matrix_set(renderman.worldCoordHolder, d, i, x[i][d]);
    }
  }

  //calculate total transformation matrix
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  1, renderman.rot, renderman.trans,
                 0, renderman.coordinateTransform);

  //printf("\ntransform matrix: \n");
  //for (unsigned int i=0; i<4; i++)
  //printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.coordinateTransform, i,0), 
  //                            gsl_matrix_get(renderman.coordinateTransform, i,1),
  //                            gsl_matrix_get(renderman.coordinateTransform, i,2), 
  //                            gsl_matrix_get(renderman.coordinateTransform, i,3));

  //perform total transformation
  gsl_blas_dgemm(CblasTrans, CblasTrans, 1, renderman.worldCoordHolder, 
                 renderman.coordinateTransform,  0, renderman.renderCoord);
  //transform box points
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, renderman.coordinateTransform,
                 renderman.boxpoints, 0, renderman.renderBoxpoints);


  //then project particles onto image plane via
  //ximage = f*xcam/zcam
  //yimage = f*ycam/zcam
  //where f is focal length
  //pick f (and camera location) so entire bounded box is displayed comfortably
  //draw box vertices
  int screenvertices[2][8];
  for (unsigned int i=0; i<8; i++)
  {
    int xdraw = mycam.zoom * gsl_matrix_get(renderman.renderBoxpoints, 0, i) 
                / gsl_matrix_get(renderman.renderBoxpoints, 2, i) + (double)xdim/2.0;
    int ydraw = mycam.zoom * gsl_matrix_get(renderman.renderBoxpoints, 1, i) 
                / gsl_matrix_get(renderman.renderBoxpoints, 2, i) + (double)ydim/2.0;
    myfilledcircle(2 + ((i == 0) ? 0 : 1), xdraw, ydraw, 2);
    screenvertices[0][i] = xdraw;
    screenvertices[1][i] = ydraw;
  }
  //draw bounding lines
  mydrawline(2, screenvertices[0][0], screenvertices[1][0], 
                screenvertices[0][1], screenvertices[1][1]); //red
  mydrawline(3, screenvertices[0][0], screenvertices[1][0], 
                screenvertices[0][2], screenvertices[1][2]); //green
  mydrawline(4, screenvertices[0][0], screenvertices[1][0], 
                screenvertices[0][3], screenvertices[1][3]); //blue

  mydrawline(6, screenvertices[0][4], screenvertices[1][4], 
                screenvertices[0][1], screenvertices[1][1]);
  mydrawline(6, screenvertices[0][4], screenvertices[1][4], 
                screenvertices[0][2], screenvertices[1][2]);
  mydrawline(6, screenvertices[0][4], screenvertices[1][4], 
                screenvertices[0][7], screenvertices[1][7]);

  mydrawline(6, screenvertices[0][5], screenvertices[1][5], 
                screenvertices[0][7], screenvertices[1][7]);
  mydrawline(6, screenvertices[0][5], screenvertices[1][5], 
                screenvertices[0][2], screenvertices[1][2]);
  mydrawline(6, screenvertices[0][5], screenvertices[1][5], 
                screenvertices[0][3], screenvertices[1][3]);

  mydrawline(6, screenvertices[0][6], screenvertices[1][6], 
                screenvertices[0][7], screenvertices[1][7]);
  mydrawline(6, screenvertices[0][6], screenvertices[1][6], 
                screenvertices[0][3], screenvertices[1][3]);
  mydrawline(6, screenvertices[0][6], screenvertices[1][6], 
                screenvertices[0][1], screenvertices[1][1]);

  //sort rendercoord by camera distance
  //assign color
  for (unsigned int i=0; i<NUM_PARTICLES; i++)
  {
    gsl_matrix_set(renderman.renderCoord, i, 3, i%3 +2);
  }
  //sort by depth
  qsort(gsl_matrix_ptr(renderman.renderCoord,0,0), NUM_PARTICLES, 
        4*sizeof(double), &compare);

  //display particles
  for (unsigned int i=0; i<NUM_PARTICLES; i++)
  {
    int xdraw = mycam.zoom * gsl_matrix_get(renderman.renderCoord, i, 0) 
                / gsl_matrix_get(renderman.renderCoord, i, 2) + (double)xdim/2.0;
    int ydraw = mycam.zoom * gsl_matrix_get(renderman.renderCoord, i, 1) 
                / gsl_matrix_get(renderman.renderCoord, i, 2) + (double)ydim/2.0;
    myfilledcircle((int) gsl_matrix_get(renderman.renderCoord, i, 3), xdraw, ydraw, 3);

   /*printf("\nworld coord: \n");
   printf("%f, %f, %f\n", x[i][0],x[i][1],x[i][2]);
   printf("\nworld coord (in holder): \n");
   printf("%f, %f, %f, %f\n", gsl_matrix_get(renderman.worldCoordHolder, 0, i),
                              gsl_matrix_get(renderman.worldCoordHolder, 1, i),
                              gsl_matrix_get(renderman.worldCoordHolder, 2, i),
                              gsl_matrix_get(renderman.worldCoordHolder, 3, i));
   printf("\ncamera coord: \n");
   printf("%f, %f, %f, %f\n", gsl_matrix_get(renderman.renderCoord, i, 0),
                              gsl_matrix_get(renderman.renderCoord, i, 1),
                              gsl_matrix_get(renderman.renderCoord, i, 2),
                              gsl_matrix_get(renderman.renderCoord, i, 3));
   printf("\ndraw coord: \n");
   printf("%d, %d\n", xdraw, ydraw);*/
  }
}

void produceAverages(double mytemp)
{
    runningAvgSum -= runningAvgBuff[runningAvgIndex]; 
    runningAvgBuff[runningAvgIndex] = mytemp;
    runningAvgSum += mytemp;
    if (++runningAvgIndex >= RUNNINGAVGLENGTH) runningAvgIndex =0;

    runningAvgPlot[runningAvgPlotIndex] = runningAvgSum / ((double)RUNNINGAVGLENGTH);
    if (++runningAvgPlotIndex >= SAMPLELENGTH -1) runningAvgPlotIndex =0;
}

int main(){
  struct timespec ts={0,100};
  int cont=0;
  int sstep=0;
  int done=0;
  int repeat=100;

  double mytemp =0;
  unsigned int tempindex =0;


  initDraw();

  init();

  AddFreedraw("Particles",&draw3d);
  StartMenu("Newton",1);
  DefineDouble("dt",&dt);
  DefineDouble("zoom", &(mycam.zoom));
  DefineDouble("pan", &(mycam.pan));
  DefineDouble("tilt", &(mycam.tilt));
  DefineDouble("roll", &(mycam.roll));

  DefineDouble("camX", gsl_vector_ptr(mycam.pose, 0));
  DefineDouble("camY", gsl_vector_ptr(mycam.pose, 1));
  DefineDouble("camZ", gsl_vector_ptr(mycam.pose, 2));
  //DefineDouble("k",&k);
  DefineGraphN_R("Temp vs Time",&tplot[0],&points,NULL);
  SetDefaultColor(2);
  DefineGraphN_R("Running Avg (50 samples)", runningAvgPlot, &points, NULL);

  StartMenu("init menu",0);
  for (int n=0; n<NUM_PARTICLES; n++)
  {
    //DefineDouble("x0", &x0[n][0]);
    //DefineDouble("y0", &x0[n][1]);
    //DefineDouble("vel x0", &v0[n][0]);
    //DefineDouble("vel y0", &v0[n][1]);
    //DefineDouble("charge",&charges[n]);
    DefineDouble("mass"  ,&masses[n]);
  }
  DefineFunction("Do init",&init);
  EndMenu();
  DefineFunction("Collision Init", &initCollision);
  DefineDouble("Curr Temp", &mytemp);
  DefineDouble("Set Temp val", &setTempVal);
  DefineFunction("do Set Temp", &setTemp);
  DefineGraph(curve2d_, "Time vs Temp and Avgs");

  DefineDouble("Set Density val", &setDensityVal);
  DefineFunction("do Set Density", &setDensity);

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
         iterate(x,v,dt);
         mytemp = findTemp();
         tplot[tempindex] =mytemp;
         if (++tempindex >= SAMPLELENGTH-1) {tempindex =0; cont =0;}
         produceAverages(mytemp);
      }
      if (cont) nanosleep(&ts,NULL);
    }
    else sleep(1);
  }
}