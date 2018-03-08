#include "renderer.h"
#include <math.h>



void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols)
{

   struct renderman_S renderman;
   renderman.camU = gsl_vector_alloc(3);
   renderman.camV = gsl_vector_alloc(3);
   renderman.camW = gsl_vector_alloc(3);
   renderman.rot = gsl_matrix_alloc(4,4);
   renderman.trans = gsl_matrix_alloc(4,4);
   renderman.boxpoints = gsl_matrix_alloc(4,8);
   renderman.renderBoxpoints = gsl_matrix_alloc(4,8);
   //have camera pose and euler angles, need resolution and size of pixels 
   //and have distance of film from pinhole as zoom

   //parameterize length and height on indices seperated in space by ray density distance
   //let raydistance be .01
   double raydistance = .01;
   double centerx = raydistance * ((double)numRaysRows)/2.0;
   double centery = raydistance * ((double)numRaysCols)/2.0;
   //start at {(centerx - raydistance*numRays/2) , (centery - raydistance*numRays/2)}
   //loop by raydistance, nested dims

   //transform from camera imaging plane to world coord

  //rotate by pan and tilt(rotation around Z/W)
  gsl_vector_set_basis(renderman.camU, 0);
  gsl_vector_set_basis(renderman.camV, 1);
  gsl_vector_set_basis(renderman.camW, 2);
  gsl_blas_drot(renderman.camU, renderman.camV, cos(cam.roll), sin(cam.roll));
  gsl_blas_drot(renderman.camU, renderman.camW, cos(cam.pan), sin(cam.pan));
  gsl_blas_drot(renderman.camV, renderman.camW, cos(cam.tilt), sin(cam.tilt));

  //form rotation and translation matrices
  for (unsigned int i=0; i<3; i++)
  {
    gsl_matrix_set(renderman.trans, i, 3, gsl_vector_get(cam.pose, i));
    gsl_matrix_set(renderman.rot,   i, 0, gsl_vector_get(renderman.camU, i));
    gsl_matrix_set(renderman.rot,   i, 1, gsl_vector_get(renderman.camV, i));
    gsl_matrix_set(renderman.rot,   i, 2, gsl_vector_get(renderman.camW, i));
  }

  //calculate total transformation matrix
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  1, renderman.rot, renderman.trans,
                 0, renderman.coordinateTransform);

  //calculate inverse of this transorm
  

  //now have start point in each loop, define each ray 
  //as starting there and having a unit vector pointing towards the point zoom distance away from
  //the centerx, centery point and normal to the imaging plane
  // (so vector from point to centerx,y + normal vector. Then normalized)

  //calculate normal point (tip of normal vector) of imaging plane with zoom length
  // (thats 0,0,cam.zoom)
  double normalx = 0;
  double normaly = 0;
  double normalz = cam.zoom;

  for (int row=0; row<numRaysRows; row++)
  {
     for (int col=0; col<numRaysCols; col++)
     {
        //find point on imaging plane in camera coord
        double filmy = row * raydistance - centery;
        double filmx = col * raydistance - centerx;
        double filmz = -cam.zoom;
        //take (normal point - film point) and normalize to get unit direction 
        // vector in camera coord

        //apply inverse tf to film point to get origin in world coord
        //apply inverse tf to unit direction vector to get udv in world coord
        
     }

  }

  gsl_vector_free (renderman.camU);
  gsl_vector_free (renderman.camV);
  gsl_vector_free (renderman.camW);
  gsl_matrix_free (renderman.rot); 
  gsl_matrix_free (renderman.trans);
  gsl_matrix_free (renderman.boxpoints);
  gsl_matrix_free (renderman.renderBoxpoints);
  
}

void render(camera cam, triangleSurface * surfs, unsigned int numSurfs/*, output image*/)
{

}