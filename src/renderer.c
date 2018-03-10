#include "renderer.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>


void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols)
{
   printf("starting routine\n");
   struct renderman_S renderman;
   renderman.camU = gsl_vector_alloc(3);
   renderman.camV = gsl_vector_alloc(3);
   renderman.camW = gsl_vector_alloc(3);
   renderman.rot = gsl_matrix_alloc(4,4);
   renderman.trans = gsl_matrix_alloc(4,4);
   renderman.boxpoints = gsl_matrix_alloc(4,8);
   renderman.renderBoxpoints = gsl_matrix_alloc(4,8);
   renderman.coordinateTransform = gsl_matrix_alloc(4,4);
   gsl_matrix * inverseTransform = gsl_matrix_alloc(4,4);
   gsl_matrix * originPoint = gsl_matrix_alloc(4,1);
   gsl_matrix * unitDirectionVector = gsl_matrix_alloc(4,1);
   gsl_matrix * result = gsl_matrix_alloc(4,1);
   
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
  
  printf("forming transform matrices\n");
  //form rotation and translation matrices
  for (unsigned int i=0; i<3; i++)
  {
    gsl_matrix_set(renderman.trans, i, 3, gsl_vector_get(cam.pose, i));
    gsl_matrix_set(renderman.rot,   i, 0, gsl_vector_get(renderman.camU, i));
    gsl_matrix_set(renderman.rot,   i, 1, gsl_vector_get(renderman.camV, i));
    gsl_matrix_set(renderman.rot,   i, 2, gsl_vector_get(renderman.camW, i));
  }
  gsl_matrix_set(renderman.rot, 3,3,1);
  for (unsigned int i=0; i<4; i++)
       {gsl_matrix_set(renderman.trans, i,i,1);}

  printf("\nrotation matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.rot, i,0),
                              gsl_matrix_get(renderman.rot, i,1),
                              gsl_matrix_get(renderman.rot, i,2), 
                              gsl_matrix_get(renderman.rot, i,3));

  printf("\ntranslation matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.trans, i,0),
                              gsl_matrix_get(renderman.trans, i,1),
                              gsl_matrix_get(renderman.trans, i,2), 
                              gsl_matrix_get(renderman.trans, i,3));

  //calculate total transformation matrix
  printf("calculating total tf\n");
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,  1, renderman.trans, renderman.rot,
                 0, renderman.coordinateTransform);

  printf("\ntransform matrix: \n");
  for (unsigned int i=0; i<4; i++)
  printf("%f, %f, %f, %f,\n", gsl_matrix_get(renderman.coordinateTransform, i,0), 
                              gsl_matrix_get(renderman.coordinateTransform, i,1),
                              gsl_matrix_get(renderman.coordinateTransform, i,2), 
                              gsl_matrix_get(renderman.coordinateTransform, i,3));

  //calculate inverse of this transorm
  printf("calculating inverse\n");
  gsl_permutation * myPerm =  gsl_permutation_alloc (4);
  int signum;
  printf("calculating inverse2\n");
  gsl_linalg_LU_decomp(renderman.coordinateTransform, myPerm, &signum);
  printf("calculating inverse3\n");
  gsl_linalg_LU_invert(renderman.coordinateTransform, myPerm, inverseTransform);
  gsl_permutation_free(myPerm);
  //now have start point in each loop, define each ray 
  //as starting there and having a unit vector pointing towards the point zoom distance away from
  //the centerx, centery point and normal to the imaging plane
  // (so vector from point to centerx,y + normal vector. Then normalized)
  printf("going to generate rays\n");
  //calculate normal point (tip of normal vector) of imaging plane with zoom length
  // (thats 0,0,0)
  double normalx = 0;
  double normaly = 0;
  double normalz = 0;

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
        double magnitude = sqrt(filmy * filmy + filmx * filmx + normalz*normalz);
        
        double udvU = filmy/magnitude;
        double udvV = filmx/magnitude;
        double udvW = cam.zoom/magnitude;
        //apply inverse tf to film point to get origin in world coord
        gsl_matrix_set(originPoint, 0,0, filmx);
        gsl_matrix_set(originPoint,1,0,filmy);
        gsl_matrix_set(originPoint,2,0,filmz);
        gsl_matrix_set(originPoint,3,0,1);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverseTransform, 
                       originPoint, 0, result);
        
        rays[row + numRaysRows * col].origin[0] = gsl_matrix_get(result, 0,0);
        rays[row + numRaysRows * col].origin[1] = gsl_matrix_get(result, 1,0);
        rays[row + numRaysRows * col].origin[2] = gsl_matrix_get(result, 2,0);
        //apply inverse tf to unit direction vector to get udv in world coord
        gsl_matrix_set(unitDirectionVector, 0,0, udvU);
        gsl_matrix_set(unitDirectionVector,1,0,udvV);
        gsl_matrix_set(unitDirectionVector,2,0,udvW);
        gsl_matrix_set(unitDirectionVector,3,0,1);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverseTransform, 
                       unitDirectionVector, 0, result);

        rays[row + numRaysRows * col].direction[0] = gsl_matrix_get(result, 0,0);
        rays[row + numRaysRows * col].direction[1] = gsl_matrix_get(result, 1,0);
        rays[row + numRaysRows * col].direction[2] = gsl_matrix_get(result, 2,0);
     }

  }
  
  printf("freeing allocated memory");
  gsl_vector_free (renderman.camU);
  gsl_vector_free (renderman.camV);
  gsl_vector_free (renderman.camW);
  gsl_matrix_free (renderman.rot); 
  gsl_matrix_free (renderman.trans);
  gsl_matrix_free (renderman.boxpoints);
  gsl_matrix_free (renderman.renderBoxpoints);
  gsl_matrix_free (inverseTransform);
  gsl_matrix_free (originPoint);
  gsl_matrix_free (unitDirectionVector);
  gsl_matrix_free (result);
  
}

void render(camera cam, triangleSurface * surfs, unsigned int numSurfs/*, output image*/)
{

}