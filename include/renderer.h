#ifndef RENDERER_H
#define RENDERER_H

#include "tracingRay.h"
#include "camera.h"
#include "triangleSurface.h"
#include <gsl/gsl_blas.h>

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
};


void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols);
void render(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols, const char* filename);

#endif