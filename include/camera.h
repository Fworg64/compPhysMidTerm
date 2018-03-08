#ifndef CAMERA_H
#define CAMERA_H
//camera struct with properties
#include <gsl/gsl_blas.h>
//functions for moving camera around and such
typedef struct
{
  gsl_vector * pose;
  double pan;
  double tilt;
  double roll;
  double zoom;

} camera;

#endif