#ifndef CAMERA_H
#define CAMERA_H
//camera struct with properties

//functions for moving camera around and such
typedef struct
{
  double pose[3]; //x,y,z in world coord
  double pan;
  double tilt;
  double roll;
  double zoom;

} camera;

#endif