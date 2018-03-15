#ifndef TRACING_RAY_H
#define TRACING_RAY_H

typedef struct
{
  double origin[3]; //point determining origin
  double direction[3]; //unit vector determining direction
  double hue, sat, lightness;
  unsigned int reflections; //number of reflections the ray has done
} tracingRay;

typedef struct
{
  tracingRay ray;
  double dist;

} intersectionResults;

typedef enum coll_E
{
  INTERSECT,
  MISS,
  ABSORBED,
} collisionState;

#endif