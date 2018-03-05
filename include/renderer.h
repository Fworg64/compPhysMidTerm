#ifndef RENDERER_H
#define RENDERER_H

#include "tracingRay.h"
#include "camera.h"

void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRaysRows, unsigned int numRaysCols);
void render(camera cam, triangularSurface * surfs, unsigned int numSurfs);

#endif