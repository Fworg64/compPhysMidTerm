#ifndef RENDERER_H
#define RENDERER_H

#include "tracingRay.h"
#include "camera.h"

void generateRaysFromCamera(camera cam, tracingRay * rays, unsigned int numRays);
void render(camera cam, triangularSurface * surfs, unsigned int numSurfs);

#endif