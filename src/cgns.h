#ifndef CGNS_H
#define CGNS_H

#include "input.h"
#include "cgnslib.h"
#include "pcgnslib.h"

struct inputConfig writeGrid(struct inputConfig cf, double *x, double *y, double *z,char * fname);

void writeSolution(struct inputConfig cf, double *x, double *y, double *z, double * v, int tdx, double time);

#endif