#include "fiesta.hpp"
#ifndef CGNS_H
#define CGNS_H

#include "input.hpp"
extern "C" {
#include "cgnslib.h"
#include "pcgnslib.h"
}

struct inputConfig writeGrid(struct inputConfig cf, double *x, double *y, double *z,char * fname);
struct inputConfig writeSPGrid(struct inputConfig cf, float *x, float *y, float *z,char * fname);

void writeSolution(struct inputConfig cf, float *x, float *y, float *z, const FS4D deviceV, int tdx, double time);
void writeRestart(struct inputConfig cf, double *x, double *y, double *z, const FS4D deviceV, int tdx, double time);

void readSolution(struct inputConfig cf, const FS4D deviceV);

#endif
