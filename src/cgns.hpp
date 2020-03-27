#include "fiesta.hpp"
#ifndef CGNS_H
#define CGNS_H

#include "input.hpp"
extern "C" {
#include "cgnslib.h"
#include "pcgnslib.h"
}

class cgnsWriter {

public:

    cgnsWriter(struct inputConfig, FS4D gridD, FS4D varD);
    struct inputConfig writeGrid(struct inputConfig cf, const FS4D gridD, const char * fname);
    struct inputConfig writeSPGrid(struct inputConfig cf, const FS4D gridD, const char * fname);

    void writeSolution(struct inputConfig cf, const FS4D gridD, const FS4D deviceV, int tdx, double time);
    void writeRestart(struct inputConfig cf, const FS4D gridD, const FS4D deviceV, int tdx, double time);

    void readSolution(struct inputConfig cf, FS4D &deviceG, FS4D &deviceV);

private:

    double *x;
    double *y;
    double *z;
    float *xsp;
    float *ysp;
    float *zsp;

    double *v;
    float *vsp;
    double * readV;

    FS4DH gridH;
    FS4DH varH;

};

#endif
