#ifndef WRITER_H
#define WRITER_H

#include "fiesta.hpp"
#include "input.hpp"
#include "Kokkos_Core.hpp"
#include <map>
#include <string>
#include "timer.hpp"

class writer
{

public:
    //writer(struct inputConfig &cf_, FS4D gridD, FS4D varD);
    writer();

    virtual void writeGrid(struct inputConfig cf, const FS4D gridD, const char * fname) = 0;
    virtual void writeSPGrid(struct inputConfig cf, const FS4D gridD, const char * fname) = 0;

    virtual void writeSolution(struct inputConfig cf, const FS4D gridD, const FS4D varD, int tdx, double time) = 0;
    virtual void writeRestart(struct inputConfig cf, const FS4D gridD, const FS4D varD, int tdx, double time) = 0;

    virtual void readSolution(struct inputConfig cf, FS4D &gridD, FS4D &varD) = 0;

//protected:
//    double *xdp;
//    double *ydp;
//    double *zdp;
//    float *xsp;
//    float *ysp;
//    float *zsp;
//
//    double *v;
//    float *vsp;
//    double *readV;
//
//    FS4DH gridH;
//    FS4DH varH;
};

#endif
