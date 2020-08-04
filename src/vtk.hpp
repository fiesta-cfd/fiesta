#include "fiesta.hpp"
#ifndef VTK_H
#define VTK_H

#include "input.hpp"
#include "rkfunction.hpp"
#include "writer.hpp"

class serialVTKWriter : public writer {

public:
  serialVTKWriter(struct inputConfig, FS4D gridD, FS4D varD);
  void writeGrid(struct inputConfig cf, const FS4D gridD, const char *fname);
  void writeSPGrid(struct inputConfig cf, const FS4D gridD, const char *fname);

  void writeSolution(struct inputConfig cf, rk_func *f, int tdx, double time);
  void writeRestart(struct inputConfig cf, rk_func *f, int tdx, double time);

  void readSolution(struct inputConfig cf, FS4D &deviceG, FS4D &deviceV);

private:
  FS4DH gridH;
  FS4DH varH;
  //    FSP2DH particlesH;
};

#endif
