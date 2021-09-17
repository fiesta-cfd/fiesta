/*
  Copyright 2019-2021 The University of New Mexico

  This file is part of FIESTA.
  
  FIESTA is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.
  
  FIESTA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with FIESTA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "xdmf.hpp"

using namespace std;

void writeXMFDataItem(FILE* xmf, string path, int ndim, size_t *dims){
  if (ndim == 1)
    fprintf(xmf, "       <DataItem Dimensions=\"%zu\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0]);
  else if (ndim == 2)
    fprintf(xmf, "       <DataItem Dimensions=\"%zu %zu\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[1]);
  else
    fprintf(xmf, "       <DataItem Dimensions=\"%zu %zu %zu\" NumberType=\"Float\" " "Precision=\"4\" Format=\"HDF\">\n",dims[0],dims[1],dims[2]);
  fprintf(xmf, "        %s\n", path.c_str());
  fprintf(xmf, "       </DataItem>\n");
}
void writeXMF(string fname, string hname, int gridType, double time,
               int ndim, size_t *in_dims, vector<double> origin, vector<double> dx, int nvt, bool writeVarx,
               vector<string> vNames, vector<string> vxNames ){

    size_t gdims[ndim];
    size_t dims[ndim];
    for (int i=0; i<ndim; ++i){
      dims[i] = in_dims[ndim-1-i];
      gdims[i] = dims[i]+1;
    }
    //for (int i : dims) gdims.push_back(i+1);

    stringstream path;

    FILE *xmf = 0;
    xmf = fopen(fname.c_str(), "w");

    // Header
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    fprintf(xmf, " <Domain>\n");


    // Grid Header
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Time Value=\"%e\" />\n",time);
    if (gridType==0){
        if (ndim == 2){
          fprintf(xmf, "     <Topology TopologyType=\"2DCORECTMesh\" Dimensions=\"%zu %zu\"/>\n", dims[0]+1, dims[1]+1);
          fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n");
          fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"2\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
          fprintf(xmf, "            %f %f\n",origin[0],origin[1]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"2\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
          fprintf(xmf, "            %f %f\n",dx[0],dx[1]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "     </Geometry>\n");
        }else{
          fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%zu %zu %zu\"/>\n", dims[0]+1, dims[1]+1, dims[2]+1);
          fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
          fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
          fprintf(xmf, "            %f %f %f\n",origin[0],origin[1],origin[2]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
          fprintf(xmf, "            %f %f %f\n",dx[0],dx[1],dx[2]);
          fprintf(xmf, "       </DataItem>\n");
          fprintf(xmf, "     </Geometry>\n");
        }
    }else{
      if (ndim == 1){
        fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%zu\"/>\n", dims[0]+1);
        fprintf(xmf, "     <Geometry GeometryType=\"X\">\n");
      }else if (ndim == 2){
        fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%zu %zu\"/>\n", dims[0]+1, dims[1]+1);
        fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
      }else{
        fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%zu %zu %zu\"/>\n", dims[0]+1, dims[1]+1, dims[2]+1);
        fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
      }

      // Grid Coordinate Arrays
      for (int d=0; d<ndim; ++d){
        //for (int d : actDims){
        path.str("");
        path << hname << ":/Grid/Dimension" << d;
        writeXMFDataItem(xmf, path.str(), ndim, gdims);
      }
      fprintf(xmf, "     </Geometry>\n");
    }

    // Momentum Vector
    fprintf(xmf, "     <Attribute Name=\"Momentum\" AttributeType=\"Vector\" " "Center=\"Cell\">\n");
    if (ndim == 1)
      fprintf(xmf, "      <DataItem Dimensions=\"%zu 1\" Function=\"JOIN($0)\" " "ItemType=\"Function\">\n",dims[0]);
    else if (ndim == 2)
      fprintf(xmf, "      <DataItem Dimensions=\"%zu %zu 2\" Function=\"JOIN($0,$1)\" " "ItemType=\"Function\">\n",dims[0],dims[1]);
    else
      fprintf(xmf, "      <DataItem Dimensions=\"%zu %zu %zu 3\" Function=\"JOIN($0,$1,$2)\" " "ItemType=\"Function\">\n", dims[0], dims[1],dims[2]);

    for (int d=0; d<ndim; ++d){
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << d;
      writeXMFDataItem(xmf, path.str(), ndim, dims);
    }
    fprintf(xmf, "      </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    // Other Variables
    fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n");
    path.str("");
    path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << ndim;
    writeXMFDataItem(xmf, path.str(), ndim, dims);
    fprintf(xmf, "     </Attribute>\n");

    for (int var = ndim+1; var < nvt; ++var) {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n",vNames[var].c_str());
      path.str("");
      path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var;
      writeXMFDataItem(xmf, path.str(), ndim, dims);
      fprintf(xmf, "     </Attribute>\n");
    }

    if (writeVarx){
      // Velocity Vector
      fprintf(xmf, "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" " "Center=\"Cell\">\n");
      if (ndim == 1)
        fprintf(xmf, "      <DataItem Dimensions=\"%zu 1\" Function=\"JOIN($0)\" " "ItemType=\"Function\">\n",dims[0]);
      else if (ndim == 2)
        fprintf(xmf, "      <DataItem Dimensions=\"%zu %zu 2\" Function=\"JOIN($0,$1)\" " "ItemType=\"Function\">\n",dims[0],dims[1]);
      else
        fprintf(xmf, "      <DataItem Dimensions=\"%zu %zu %zu 3\" Function=\"JOIN($0,$1,$2)\" " "ItemType=\"Function\">\n", dims[0], dims[1],dims[2]);

      for (int d=0; d<ndim; ++d){
        path.str("");
        path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << nvt+d;
        writeXMFDataItem(xmf, path.str(), ndim, dims);
      }
      fprintf(xmf, "      </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");

      // Other Extra Variables
      for (size_t var = ndim; var < vxNames.size(); ++var) {
        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" " "Center=\"Cell\">\n",vxNames[var].c_str());
        path.str("");
        path << hname << ":/Solution/Variable" << setw(2) << setfill('0') << var+nvt;
        writeXMFDataItem(xmf, path.str(), ndim, dims);
        fprintf(xmf, "     </Attribute>\n");
      }
    }

    // End Field Variables
    fprintf(xmf, "   </Grid>\n");

    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}
