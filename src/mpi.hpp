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

#ifndef MPI_INIT_H
#define MPI_INIT_H

#include "Kokkos_Core.hpp"
#include "kokkosTypes.hpp"
#include "input.hpp"
#include "mpi.h"

void mpi_init(struct inputConfig &c);

class mpiHaloExchange {
  public:
    virtual void sendHalo(MPI_Request reqs[]) = 0;
    virtual void receiveHalo(MPI_Request reqs[]) = 0;
    virtual void unpackHalo() { };
    void packFace(std::vector<int> ion, FS4D &var, FS4D &buff);
    void unpackFace(std::vector<int> ion, FS4D &var, FS4D &buff);

    void pack(std::vector<int> ion, FS4D &var, FS4D &buff);
    void unpack(std::vector<int> ion, FS4D &var, FS4D &buff);

    FS4D &deviceV;
    struct inputConfig &cf;
    virtual void haloExchange();

    mpiHaloExchange(struct inputConfig &c, FS4D &v) : cf(c), deviceV(v) { }; 
};

class directHaloExchange : public mpiHaloExchange 
{
  public:
    MPI_Datatype leftRecvSubArray, rightRecvSubArray;
    MPI_Datatype bottomRecvSubArray, topRecvSubArray;
    MPI_Datatype backRecvSubArray, frontRecvSubArray;

    MPI_Datatype leftSendSubArray, rightSendSubArray;
    MPI_Datatype bottomSendSubArray, topSendSubArray;
    MPI_Datatype backSendSubArray, frontSendSubArray;

    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);

    directHaloExchange(inputConfig &c, FS4D &v);
};

class packedHaloExchange : public mpiHaloExchange 
{

  public:
    FS4D leftSend, leftRecv;
    FS4D rightSend, rightRecv;
    FS4D bottomSend, bottomRecv;
    FS4D topSend, topRecv;
    FS4D backSend, backRecv;
    FS4D frontSend, frontRecv;

    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);
    virtual void unpackHalo();

    packedHaloExchange(inputConfig &c, FS4D &v);
};

class copyHaloExchange : public packedHaloExchange 
{
  public:
    FS4DH leftSend_H, leftRecv_H;
    FS4DH rightSend_H, rightRecv_H;
    FS4DH bottomSend_H, bottomRecv_H;
    FS4DH topSend_H, topRecv_H;
    FS4DH backSend_H, backRecv_H;
    FS4DH frontSend_H, frontRecv_H;

    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);
    virtual void unpackHalo();

    copyHaloExchange(inputConfig &c, FS4D &v);
};

class orderedHaloExchange : public mpiHaloExchange 
{
  public:
    FS4D leftSend, leftRecv;
    FS4D rightSend, rightRecv;
    FS4D bottomSend, bottomRecv;
    FS4D topSend, topRecv;
    FS4D backSend, backRecv;
    FS4D frontSend, frontRecv;

    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);
    virtual void unpackHalo();
    virtual void haloExchange();
    //void pack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);
    //void unpack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);

    orderedHaloExchange(inputConfig &c, FS4D &v);
};

class orderedHostHaloExchange : public mpiHaloExchange 
{
  public:
    FS4D leftSend, leftRecv;
    FS4D rightSend, rightRecv;
    FS4D bottomSend, bottomRecv;
    FS4D topSend, topRecv;
    FS4D backSend, backRecv;
    FS4D frontSend, frontRecv;
    FS4DH leftSend_H, leftRecv_H;
    FS4DH rightSend_H, rightRecv_H;
    FS4DH bottomSend_H, bottomRecv_H;
    FS4DH topSend_H, topRecv_H;
    FS4DH backSend_H, backRecv_H;
    FS4DH frontSend_H, frontRecv_H;

    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);
    virtual void unpackHalo();
    virtual void haloExchange();
    //void pack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);
    //void unpack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);

    orderedHostHaloExchange(inputConfig &c, FS4D &v);
};

class unorderedHaloExchange : public mpiHaloExchange 
{
  public:
    virtual void sendHalo(MPI_Request reqs[]);
    virtual void receiveHalo(MPI_Request reqs[]);
    virtual void unpackHalo();
    virtual void haloExchange();
    //void pack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);
    //void unpack(const int ih, const int jh, const int kh, FS4D &var, FS4D &buff);

    unorderedHaloExchange(inputConfig &c, FS4D &v);
  private:
    FS4D sendBuffers[26];
    FS4D recvBuffers[26];
    int buffSize[26];
};
#endif

