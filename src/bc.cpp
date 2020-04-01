#include "Kokkos_Core.hpp"
#include "lsdebug.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include <mpi.h>
#include "mesh/mesh.h"

struct bc_L {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nlft;
    Kokkos::View<int*> nrht;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    //bc_L(int n_, int ng_, int bc_type_, Kokkos::View<double**> u_, int 26_, int 46_, Kokkos::View<int*> nlft_, Kokkos::View<int*> nrht_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_), 26(26_), 46(46_), nlft(nlft_), nrht(nrht_) {}
    bc_L(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nlft_, Kokkos::View<int*> nrht_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nlft(nlft_), nrht(nrht_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int j, const int k, const int v) const {
        int n = ng;
        int idx = n + j * 26 + k * 26 * 46;
        int idx_bc = nlft(idx);
        for (int i=0; i<ng; ++i) {
            //int idx_bc = n + i + j * 26 + k * 26 * 46; // index for the boundary cell
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = nlft(idx_bc);
            idx = nrht(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = nlft(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int i=0; i<ng; ++i) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = nlft(idx_bc);
                idx = nrht(idx);
            }
        }
    }
};

struct bc_R {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nlft;
    Kokkos::View<int*> nrht;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    //bc_R(int n_, int ng_, int bc_type_, Kokkos::View<double**> u_, int 26_, int 46_, Kokkos::View<int*> nlft_, Kokkos::View<int*> nrht_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_), 26(26_), 46(46_), nlft(nlft_), nrht(nrht_) {}
    bc_R(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nlft_, Kokkos::View<int*> nrht_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nlft(nlft_), nrht(nrht_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int j, const int k, const int v) const {
        int n = 23;
        int idx = n - 1 + j * 26 + k * 26 * 46;
        int idx_bc = nrht(idx);
        for (int i=0; i<ng; ++i) {
            //int idx_bc = n + i + j * 26 + k * 26 * 46; // index for the boundary cell
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = nrht(idx_bc);
            idx = nlft(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = nrht(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int i=0; i<ng; ++i) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = nrht(idx_bc);
                idx = nlft(idx);
            }
        }
    }
};

struct bc_B {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nbot;
    Kokkos::View<int*> ntop;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    //bc_B(int n_, int ng_, int bc_type_, Kokkos::View<double**> u_, int 26_, int 46_, Kokkos::View<int*> nbot_, Kokkos::View<int*> ntop_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_), 26(26_), 46(46_), nbot(nbot_), ntop(ntop_) {}
    bc_B(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nbot_, Kokkos::View<int*> ntop_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nbot(nbot_), ntop(ntop_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int i, const int k, const int v) const {
        int n = ng;
        int idx = i + n * 26 + k * 26 * 46;
        int idx_bc = nbot(idx); // index for the boundary cell
        for (int j=0; j<ng; ++j) {
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = nbot(idx_bc);
            idx = ntop(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = nbot(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int j=0; j<ng; ++j) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = nbot(idx_bc);
                idx = ntop(idx);
            }
        }
    }
};

struct bc_T {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nbot;
    Kokkos::View<int*> ntop;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    //bc_T(int n_, int ng_, int bc_type_, Kokkos::View<double**> u_, int 26_, int 46_, Kokkos::View<int*> nbot_, Kokkos::View<int*> ntop_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_), 26(26_), 46(46_), nbot(nbot_), ntop(ntop_) {}
    bc_T(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nbot_, Kokkos::View<int*> ntop_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nbot(nbot_), ntop(ntop_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int i, const int k, const int v) const {
        int n = 43;
        int idx = i + (n - 1) * 26 + k * 26 * 46;
        int idx_bc = ntop(idx); // index for the boundary cell
        for (int j=0; j<ng; ++j) {
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = ntop(idx_bc);
            idx = nbot(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = ntop(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int j=0; j<ng; ++j) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = ntop(idx_bc);
                idx = nbot(idx);
            }
        }
    }
};

struct bc_H {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nbak;
    Kokkos::View<int*> nfrt;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    //bc_H(int n_, int ng_, int bc_type_, Kokkos::View<double**> u_, int 26_, int 46_, Kokkos::View<int*> nbak_, Kokkos::View<int*> nfrt_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_), 26(26_), 46(46_), nbak(nbak_), nfrt(nfrt_) {}
    bc_H(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nbak_, Kokkos::View<int*> nfrt_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nbak(nbak_), nfrt(nfrt_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int i, const int j, const int v) const {
        int n = ng;
        int idx = i + j * 26 + n * 26 * 46;
        int idx_bc = nbak(idx); // index for the boundary cell
        for (int k=0; k<ng; ++k) {
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = nbak(idx_bc);
            idx = nfrt(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = nbak(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int k=0; k<ng; ++k) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = nbak(idx_bc);
                idx = nfrt(idx);
            }
        }
    }
};

struct bc_F {
    Kokkos::View<double**> u;
    Kokkos::View<int*> nbak;
    Kokkos::View<int*> nfrt;
    Kokkos::View<int*> cell_type;
    int ng, bc_type;

    bc_F(int ng_, int bc_type_, Kokkos::View<double**> u_, Kokkos::View<int*> nbak_, Kokkos::View<int*> nfrt_, Kokkos::View<int*> cell_type_) : ng(ng_), bc_type(bc_type_), u(u_), nbak(nbak_), nfrt(nfrt_), cell_type(cell_type_) {}

    KOKKOS_INLINE_FUNCTION
    /*
    void operator()(const int i, const int j, const int v) const {
        int n = ng;
        int idx = i + j * 26 + (n - 1) * 26 * 46;
        int idx_bc = nfrt(idx); // index for the boundary cell
        for (int k=0; k<ng; ++k) {
            if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
            else u(idx_bc,v) = u(idx,v);
            idx_bc = nfrt(idx_bc);
            idx = nbak(idx);
        }
    }
    */
    void operator()(const int n, const int v) const {
        int idx = n;
        int idx_bc = nfrt(idx); // index for the boundary cell
        //if (cell_type(idx) == BORDER_CELL_ && cell_type(idx_bc) == BOUNDARY_CELL_) {
        if (cell_type(idx) == REAL_CELL && cell_type(idx_bc) != REAL_CELL) {
            for (int k=0; k<ng; ++k) {
                if (v==0 && bc_type==1) u(idx_bc,v) = -u(idx,v);
                else u(idx_bc,v) = u(idx,v);
                idx_bc = nfrt(idx_bc);
                idx = nbak(idx);
            }
        }
    }
};

void applyBCs(struct inputConfig cf, Kokkos::View<double**> &u, class mpiBuffers &m, Kokkos::View<int*> nlft, Kokkos::View<int*> nrht, Kokkos::View<int*> nbot, Kokkos::View<int*> ntop, Kokkos::View<int*> nbak, Kokkos::View<int*> nfrt, Kokkos::View<int*> cell_type, int ncells){

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_bl;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy_b2;

    //haloExchange(cf,u,m);

    int cv = 0;
    if (cf.ceq == 1)
        cv = 5;

    if (cf.xMinus < 0){
        //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngj,cf.ngk,cf.nv+cv}), bc_L(cf.ng,cf.bcL,u,nlft,nrht,cell_type));
        Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_L(cf.ng,cf.bcL,u,nlft,nrht,cell_type));
    }

    if (cf.xPlus < 0){
        //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngj,cf.ngk,cf.nv+cv}), bc_R(cf.ng,cf.bcR,u,nlft,nrht,cell_type));
        Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_R(cf.ng,cf.bcR,u,nlft,nrht,cell_type));
    }

    if (cf.yMinus < 0){
        //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngk,cf.nv+cv}), bc_B(cf.ng,cf.bcB,u,nbot,ntop,cell_type));
        Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_B(cf.ng,cf.bcB,u,nbot,ntop,cell_type));
    }

    if (cf.yPlus < 0){
        //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngk,cf.nv+cv}), bc_T(cf.ng,cf.bcT,u,nbot,ntop,cell_type));
        Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_T(cf.ng,cf.bcT,u,nbot,ntop,cell_type));
    }

    if (cf.ndim == 3){
        if (cf.zMinus < 0){
            //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngj,cf.nv+cv}), bc_H(cf.ng,cf.bcH,u,nbak,nfrt,cell_type));
            Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_H(cf.ng,cf.bcH,u,nbak,nfrt,cell_type));
        }

        if (cf.zPlus < 0){
            //Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngj,cf.nv+cv}), bc_F(cf.ng,cf.bcF,u,nbak,nfrt,cell_type));
            Kokkos::parallel_for(policy_b2({0,0},{ncells,cf.nv+cv}), bc_F(cf.ng,cf.bcF,u,nbak,nfrt,cell_type));
        }
    }
}
