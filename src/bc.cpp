#include "Kokkos_Core.hpp"
#include "lsdebug.hpp"
#include "input.hpp"
#include "mpi.hpp"
#include <mpi.h>

struct bc_L {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_L(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int j, const int k, const int v) const {
        for (int i=0; i<ng; ++i)
            if (v==0 && bc_type==1) u(n-i-1,j,k,v) = -u(n+i,j,k,v);
            else u(n-i-1,j,k,v) = u(n+i,j,k,v);
    }
};

struct bc_R {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_R(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int j, const int k, const int v) const {
        for (int i=0; i<ng; ++i)
            if (v==0 && bc_type==1) u(n+i,j,k,v) = -u(n-i-1,j,k,v);
            else u(n+i,j,k,v) = u(n-i-1,j,k,v);
    }
};

struct bc_B {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_B(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int k, const int v) const {
        for (int j=0; j<ng; ++j)
            if (v==1 && bc_type==1) u(i,n-j-1,k,v) = -u(i,n+j,k,v);
            else u(i,n-j-1,k,v) = u(i,n+j,k,v);
    }
};

struct bc_T {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_T(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int k, const int v) const {
        for (int j=0; j<ng; ++j)
            if (v==1 && bc_type==1) u(i,n+j,k,v) = -u(i,n-j-1,k,v);
            else u(i,n+j,k,v) = u(i,n-j-1,k,v);
    }
};

struct bc_H {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_H(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {
        for (int k=0; k<ng; ++k)
            if (v==2 && bc_type==1) u(i,j,n-k-1,v) = -u(i,j,n+k,v);
            else u(i,j,n-k-1,v) = u(i,j,n+k,v);
    }
};

struct bc_F {
    Kokkos::View<double****> u;
    int n, ng, bc_type;

    bc_F(int n_, int ng_, int bc_type_, Kokkos::View<double****> u_) : n(n_), ng(ng_), bc_type(bc_type_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {
        for (int k=0; k<ng; ++k)
            if (v==2 && bc_type==1) u(i,j,n+k,v) = -u(i,j,n-k-1,v);
            else u(i,j,n+k,v) = u(i,j,n-k-1,v);
    }
};

void applyBCs(struct inputConfig cf, Kokkos::View<double****> &u){

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_bl;

    haloExchange(cf,u);

    if (cf.xMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngj,cf.ngk,cf.nv+5}), bc_L(cf.ng,cf.ng,cf.bcL,u));
    }

    if (cf.xPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngj,cf.ngk,cf.nv+5}), bc_R(cf.ng+cf.nci,cf.ng,cf.bcR,u));
    }

    if (cf.yMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngk,cf.nv+5}), bc_B(cf.ng,cf.ng,cf.bcB,u));
    }

    if (cf.yPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngk,cf.nv+5}), bc_T(cf.ng+cf.ncj,cf.ng,cf.bcT,u));
    }

    if (cf.zMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngj,cf.nv+5}), bc_H(cf.ng,cf.ng,cf.bcH,u));
    }

    if (cf.zPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ngi,cf.ngj,cf.nv+5}), bc_F(cf.ng+cf.nck,cf.ng,cf.bcF,u));
    }
}
