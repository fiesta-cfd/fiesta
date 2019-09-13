#include "Kokkos_Core.hpp"
#include "lsdebug.hpp"
#include "input.hpp"
#include "mpi_init.hpp"
#include <mpi.h>

struct bc_L {
    Kokkos::View<double****> u;
    int n, ng;

    bc_L(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int j, const int k, const int v) const {
        for (int i=0; i<ng; ++i)
            u(n-i-1,j,k,v) = u(n+i,j,k,v);
    }
};

struct bc_R {
    Kokkos::View<double****> u;
    int n, ng;

    bc_R(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int j, const int k, const int v) const {
        for (int i=0; i<ng; ++i)
            u(n+i+1,j,k,v) = u(n-i,j,k,v);
    }
};

struct bc_B {
    Kokkos::View<double****> u;
    int n, ng;

    bc_B(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int k, const int v) const {
        for (int j=0; j<ng; ++j)
            u(i,n-j-1,k,v) = u(i,n+j,k,v);
    }
};

struct bc_T {
    Kokkos::View<double****> u;
    int n, ng;

    bc_T(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int k, const int v) const {
        for (int j=0; j<ng; ++j)
            u(i,n+j+1,k,v) = u(i,n-j,k,v);
    }
};

struct bc_H {
    Kokkos::View<double****> u;
    int n, ng;

    bc_H(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {
        for (int k=0; k<ng; ++k)
            u(i,j,n-k-1,v) = u(i,j,n+k,v);
    }
};

struct bc_F {
    Kokkos::View<double****> u;
    int n, ng;

    bc_F(int n_, int ng_, Kokkos::View<double****> u_) : n(n_), ng(ng_), u(u_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int v) const {
        for (int k=0; k<ng; ++k)
            u(i,j,n+k+1,v) = u(i,j,n-k,v);
    }
};

void applyBCs(struct inputConfig cf, Kokkos::View<double****> &u){

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_bl;

    haloExchange(cf,u);

    if (cf.xMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ncj,cf.nck,cf.nv}), bc_L(cf.ng,cf.ng,u));
    }

    if (cf.xPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.ncj,cf.nck,cf.nv}), bc_R(cf.ng+cf.nci,cf.ng,u));
    }

    if (cf.yMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.nci,cf.nck,cf.nv}), bc_B(cf.ng,cf.ng,u));
    }

    if (cf.yPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.nci,cf.nck,cf.nv}), bc_T(cf.ng+cf.ncj,cf.ng,u));
    }

    if (cf.zMinus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.nci,cf.ncj,cf.nv}), bc_H(cf.ng,cf.ng,u));
    }

    if (cf.zPlus < 0){
        Kokkos::parallel_for(policy_bl({0,0,0},{cf.nci,cf.ncj,cf.nv}), bc_F(cf.ng+cf.nck,cf.ng,u));
    }
}
