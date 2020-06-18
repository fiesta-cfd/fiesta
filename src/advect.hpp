#ifndef ADVECT_H
#define ADVECT_H

struct advect2D {
    FS4D dvar;
    FS2D fluxx,fluxy;
    FS1D cd;
    int v;

    advect2D (FS4D d_, FS2D fx_, FS2D fy_, FS1D cd_, int v_)
        : dvar(d_), fluxx(fx_), fluxy(fy_), cd(cd_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        dvar(i,j,0,v) = -( (fluxx(i,j) - fluxx(i-1,j))
                          +(fluxy(i,j) - fluxy(i,j-1)) );
    }
};
#endif
