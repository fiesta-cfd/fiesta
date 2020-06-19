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

struct advect3D {
    
    FS4D dvar;
    FS3D wenox;
    FS3D wenoy;
    FS3D wenoz;
    int v;

    advect3D (FS4D dvar_, FS3D wenox_, FS3D wenoy_, FS3D wenoz_, int v_)
        : dvar(dvar_), wenox(wenox_), wenoy(wenoy_), wenoz(wenoz_), v(v_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j, const int k) const {

        dvar(i,j,k,v) = -( (wenox(i,j,k) - wenox(i-1,j,k))
                          +(wenoy(i,j,k) - wenoy(i,j-1,k))
                          +(wenoz(i,j,k) - wenoz(i,j,k-1)) );
    }
};

#endif
