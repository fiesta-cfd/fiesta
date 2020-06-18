#ifndef presgrad_H
#define presgrad_H
struct applyPressureGradient2D {
    
    FS4D dvar;
    FS2D p;
    FS1D cd;

    applyPressureGradient2D (FS4D dvar_, FS2D p_, FS1D cd_)
        : dvar(dvar_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        double dx = cd(1);
        double dy = cd(2);
        // calculate pressure gradient across cell in each direction using 4th order
        // central difference
        double dxp = ( p(i-2,j) - 8.0*p(i-1,j) + 8.0*p(i+1,j) - p(i+2,j) )/(12.0*dx);
        double dyp = ( p(i,j-2) - 8.0*p(i,j-1) + 8.0*p(i,j+1) - p(i,j+2) )/(12.0*dy);

        //apply pressure gradient term to right hand side of Euler equation dV/dt = ...
        double p1,p2;
        p1 = dvar(i,j,0,0) - dxp;
        p2 = dvar(i,j,0,1) - dyp;

        dvar(i,j,0,0) = p1;
        dvar(i,j,0,1) = p2;
    }
};
struct applyGenPressureGradient2D {
    
    FS4D dvar;
    FS2D p;
    FS4D m;
    FS1D cd;

    applyGenPressureGradient2D (FS4D dvar_, FS4D m_, FS2D p_, FS1D cd_)
        : dvar(dvar_), m(m_), p(p_), cd(cd_) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int j) const {
        
        double dxipxix;
        double detpetx;
        double dxipxiy;
        double detpety;

        dxipxix = ( p(i-2,j)*m(i-2,j,0,0) - 8.0*p(i-1,j)*m(i-1,j,0,0) + 8.0*p(i+1,j)*m(i+1,j,0,0) - p(i+2,j)*m(i+2,j,0,0) )/(12.0);
        detpetx = ( p(i,j-2)*m(i,j-2,1,0) - 8.0*p(i,j-1)*m(i,j-1,1,0) + 8.0*p(i,j+1)*m(i,j+1,1,0) - p(i,j+2)*m(i,j+2,1,0) )/(12.0);
        dxipxiy = ( p(i-2,j)*m(i-2,j,0,1) - 8.0*p(i-1,j)*m(i-1,j,0,1) + 8.0*p(i+1,j)*m(i+1,j,0,1) - p(i+2,j)*m(i+2,j,0,1) )/(12.0);
        detpety = ( p(i,j-2)*m(i,j-2,1,1) - 8.0*p(i,j-1)*m(i,j-1,1,1) + 8.0*p(i,j+1)*m(i,j+1,1,1) - p(i,j+2)*m(i,j+2,1,1) )/(12.0);

        dvar(i,j,0,0) = ( dvar(i,j,0,0) - (dxipxix + detpetx) );
        dvar(i,j,0,1) = ( dvar(i,j,0,1) - (dxipxiy + detpety) );
    }
};
#endif
