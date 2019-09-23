#include "input.hpp"
#include "mpi_init.hpp"
#include "cgns.hpp"
#include "bc.hpp"
#include "Kokkos_Core.hpp"
#include <mpi.h>
#include "lsdebug.hpp"

//template < class ViewType >
struct rhs_func {

    Kokkos::View<double****> mvar;
    Kokkos::View<double****> mdvar;
    Kokkos::View<double*> mcd;

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_f;
    struct inputConfig cf;
    double dt;

    rhs_func(double dt_, struct inputConfig &cf_, const Kokkos::View<double****> & u_, Kokkos::View<double****> & k_, Kokkos::View<double*> & cd_)
        : dt(dt_), cf(cf_) , mvar(u_), mdvar(k_) , mcd(cd_){}

    void operator()() const {
        Kokkos::View<double****> dvar = mdvar;
        Kokkos::View<double****> var = mvar;
        Kokkos::View<double*> cd = mcd;

        Kokkos::parallel_for("Loopf", policy_f({cf.ng,cf.ng,cf.ng},{cf.ngi-cf.ng, cf.ngj-cf.ng, cf.ngk-cf.ng}),
        KOKKOS_LAMBDA __device__ (const int i, const int j, const int k) {
            
            int ns = (int)cd(0);

            double eps = 0.000001;

            double Cp = 0;
            double Cv = 0;
            double gammas,Rs;

            double ur,ul,vr,vl,wr,wl;
            double dxp,dyp,dzp;

            double f1,f2,f3,f4,f5;
            double b1,b2,b3, w1,w2,w3, p1,p2,p3;
            double wxl,wxr, wyl,wyr, wzl,wzr;

            double rho [19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double uvel[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double vvel[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double wvel[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double pres[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double gamma[19] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double mf[2][19];

            for (int s=0; s<ns; ++s){
                rho[0]  = rho[0]  + var(i  ,j  ,k  ,4+s);
                rho[1]  = rho[1]  + var(i-1,j  ,k  ,4+s);
                rho[2]  = rho[2]  + var(i+1,j  ,k  ,4+s);
                rho[3]  = rho[3]  + var(i  ,j-1,k  ,4+s);
                rho[4]  = rho[4]  + var(i  ,j+1,k  ,4+s);
                rho[5]  = rho[5]  + var(i  ,j  ,k-1,4+s);
                rho[6]  = rho[6]  + var(i  ,j  ,k+1,4+s);
                rho[7]  = rho[7]  + var(i-2,j  ,k  ,4+s);
                rho[8]  = rho[8]  + var(i+2,j  ,k  ,4+s);
                rho[9]  = rho[9]  + var(i  ,j-2,k  ,4+s);
                rho[10] = rho[10] + var(i  ,j+2,k  ,4+s);
                rho[11] = rho[11] + var(i  ,j  ,k-2,4+s);
                rho[12] = rho[12] + var(i  ,j  ,k+2,4+s);
                rho[13] = rho[13] + var(i-3,j  ,k  ,4+s);
                rho[14] = rho[14] + var(i+3,j  ,k  ,4+s);
                rho[15] = rho[15] + var(i  ,j-3,k  ,4+s);
                rho[16] = rho[16] + var(i  ,j+3,k  ,4+s);
                rho[17] = rho[17] + var(i  ,j  ,k-3,4+s);
                rho[18] = rho[18] + var(i  ,j  ,k+3,4+s);
            }

            for (int s=0; s<ns; ++s){
                mf[s][0]  = var(i  ,j  ,k  ,4+s)/rho[0] ;
                mf[s][1]  = var(i-1,j  ,k  ,4+s)/rho[1] ;
                mf[s][2]  = var(i+1,j  ,k  ,4+s)/rho[2] ;
                mf[s][3]  = var(i  ,j-1,k  ,4+s)/rho[3] ;
                mf[s][4]  = var(i  ,j+1,k  ,4+s)/rho[4] ;
                mf[s][5]  = var(i  ,j  ,k-1,4+s)/rho[5] ;
                mf[s][6]  = var(i  ,j  ,k+1,4+s)/rho[6] ;
                mf[s][7]  = var(i-2,j  ,k  ,4+s)/rho[7] ;
                mf[s][8]  = var(i+2,j  ,k  ,4+s)/rho[8] ;
                mf[s][9]  = var(i  ,j-2,k  ,4+s)/rho[9] ;
                mf[s][10] = var(i  ,j+2,k  ,4+s)/rho[10];
                mf[s][11] = var(i  ,j  ,k-2,4+s)/rho[11];
                mf[s][12] = var(i  ,j  ,k+2,4+s)/rho[12];
                mf[s][13] = var(i-3,j  ,k  ,4+s)/rho[13];
                mf[s][14] = var(i+3,j  ,k  ,4+s)/rho[14];
                mf[s][15] = var(i  ,j-3,k  ,4+s)/rho[15];
                mf[s][16] = var(i  ,j+3,k  ,4+s)/rho[16];
                mf[s][17] = var(i  ,j  ,k-3,4+s)/rho[17];
                mf[s][18] = var(i  ,j  ,k+3,4+s)/rho[18];
            }

            uvel[0]  = var(i  ,j  ,k  ,0)/rho[0];
            uvel[1]  = var(i-1,j  ,k  ,0)/rho[1];
            uvel[2]  = var(i+1,j  ,k  ,0)/rho[2];
            uvel[3]  = var(i  ,j-1,k  ,0)/rho[3];
            uvel[4]  = var(i  ,j+1,k  ,0)/rho[4];
            uvel[5]  = var(i  ,j  ,k-1,0)/rho[5];
            uvel[6]  = var(i  ,j  ,k+1,0)/rho[6];
            uvel[7]  = var(i-2,j  ,k  ,0)/rho[7];
            uvel[8]  = var(i+2,j  ,k  ,0)/rho[8];
            uvel[9]  = var(i  ,j-2,k  ,0)/rho[9];
            uvel[10] = var(i  ,j+2,k  ,0)/rho[10];
            uvel[11] = var(i  ,j  ,k-2,0)/rho[11];
            uvel[12] = var(i  ,j  ,k+2,0)/rho[12];
            uvel[13] = var(i-3,j  ,k  ,0)/rho[13];
            uvel[14] = var(i+3,j  ,k  ,0)/rho[14];
            uvel[15] = var(i  ,j-3,k  ,0)/rho[15];
            uvel[16] = var(i  ,j+3,k  ,0)/rho[16];
            uvel[17] = var(i  ,j  ,k-3,0)/rho[17];
            uvel[18] = var(i  ,j  ,k+3,0)/rho[18];

            vvel[0]  = var(i  ,j  ,k  ,1)/rho[0];
            vvel[1]  = var(i-1,j  ,k  ,1)/rho[1];
            vvel[2]  = var(i+1,j  ,k  ,1)/rho[2];
            vvel[3]  = var(i  ,j-1,k  ,1)/rho[3];
            vvel[4]  = var(i  ,j+1,k  ,1)/rho[4];
            vvel[5]  = var(i  ,j  ,k-1,1)/rho[5];
            vvel[6]  = var(i  ,j  ,k+1,1)/rho[6];
            vvel[7]  = var(i-2,j  ,k  ,1)/rho[7];
            vvel[8]  = var(i+2,j  ,k  ,1)/rho[8];
            vvel[9]  = var(i  ,j-2,k  ,1)/rho[9];
            vvel[10] = var(i  ,j+2,k  ,1)/rho[10];
            vvel[11] = var(i  ,j  ,k-2,1)/rho[11];
            vvel[12] = var(i  ,j  ,k+2,1)/rho[12];
            vvel[13] = var(i-3,j  ,k  ,1)/rho[13];
            vvel[14] = var(i+3,j  ,k  ,1)/rho[14];
            vvel[15] = var(i  ,j-3,k  ,1)/rho[15];
            vvel[16] = var(i  ,j+3,k  ,1)/rho[16];
            vvel[17] = var(i  ,j  ,k-3,1)/rho[17];
            vvel[18] = var(i  ,j  ,k+3,1)/rho[18];

            wvel[0]  = var(i  ,j  ,k  ,2)/rho[0];
            wvel[1]  = var(i-1,j  ,k  ,2)/rho[1];
            wvel[2]  = var(i+1,j  ,k  ,2)/rho[2];
            wvel[3]  = var(i  ,j-1,k  ,2)/rho[3];
            wvel[4]  = var(i  ,j+1,k  ,2)/rho[4];
            wvel[5]  = var(i  ,j  ,k-1,2)/rho[5];
            wvel[6]  = var(i  ,j  ,k+1,2)/rho[6];
            wvel[7]  = var(i-2,j  ,k  ,2)/rho[7];
            wvel[8]  = var(i+2,j  ,k  ,2)/rho[8];
            wvel[9]  = var(i  ,j-2,k  ,2)/rho[9];
            wvel[10] = var(i  ,j+2,k  ,2)/rho[10];
            wvel[11] = var(i  ,j  ,k-2,2)/rho[11];
            wvel[12] = var(i  ,j  ,k+2,2)/rho[12];
            wvel[13] = var(i-3,j  ,k  ,2)/rho[13];
            wvel[14] = var(i+3,j  ,k  ,2)/rho[14];
            wvel[15] = var(i  ,j-3,k  ,2)/rho[15];
            wvel[16] = var(i  ,j+3,k  ,2)/rho[16];
            wvel[17] = var(i  ,j  ,k-3,2)/rho[17];
            wvel[18] = var(i  ,j  ,k+3,2)/rho[18];
            
            // need to weight by mass fractions.  this means another array for gamma at each stencil point
            for(int idx=0; idx<19; ++idx){
                Cp = 0;
                Cv = 0;
                for (int s=0; s<ns; ++s){
                    gammas = cd(4+2*s);
                    Rs = cd(4+2*s+1);
                    
                    Cp = Cp + mf[s][idx]*( (gammas*Rs)/(gammas-1) );
                    Cv = Cv + mf[s][idx]*( Rs/(gammas-1) );
                }

                gamma[idx] = Cp/Cv;
            }

            pres[0]  = (gamma[0] -1)*(var(i  ,j  ,k  ,3)-0.5*rho[0]  *(uvel[0] *uvel[0]  + vvel[0] *vvel[0]  + wvel[0] *wvel[0] ));
            pres[1]  = (gamma[1] -1)*(var(i-1,j  ,k  ,3)-0.5*rho[1]  *(uvel[1] *uvel[1]  + vvel[1] *vvel[1]  + wvel[1] *wvel[1] ));
            pres[2]  = (gamma[2] -1)*(var(i+1,j  ,k  ,3)-0.5*rho[2]  *(uvel[2] *uvel[2]  + vvel[2] *vvel[2]  + wvel[2] *wvel[2] ));
            pres[3]  = (gamma[3] -1)*(var(i  ,j-1,k  ,3)-0.5*rho[3]  *(uvel[3] *uvel[3]  + vvel[3] *vvel[3]  + wvel[3] *wvel[3] ));
            pres[4]  = (gamma[4] -1)*(var(i  ,j+1,k  ,3)-0.5*rho[4]  *(uvel[4] *uvel[4]  + vvel[4] *vvel[4]  + wvel[4] *wvel[4] ));
            pres[5]  = (gamma[5] -1)*(var(i  ,j  ,k-1,3)-0.5*rho[5]  *(uvel[5] *uvel[5]  + vvel[5] *vvel[5]  + wvel[5] *wvel[5] ));
            pres[6]  = (gamma[6] -1)*(var(i  ,j  ,k+1,3)-0.5*rho[6]  *(uvel[6] *uvel[6]  + vvel[6] *vvel[6]  + wvel[6] *wvel[6] ));
            pres[7]  = (gamma[7] -1)*(var(i-2,j  ,k  ,3)-0.5*rho[7]  *(uvel[7] *uvel[7]  + vvel[7] *vvel[7]  + wvel[7] *wvel[7] ));
            pres[8]  = (gamma[8] -1)*(var(i+2,j  ,k  ,3)-0.5*rho[8]  *(uvel[8] *uvel[8]  + vvel[8] *vvel[8]  + wvel[8] *wvel[8] ));
            pres[9]  = (gamma[9] -1)*(var(i  ,j-2,k  ,3)-0.5*rho[9]  *(uvel[9] *uvel[9]  + vvel[9] *vvel[9]  + wvel[9] *wvel[9] ));
            pres[10] = (gamma[10]-1)*(var(i  ,j+2,k  ,3)-0.5*rho[10] *(uvel[10]*uvel[10] + vvel[10]*vvel[10] + wvel[10]*wvel[10]));
            pres[11] = (gamma[11]-1)*(var(i  ,j  ,k-2,3)-0.5*rho[11] *(uvel[11]*uvel[11] + vvel[11]*vvel[11] + wvel[11]*wvel[11]));
            pres[12] = (gamma[12]-1)*(var(i  ,j  ,k+2,3)-0.5*rho[12] *(uvel[12]*uvel[12] + vvel[12]*vvel[12] + wvel[12]*wvel[12]));
            pres[13] = (gamma[13]-1)*(var(i-3,j  ,k  ,3)-0.5*rho[13] *(uvel[13]*uvel[13] + vvel[13]*vvel[13] + wvel[13]*wvel[13]));
            pres[14] = (gamma[14]-1)*(var(i+3,j  ,k  ,3)-0.5*rho[14] *(uvel[14]*uvel[14] + vvel[14]*vvel[14] + wvel[14]*wvel[14]));
            pres[15] = (gamma[15]-1)*(var(i  ,j-3,k  ,3)-0.5*rho[15] *(uvel[15]*uvel[15] + vvel[15]*vvel[15] + wvel[15]*wvel[15]));
            pres[16] = (gamma[16]-1)*(var(i  ,j+3,k  ,3)-0.5*rho[16] *(uvel[16]*uvel[16] + vvel[16]*vvel[16] + wvel[16]*wvel[16]));
            pres[17] = (gamma[17]-1)*(var(i  ,j  ,k-3,3)-0.5*rho[17] *(uvel[17]*uvel[17] + vvel[17]*vvel[17] + wvel[17]*wvel[17]));
            pres[18] = (gamma[18]-1)*(var(i  ,j  ,k+3,3)-0.5*rho[18] *(uvel[18]*uvel[18] + vvel[18]*vvel[18] + wvel[18]*wvel[18]));

            ul = (-uvel[2]  + 7.0*uvel[0] + 7.0*uvel[1] - uvel[7])/12.0;
            ur = (-uvel[8]  + 7.0*uvel[2] + 7.0*uvel[0] - uvel[1])/12.0;
            
            vl = (-vvel[4]  + 7.0*vvel[0] + 7.0*vvel[3] - vvel[9])/12.0;
            vr = (-vvel[10] + 7.0*vvel[4] + 7.0*vvel[0] - vvel[3])/12.0;
            
            wl = (-wvel[6]  + 7.0*wvel[0] + 7.0*wvel[5] - wvel[11])/12.0;
            wr = (-wvel[12] + 7.0*wvel[6] + 7.0*wvel[0] - wvel[5] )/12.0;

            dxp = (pres[7]  - 8.0*pres[1] + 8.0*pres[2] - pres[8] )/(12.0*cd(1));
            dyp = (pres[9]  - 8.0*pres[3] + 8.0*pres[4] - pres[10])/(12.0*cd(2));
            dzp = (pres[11] - 8.0*pres[5] + 8.0*pres[6] - pres[12])/(12.0*cd(3));

            //printf("%f, %f, %f\n",dxp,dyp,dzp);

            //can maybe reduce local variable usage by converting momentumns to velities in expressions at the last moment
            //and by using the same array space for pressure and density, we only need one at a time, test this first

            for (int vv=0; vv <4+ns; ++vv){
                if (vv == 3){
                    if (ul < 0.0){
                        f1 = var(i+2,j,k,vv) + pres[8];
                        f2 = var(i+1,j,k,vv) + pres[2];
                        f3 = var(i  ,j,k,vv) + pres[0];
                        f4 = var(i-1,j,k,vv) + pres[1];
                        f5 = var(i-2,j,k,vv) + pres[7];
                    }else{              
                        f1 = var(i-3,j,k,vv) + pres[13];
                        f2 = var(i-2,j,k,vv) + pres[7];
                        f3 = var(i-1,j,k,vv) + pres[1];
                        f4 = var(i  ,j,k,vv) + pres[0];
                        f5 = var(i+1,j,k,vv) + pres[2];
                    }
                }else{
                    if (ul < 0.0){
                        f1 = var(i+2,j,k,vv);
                        f2 = var(i+1,j,k,vv);
                        f3 = var(i  ,j,k,vv);
                        f4 = var(i-1,j,k,vv);
                        f5 = var(i-2,j,k,vv);
                    }else{              
                        f1 = var(i-3,j,k,vv);
                        f2 = var(i-2,j,k,vv);
                        f3 = var(i-1,j,k,vv);
                        f4 = var(i  ,j,k,vv);
                        f5 = var(i+1,j,k,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wxl = ul*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                if (vv == 3){
                    if (ur < 0.0){
                        f1 = var(i+3,j,k,vv) + pres[14];
                        f2 = var(i+2,j,k,vv) + pres[8];
                        f3 = var(i+1,j,k,vv) + pres[2];
                        f4 = var(i  ,j,k,vv) + pres[0];
                        f5 = var(i-1,j,k,vv) + pres[1];
                    }else{              
                        f1 = var(i-2,j,k,vv) + pres[7];
                        f2 = var(i-1,j,k,vv) + pres[1];
                        f3 = var(i  ,j,k,vv) + pres[0];
                        f4 = var(i+1,j,k,vv) + pres[2];
                        f5 = var(i+2,j,k,vv) + pres[8];
                    }
                }else{
                    if (ur < 0.0){
                        f1 = var(i+3,j,k,vv);
                        f2 = var(i+2,j,k,vv);
                        f3 = var(i+1,j,k,vv);
                        f4 = var(i  ,j,k,vv);
                        f5 = var(i-1,j,k,vv);
                    }else{              
                        f1 = var(i-2,j,k,vv);
                        f2 = var(i-1,j,k,vv);
                        f3 = var(i  ,j,k,vv);
                        f4 = var(i+1,j,k,vv);
                        f5 = var(i+2,j,k,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wxr = ur*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                if (vv == 3){
                    if (vl < 0.0){
                        f1 = var(i,j+2,k,vv) + pres[10];
                        f2 = var(i,j+1,k,vv) + pres[4];
                        f3 = var(i,j  ,k,vv) + pres[0];
                        f4 = var(i,j-1,k,vv) + pres[3];
                        f5 = var(i,j-2,k,vv) + pres[9];
                    }else{              
                        f1 = var(i,j-3,k,vv) + pres[15];
                        f2 = var(i,j-2,k,vv) + pres[9];
                        f3 = var(i,j-1,k,vv) + pres[3];
                        f4 = var(i,j  ,k,vv) + pres[0];
                        f5 = var(i,j+1,k,vv) + pres[4];
                    }
                }else{
                    if (vl < 0.0){
                        f1 = var(i,j+2,k,vv);
                        f2 = var(i,j+1,k,vv);
                        f3 = var(i,j  ,k,vv);
                        f4 = var(i,j-1,k,vv);
                        f5 = var(i,j-2,k,vv);
                    }else{              
                        f1 = var(i,j-3,k,vv);
                        f2 = var(i,j-2,k,vv);
                        f3 = var(i,j-1,k,vv);
                        f4 = var(i,j  ,k,vv);
                        f5 = var(i,j+1,k,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wyl = vl*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                if (vv == 3){
                    if (vr < 0.0){
                        f1 = var(i,j+3,k,vv) + pres[16];
                        f2 = var(i,j+2,k,vv) + pres[10];
                        f3 = var(i,j+1,k,vv) + pres[4];
                        f4 = var(i,j  ,k,vv) + pres[0];
                        f5 = var(i,j-1,k,vv) + pres[3];
                    }else{              
                        f1 = var(i,j-2,k,vv) + pres[9];
                        f2 = var(i,j-1,k,vv) + pres[3];
                        f3 = var(i,j  ,k,vv) + pres[0];
                        f4 = var(i,j+1,k,vv) + pres[4];
                        f5 = var(i,j+2,k,vv) + pres[10];
                    }
                }else{
                    if (vr < 0.0){
                        f1 = var(i,j+3,k,vv);
                        f2 = var(i,j+2,k,vv);
                        f3 = var(i,j+1,k,vv);
                        f4 = var(i,j  ,k,vv);
                        f5 = var(i,j-1,k,vv);
                    }else{              
                        f1 = var(i,j-2,k,vv);
                        f2 = var(i,j-1,k,vv);
                        f3 = var(i,j  ,k,vv);
                        f4 = var(i,j+1,k,vv);
                        f5 = var(i,j+2,k,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wyr = vr*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                if (vv == 3){
                    if (wl < 0.0){
                        f1 = var(i,j,k+2,vv) + pres[12];
                        f2 = var(i,j,k+1,vv) + pres[6];
                        f3 = var(i,j,k  ,vv) + pres[0];
                        f4 = var(i,j,k-1,vv) + pres[5];
                        f5 = var(i,j,k-2,vv) + pres[11];
                    }else{              
                        f1 = var(i,j,k-3,vv) + pres[17];
                        f2 = var(i,j,k-2,vv) + pres[11];
                        f3 = var(i,j,k-1,vv) + pres[5];
                        f4 = var(i,j,k  ,vv) + pres[0];
                        f5 = var(i,j,k+1,vv) + pres[6];
                    }
                }else{
                    if (wl < 0.0){
                        f1 = var(i,j,k+2,vv);
                        f2 = var(i,j,k+1,vv);
                        f3 = var(i,j,k  ,vv);
                        f4 = var(i,j,k-1,vv);
                        f5 = var(i,j,k-2,vv);
                    }else{              
                        f1 = var(i,j,k-3,vv);
                        f2 = var(i,j,k-2,vv);
                        f3 = var(i,j,k-1,vv);
                        f4 = var(i,j,k  ,vv);
                        f5 = var(i,j,k+1,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wzl = wl*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                if (vv == 3){
                    if (wr < 0.0){
                        f1 = var(i,j,k+3,vv) + pres[18];
                        f2 = var(i,j,k+2,vv) + pres[12];
                        f3 = var(i,j,k+1,vv) + pres[6];
                        f4 = var(i,j,k  ,vv) + pres[0];
                        f5 = var(i,j,k-1,vv) + pres[5];
                    }else{              
                        f1 = var(i,j,k-2,vv) + pres[11];
                        f2 = var(i,j,k-1,vv) + pres[5];
                        f3 = var(i,j,k  ,vv) + pres[0];
                        f4 = var(i,j,k+1,vv) + pres[6];
                        f5 = var(i,j,k+2,vv) + pres[12];
                    }
                }else{
                    if (wr < 0.0){
                        f1 = var(i,j,k+3,vv);
                        f2 = var(i,j,k+2,vv);
                        f3 = var(i,j,k+1,vv);
                        f4 = var(i,j,k  ,vv);
                        f5 = var(i,j,k-1,vv);
                    }else{              
                        f1 = var(i,j,k-2,vv);
                        f2 = var(i,j,k-1,vv);
                        f3 = var(i,j,k  ,vv);
                        f4 = var(i,j,k+1,vv);
                        f5 = var(i,j,k+2,vv);
                    }
                }

                b1 = (13/12)*pow((f1-2.0*f2+f3),2.0) + (0.25)*pow((f1-4.0*f2+3.0*f3),2.0);
                b2 = (13/12)*pow((f2-2.0*f3+f4),2.0) + (0.25)*pow((f2-f4),2.0);
                b3 = (13/12)*pow((f3-2.0*f4+f5),2.0) + (0.25)*pow((3.0*f3-4.0*f4+f5),2.0);

                w1 = (0.1)/pow((eps+b1),2.0);
                w2 = (0.6)/pow((eps+b2),2.0);
                w3 = (0.3)/pow((eps+b3),2.0);

                p1 = ( 1.0/3.0)*f1 + (-7.0/6.0)*f2 + (11.0/6.0)*f3;
                p2 = (-1.0/6.0)*f2 + ( 5.0/6.0)*f3 + ( 1.0/3.0)*f4;
                p3 = ( 1.0/3.0)*f3 + ( 5.0/6.0)*f4 + (-1.0/6.0)*f5;

                wzr = wr*(w1*p1+w2*p2+w3*p3)/(w1+w2+w3);

                dvar(i,j,k,vv) = -( (wxr-wxl)/cd(1) + (wyr-wyl)/cd(2) +(wzr-wzl)/cd(3) );
                //dvar(i,j,k,vv) = -( (wxr-wxl)/cd(1) + (wyr-wyl)/cd(2)  );
                //dvar(i,j,k,vv) = -( (wxr-wxl)/cd(1) );
            }

            dvar(i,j,k,0) = dvar(i,j,k,0) - dxp;
            dvar(i,j,k,1) = dvar(i,j,k,1) - dyp;
            dvar(i,j,k,2) = dvar(i,j,k,2) - dzp;
            //dvar(i,j,k,1) = 0;
            //dvar(i,j,k,2) = 0;
        });

    }
};

void fnExit1(void){
    Kokkos::finalize();
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);
    Kokkos::initialize(argc, argv);

    atexit(fnExit1);

    int idx;

    struct inputConfig cf;

    cf = executeConfiguration(argv[1]);

    cf = mpi_init(cf);

    Kokkos::View<double****> myV("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> tmp("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> K1("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double****> K2("testView",cf.ngi,cf.ngj,cf.ngk,cf.nv);
    Kokkos::View<double*> cd("deviceCF",4+cf.ns*2);
    typename Kokkos::View<double*>::HostMirror hostcd = Kokkos::create_mirror_view(cd);
    Kokkos::deep_copy(hostcd, cd);
    hostcd(0) = cf.ns;
    hostcd(1) = cf.dx;
    hostcd(2) = cf.dy;
    hostcd(3) = cf.dz;
    
    int sdx = 4;
    for (int s=0; s<cf.ns; ++s){
        hostcd(sdx) = cf.gamma[s];
        hostcd(sdx+1) = cf.R/cf.M[s];
        sdx += 2;
    }
    Kokkos::deep_copy(cd,hostcd);

    if (cf.restart == 0)
        loadInitialConditions(cf,myV);
    
    if (cf.rank == 0){
        printf("loboSHOK - %s\n",cf.inputFname);
        printf("Restart File: %d | %s\n",cf.restart,cf.sfName);
        printf("-----------------------\n");
        printf("Running %d processes as (%d,%d,%d)\n",cf.numProcs,cf.xProcs,cf.yProcs,cf.zProcs);
        printf("nt = %d, dt = %f\n",cf.nt,cf.dt);
        printf("glbl_ni = %d, dx = %f\n",cf.glbl_ni,cf.dx);
        printf("glbl_nj = %d, dy = %f\n",cf.glbl_nj,cf.dy);
        printf("glbl_nk = %d, dz = %f\n",cf.glbl_nk,cf.dz);
        printf("Number of Species = %d:\n",cf.ns);
        for (int s=0; s<cf.ns; ++s)
            printf("    Species %d, Gamma = %4.2f, M = %6.4f\n",s+1,cf.gamma[s],cf.M[s]);
        printf("-----------------------\n");
    }
    

    /* allocate grid coordinate and flow variables */
    double *x = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *y = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    double *z = (double*)malloc(cf.ni*cf.nj*cf.nk*sizeof(double));
    

    /* calculate rank local grid coordinates */
    for (int k=0; k<cf.nk; ++k){
        for (int j=0; j<cf.nj; ++j){
            for (int i=0; i<cf.ni; ++i){
                idx = (cf.ni*cf.nj)*k+cf.ni*j+i;
                x[idx] = (cf.iStart + i)*cf.dx;
                y[idx] = (cf.jStart + j)*cf.dy;
                z[idx] = (cf.kStart + k)*cf.dz;
            }
        }
    }
    
    //typedef typename Kokkos::View<double****> ViewType;
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy_1;
    //Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy_1({0,0,0},{cf.nci, cf.ncj, cf.nck});
    
    double time = cf.time;
    int tstart = cf.tstart;
    
    MPI_Barrier(cf.comm);

    if (cf.restart == 1){
        readSolution(cf,myV);
    }else{
        writeSolution(cf,x,y,z,myV,0,0.00);
    }

    for (int t=tstart; t<cf.nt; ++t){
        time = time + cf.dt;

//        applyBCs(cf,myV);
//        rhs_func f1(cf.dt,cf,myV,K1, cd);
//        f1();
//
//        //myV = myV + K2
//        Kokkos::parallel_for("Loop2", policy_1({cf.ng,cf.ng,cf.ng,0},{cf.ngi-cf.ng,cf.ngj-cf.ng,cf.ngk-cf.ng,cf.nv}),
//               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
//            myV(i,j,k,v) = myV(i,j,k,v) + cf.dt*K1(i,j,k,v);
//        });

        //tmp = myV
        Kokkos::deep_copy(tmp,myV);

        //K1 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        rhs_func f1(cf.dt,cf,tmp,K1, cd);
        f1();
        
        //tmp = myV + k1/2
        Kokkos::parallel_for("Loop1", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            tmp(i,j,k,v) = myV(i,j,k,v) + cf.dt*K1(i,j,k,v)/2;
        });

        //K2 = dt*f(tmp) dt is member data to f()
        applyBCs(cf,tmp);
        rhs_func f2(cf.dt,cf,tmp,K2, cd);
        f2();

        //myV = myV + K2
        Kokkos::parallel_for("Loop2", policy_1({0,0,0,0},{cf.ngi, cf.ngj, cf.ngk, cf.nv}),
               KOKKOS_LAMBDA __device__ (const int i, const int j, const int k, const int v) {
            myV(i,j,k,v) = myV(i,j,k,v) + cf.dt*K2(i,j,k,v);
        });
        
        if (cf.rank==0){
            if ((t+1) % cf.out_freq == 0)
                printf("%d/%d, %f\n",t+1,cf.nt,time);
        }
        if ((t+1) % cf.write_freq == 0)
            writeSolution(cf,x,y,z,myV,t+1,time);
    }

    //if (cf.rank==0) printf("\n");
    //MPI_Barrier(cf.comm);
    //MYDBGR
    
    MPI_Finalize();
    //Kokkos::finalize();
    return 0;
}
