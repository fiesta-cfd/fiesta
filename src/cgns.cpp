#include "fiesta.hpp"
#include "cgns.hpp"
#include "debug.hpp"
#include <iostream>
#include <fstream>

struct inputConfig writeGrid(struct inputConfig cf, double *x, double *y, double *z,char * fname){

    //int icelldim, iphysdim, index_coord, index_field;
    cgsize_t isize[3][3];

    int Cx, Cy, Cz;
    /* calculate rank local start and end indices for CGNS file shapes */
    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t end[3] = {cf.iEnd+1,cf.jEnd+1,cf.kEnd+1};

    /* specify overall file shape for cgns */
    isize[0][0] = cf.glbl_ni;
    isize[0][1] = cf.glbl_nj;
    isize[0][2] = cf.glbl_nk;

    isize[1][0] = isize[0][0] - 1;
    isize[1][1] = isize[0][1] - 1;
    isize[1][2] = isize[0][2] - 1;
    
    isize[2][0] = 0;
    isize[2][1] = 0;
    isize[2][2] = 0;

    cgp_mpi_comm(cf.comm);
    
    /* create cgns file and write grid coordinates in parallel*/
    if (cgp_open(fname, CG_MODE_WRITE, &cf.cF) ||
        cg_base_write(cf.cF,"Base",3,3,&cf.cB) ||
        cg_zone_write(cf.cF,cf.cB,"Zone 1",*isize,CG_Structured,&cf.cZ))
        cgp_error_exit();

    if (cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealDouble, "CoordinateX", &Cx) ||
        cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealDouble, "CoordinateY", &Cy) ||
        cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealDouble, "CoordinateZ", &Cz))
        cgp_error_exit();

    if (cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cx, start, end, x) ||
        cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cy, start, end, y) ||
        cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cz, start, end, z))
        cgp_error_exit();

    cgp_close(cf.cF);

    return cf;
}
struct inputConfig writeSPGrid(struct inputConfig cf, float *x, float *y, float *z, char * fname){

    //int icelldim, iphysdim, index_coord, index_field;
    cgsize_t isize[3][3];

    int Cx, Cy, Cz;
    /* calculate rank local start and end indices for CGNS file shapes */
    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t end[3] = {cf.iEnd+1,cf.jEnd+1,cf.kEnd+1};

    /* specify overall file shape for cgns */
    isize[0][0] = cf.glbl_ni;
    isize[0][1] = cf.glbl_nj;
    isize[0][2] = cf.glbl_nk;

    isize[1][0] = isize[0][0] - 1;
    isize[1][1] = isize[0][1] - 1;
    isize[1][2] = isize[0][2] - 1;
    
    isize[2][0] = 0;
    isize[2][1] = 0;
    isize[2][2] = 0;

    cgp_mpi_comm(cf.comm);
    
    /* create cgns file and write grid coordinates in parallel*/
    if (cgp_open(fname, CG_MODE_WRITE, &cf.cF) ||
        cg_base_write(cf.cF,"Base",3,3,&cf.cB) ||
        cg_zone_write(cf.cF,cf.cB,"Zone 1",*isize,CG_Structured,&cf.cZ))
        cgp_error_exit();

    if (cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealSingle, "CoordinateX", &Cx) ||
        cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealSingle, "CoordinateY", &Cy) ||
        cgp_coord_write(cf.cF,cf.cB,cf.cZ,CG_RealSingle, "CoordinateZ", &Cz))
        cgp_error_exit();

    if (cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cx, start, end, x) ||
        cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cy, start, end, y) ||
        cgp_coord_write_data(cf.cF,cf.cB,cf.cZ,Cz, start, end, z))
        cgp_error_exit();

    cgp_close(cf.cF);

    return cf;
}

void writeRestart(struct inputConfig cf, double *x, double *y, double *z, const FS4D deviceV, int tdx, double time){
    int index_sol, index_flow;

    char dName[32];

    char solname[32];
    char fsname[32];
    cgsize_t dims = 1;
    cgsize_t idims[2];
    idims[0] = 32;
    idims[1] = 1;

    snprintf(fsname,32,"restart-%06d.cgns",tdx);
    if (cf.rank == 0){
        printf("Writing Restart %s...\n",fsname);
    }

    cf = writeGrid(cf, x,y,z,fsname);
    snprintf(solname,32,"FS");

    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t endc[3] = {cf.iEnd,cf.jEnd,cf.kEnd};

    FS4DH hostV = Kokkos::create_mirror_view(deviceV);
    Kokkos::deep_copy(hostV,deviceV);

    /* open cgns file and write cell centered flow variable */
    if (cgp_open(fsname, CG_MODE_MODIFY, &cf.cF))
        cgp_error_exit();
    
    if (cg_sol_write(cf.cF,cf.cB,cf.cZ,solname,CG_CellCenter, &index_sol))
        cgp_error_exit();

    double *v = (double*)malloc(cf.nci*cf.ncj*cf.nck*sizeof(double));

    int idx;

    //write momentum x
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,0);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,0);
            }
        }
    }

    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,"MomentumX",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write momentum y
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,1);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,1);
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,"MomentumY",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write momentum z
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,2);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = 0;
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,"MomentumZ",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write energy
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,cf.ndim);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,cf.ndim);
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,"EnergyInternal",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write densities
    for (int vn=cf.ndim+1; vn<cf.nv; ++vn){
        snprintf(dName,32,"SpeciesDensity%d",vn-cf.ndim);
        if (cf.ndim == 3){
            for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        int kk = k - cf.ng;
                        idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                        v[idx] = hostV(i,j,k,vn);
                    }
                }
            }
        }else{
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    idx = cf.nci*jj+ii;
                    v[idx] = hostV(i,j,0,vn);
                }
            }
        }
        
        if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,dName,&index_flow))
            cgp_error_exit();

        if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
            cgp_error_exit();
    }

    if (cf.ceq != 0){
        //write cequation variables
        for (int vn=cf.nv; vn<cf.nv+5; ++vn){
            snprintf(dName,32,"C%d",vn-cf.nv);
            if (cf.ndim == 3){
                for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                    for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                        for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                            int ii = i - cf.ng;
                            int jj = j - cf.ng;
                            int kk = k - cf.ng;
                            idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                            v[idx] = hostV(i,j,k,vn);
                        }
                    }
                }
            }else{
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        idx = cf.nci*jj+ii;
                        v[idx] = hostV(i,j,0,vn);
                    }
                }
            }
            
            if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealDouble,dName,&index_flow))
                cgp_error_exit();

            if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
                cgp_error_exit();
        }
    }

    cg_biter_write(cf.cF,cf.cB,"TimeIterValues",1);
    cg_goto(cf.cF,cf.cB,"BaseIterativeData_t",1,"end");
    cg_array_write("TimeValues",CG_RealDouble,1,&dims,&time);
    cg_ziter_write(cf.cF,cf.cB,cf.cZ,"ZoneIterativeData");
    cg_goto(cf.cF,cf.cB,"Zone_t",cf.cZ,"ZoneIterativeData_t",1,"end");
    cg_array_write("FlowSolutionPointers",CG_Character,2,idims,solname);
    cg_simulation_type_write(cf.cF,cf.cB,CG_TimeAccurate);
    
    cgp_close(cf.cF);
}

void writeSolution(struct inputConfig cf, float *x, float *y, float *z, const FS4D deviceV, int tdx, double time){
    int index_sol, index_flow;

    char dName[32];

    char solname[32];
    char fsname[32];
    cgsize_t dims = 1;
    cgsize_t idims[2];
    idims[0] = 32;
    idims[1] = 1;

    snprintf(fsname,32,"sol-%06d.cgns",tdx);
    if (cf.rank == 0){
        printf("Writing Solution %s...\n",fsname);
    }

    cf = writeSPGrid(cf, x,y,z,fsname);
    snprintf(solname,32,"FS");

    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t endc[3] = {cf.iEnd,cf.jEnd,cf.kEnd};

    FS4DH hostV = Kokkos::create_mirror_view(deviceV);
    Kokkos::deep_copy(hostV,deviceV);

    /* open cgns file and write cell centered flow variable */
    if (cgp_open(fsname, CG_MODE_MODIFY, &cf.cF))
        cgp_error_exit();
    
    if (cg_sol_write(cf.cF,cf.cB,cf.cZ,solname,CG_CellCenter, &index_sol))
        cgp_error_exit();

    float *v = (float*)malloc(cf.nci*cf.ncj*cf.nck*sizeof(float));

    int idx;

    //write momentum x
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,0);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,0);
            }
        }
    }

    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,"MomentumX",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write momentum y
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,1);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,1);
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,"MomentumY",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write momentum z
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,2);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = 0;
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,"MomentumZ",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write energy
    if (cf.ndim == 3){
        for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    int kk = k - cf.ng;
                    idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                    v[idx] = hostV(i,j,k,cf.ndim);
                }
            }
        }
    }else{
        for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
            for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                int ii = i - cf.ng;
                int jj = j - cf.ng;
                idx = cf.nci*jj+ii;
                v[idx] = hostV(i,j,0,cf.ndim);
            }
        }
    }
    
    if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,"EnergyInternal",&index_flow))
        cgp_error_exit();

    if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
        cgp_error_exit();

    //write densities
    for (int vn=cf.ndim+1; vn<cf.nv; ++vn){
        snprintf(dName,32,"SpeciesDensity%d",vn-cf.ndim);
        if (cf.ndim == 3){
            for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        int kk = k - cf.ng;
                        idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                        v[idx] = hostV(i,j,k,vn);
                    }
                }
            }
        }else{
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    idx = cf.nci*jj+ii;
                    v[idx] = hostV(i,j,0,vn);
                }
            }
        }
        
        if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,dName,&index_flow))
            cgp_error_exit();

        if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
            cgp_error_exit();
    }

    if (cf.ceq != 0){
    //write cequation variables
        for (int vn=cf.nv; vn<cf.nv+5; ++vn){
            snprintf(dName,32,"C%d",vn-cf.nv);
            if (cf.ndim == 3){
                for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                    for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                        for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                            int ii = i - cf.ng;
                            int jj = j - cf.ng;
                            int kk = k - cf.ng;
                            idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                            v[idx] = hostV(i,j,k,vn);
                        }
                    }
                }
            }else{
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        idx = cf.nci*jj+ii;
                        v[idx] = hostV(i,j,0,vn);
                    }
                }
            }
            
            if (cgp_field_write(cf.cF,cf.cB,cf.cZ,index_sol,CG_RealSingle,dName,&index_flow))
                cgp_error_exit();

            if (cgp_field_write_data(cf.cF,cf.cB,cf.cZ,index_sol,index_flow,start,endc,v))
                cgp_error_exit();
        }
    }

    cg_biter_write(cf.cF,cf.cB,"TimeIterValues",1);
    cg_goto(cf.cF,cf.cB,"BaseIterativeData_t",1,"end");
    cg_array_write("TimeValues",CG_RealDouble,1,&dims,&time);
    cg_ziter_write(cf.cF,cf.cB,cf.cZ,"ZoneIterativeData");
    cg_goto(cf.cF,cf.cB,"Zone_t",cf.cZ,"ZoneIterativeData_t",1,"end");
    cg_array_write("FlowSolutionPointers",CG_Character,2,idims,solname);
    cg_simulation_type_write(cf.cF,cf.cB,CG_TimeAccurate);
    
    cgp_close(cf.cF);

    if (cf.numProcs == 1){
        std::ofstream myfile;
        myfile.open("output.txt");
        for (int i=cf.ng; i<cf.nci; ++i){
           int  idx = (cf.nci*cf.ncj)*0+cf.nci*3+i;
//            myfile << x[idx] << ", " << hostV(i,3,0,0) << ", " << hostV(i,3,0,1) << ", " << hostV(i,3,0,2) << ", " << hostV(i,3,0,3) << std::endl;
            myfile << x[idx] << ", " << hostV(i,3,0,3) << std::endl;
        }
        myfile << std::endl;
        myfile.close();
    }

}

void readSolution(struct inputConfig cf, const FS4D deviceV){
    FS4DH hostV = Kokkos::create_mirror_view(deviceV);
    Kokkos::deep_copy(hostV,deviceV);

    /* open cgns file for reading restart information */
    if (cgp_open(cf.sfName, CG_MODE_MODIFY, &cf.cF))
        cgp_error_exit();

    int idx,vv;
    
    cgsize_t start[3] = {cf.iStart+1,cf.jStart+1,cf.kStart+1};
    cgsize_t endc[3] = {cf.iEnd,cf.jEnd,cf.kEnd};

    double * readV;
    readV = (double *)malloc(cf.nci*cf.ncj*cf.nck*cf.nv*sizeof(double));

    for (int v=0; v<cf.nv; ++v){
        if (cf.ndim == 2 && v > 1)
            vv = v+1;
        else
            vv = v;
        cgp_field_read_data(cf.cF,1,1,1,vv+1,start,endc,readV);
        if (cf.ndim == 3){
            for (int k=cf.ng; k<cf.nck+cf.ng; ++k){
                for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                    for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                        int ii = i - cf.ng;
                        int jj = j - cf.ng;
                        int kk = k - cf.ng;
                        idx = (cf.nci*cf.ncj)*kk+cf.nci*jj+ii;
                        hostV(i,j,k,v) = readV[idx];
                    }
                }
            }
        }else{
            for (int j=cf.ng; j<cf.ncj+cf.ng; ++j){
                for (int i=cf.ng; i<cf.nci+cf.ng; ++i){
                    int ii = i - cf.ng;
                    int jj = j - cf.ng;
                    idx = cf.nci*jj+ii;
                    hostV(i,j,0,v) = readV[idx];
                }
            }
        }
    }

    cgp_close(cf.cF);
    
    Kokkos::deep_copy(deviceV,hostV);
    free(readV);
}
