#ifndef LSDEBUG_H
#define LSDEBUG_H
#define MYDBG printf("%s:%d\n",__FILE__,__LINE__);
#define MYDBGR printf("%d - %s:%d\n",cf.rank,__FILE__,__LINE__);
#define MYDBG0 if(cf.rank==0) printf("%s:%d\n",__FILE__,__LINE__);
enum cell_type
{
    BOUNDARY_CELL_ = 0,      // Ghost or Halo cell
    REAL_CELL_ = 1,          // Real cell
    BORDER_CELL_ = 2          // Border cell
};
/*
enum border_type
{
    LBORDER_CELL = 0,         // Real cell just inside the boundary cells. Used for Boundary conditions. Left
    RBORDER_CELL = 1,         // Real cell just inside the boundary cells. Used for Boundary conditions. Right
    BBORDER_CELL = 2,         // Real cell just inside the boundary cells. Used for Boundary conditions. Bot
    TBORDER_CELL = 3,         // Real cell just inside the boundary cells. Used for Boundary conditions. Top
    HBORDER_CELL = 4,         // Real cell just inside the boundary cells. Used for Boundary conditions. Back
    FBORDER_CELL = 5         // Real cell just inside the boundary cells. Used for Boundary conditions. Front
};
enum corner_type
{
    NOT_CORN = 0,
    LB_CORN = 1,
    LT_CORN = 2,
    RB_CORN = 3,
    RT_CORN = 4
    // 3D
    //NOT_CORN = 0,
    //LBH_CORN = 1,
    //LTH_CORN = 2,
    //RBH_CORN = 3,
    //RTH_CORN = 4
    //LBF_CORN = 5,
    //LTF_CORN = 6,
    //RBF_CORN = 7,
    //RTF_CORN = 8
};
*/
#endif
