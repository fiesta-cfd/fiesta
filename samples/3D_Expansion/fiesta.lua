--
-- 3D Ideal Expansion

title = "3D Idealized Expansion"

--Output Options
--pathName = "/users/beromer/scratch/3D_Test"  --Path where hdf, xmf and vtk solution files will be written
write_freq = 500                      --Solution Write Interval
stat_freq = 100

--Restart Options
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.h5"   --Restart File Name

--MPI Options
mpi = "gpu-aware"		      --Use CUDA-aware MPI

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
species_names={"Air"}
gamma = {1.409}                        --Array of Species Ratio of Specific Heats
M = {0.02897}                         --Array of Species Molar Masses [kg/mol]
mu = {2.928e-5}
visc=0
scheme="weno5"
grid="cartesian"
buoyancy = 0

--Time
nt = 10000                             --Time Step at which to end simulation
dt = 1e-7                             --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0
Lz = 10.0

--Number of cells
ndim = 3
ni = 1000           --Simulation size in x direction
nj = 1000           --Simulation Size in y direction
nk = 1000

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 4
procsy = 4
procsz = 4

--Boundary Conditions (0: Freeflow, 1: Reflective)
bcXmin = 1
bcXmax = 0
bcYmin = 0
bcYmax = 0
bcZmin = 0
bcZmax = 0

--Periodic Boundary Conditions
xPer = 0
yPer = 0
zPer = 0

ceq = 0
noise = 0
particle = 0

function g(i,j,k,v)
    if v==0 then return dx*i end
    if v==1 then return dy*j end
    if v==2 then return dz*k end
end

function f(x,y,z,v)
    --[[
    
    Apply initial conditions:  this function is called for every cell (but not ghost cells) and for each flow variable.
    Flow variables are indexed as follows:

    x-momentum:        0
    y-momentum:        1
    z-momentum:        2
    energy:            3
    species density 1: 4

    --]]
    

    ---2D ideal expansion---
    
    local xc = Lx/2 -- x coordinate of expansion center
    local yc = Ly/2 -- y coordinate of expansion center
    local zc = Lz/2 -- y coordinate of expansion center
    
    --find position and distance from expansion center
    local xdist = math.abs(xc -( x ))         --X distance from center
    local ydist = math.abs(yc -( y ))         --Y distance rom center
    local zdist = math.abs(zc -( z ))         --Y distance rom center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv = R/(M[1]*(gamma[1]-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 300 *Cv*1.0

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 3000*Cv*1000.0

    if rdist < 0.2 then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Ehot   end
        if v==4 then return 1000.0 end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eair   end
        if v==4 then return 1.0000 end
    end
end
