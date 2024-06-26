--
-- 3D Ideal Expansion
--

title = "3D Idealized Expansion"

--Restart and Output Options
out_freq = 5                          --Screen Output Interval
restart_freq = 0                      --Restart Write Interval
write_freq = 10                       --Solution Write Interval
stat_freq = 10
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
gamma = {1.400}                  --Array of Species Ratio of Specific Heats
M = {0.02897}                 --Array of Species Molar Masses [kg/mol]
mu = {2.928e-5}
visc=0
scheme="weno5"
grid="cartesian"
buoyancy=1

--Time
nt = 20                               --Time Step at which to end simulation
dt = 1e-6                             --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0
Lz = 10.0

--Number of cells
ndim = 3
ni = 100            --Simulation size in x direction
nj = 100            --Simulation Size in y direction
nk = 100            --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 1
procsy = 1
procsz = 1

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
bcXmin = 1
bcXmax = 0
bcYmin = 1
bcYmax = 0
bcZmin = 1
bcZmax = 0

xPer = 0
yPer = 0
zPer = 0

ceq = 0             --Enable/Disable Cequation
noise = 0
particle = 0

function g(i,j,k,v)
    if v==0 then return dx*i end
    if v==1 then return dy*j end
    if v==2 then return dz*k end
end

function f(i,j,k,v)
    --[[
    
    Apply initial conditions:  this function is called for every cell (but not ghost cells) and for each flow variable.
    Flow variables are indexed as follows:

    x-momentum:        0
    y-momentum:        1
    z-momentum:        2
    energy:            3
    species density 1: 4
    species density 2: 5
    species density n: n+3

    --]]
    

    ---3D ideal expansion---
    
    --find position and distance from expansion center (0,0,0)
    local xdist = math.abs(0 -( (dx/2)+dx*i ))
    local ydist = math.abs(0 -( (dy/2)+dy*j ))
    local zdist = math.abs(0 -( (dz/2)+dz*k ))
    --find radial distance from center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)

    --calcualte heat capacity of ambient and hot regions
    local Cv_air = R/(M[1]*(gamma[1]-1))
    local Cv_hot = R/(M[1]*(gamma[1]-1))
    
    --calculate internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 300 *Cv_air*1.0
    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 3000*Cv_hot*1000.0

    --setup hot region
    if rdist < 1.0 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Ehot   end
        if v==4 then return 1000.0 end
        if v==5 then return 0      end
    --ambient conditions
    else
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eair   end
        if v==4 then return 1.0000 end
        if v==5 then return 0      end
    end
end
