--
-- 2D Ideal Expansion

title = "2D Idealized Expansion"

--Restart and Output Options
out_freq = 0                          --Screen Output Interval
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
species_names={"Air"}
gamma = {1.40}                        --Array of Species Ratio of Specific Heats
M = {0.02897}                         --Array of Species Molar Masses [kg/mol]
mu = {2.928e-5}
visc=0
scheme="weno5"
grid="cartesian"
buoyancy = 0

--Time
nt = 10                               --Time Step at which to end simulation
dt = 1e-6                             --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0

--Number of cells
ndim = 2
ni = 200            --Simulation size in x direction
nj = 200            --Simulation Size in y direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj

--MPI Processors
procsx = 1
procsy = 1

--Boundary Conditions (0: Freeflow, 1: Reflective)
bcXmin = 1
bcXmax = 0
bcYmin = 0
bcYmax = 0

--Periodic Boundary Conditions
xPer = 0
yPer = 0

ceq = 0
noise = 0
particle = 0

function g(i,j,k,v)
    if v==0 then return dx*i end
    if v==1 then return dy*j end
end

function f(x,y,z,v)
    --[[
    
    Apply initial conditions:  this function is called for every cell (but not ghost cells) and for each flow variable.
    Flow variables are indexed as follows:

    x-momentum:        0
    y-momentum:        1
    energy:            2
    species density 1: 3

    --]]
    

    ---2D ideal expansion---
    
    local xcenter = 0    -- x coordinate of expansion center
    local ycenter = Ly/2 -- y coordinate of expansion center
    
    --find position and distance from expansion center
    local xdist = math.abs(xcenter -( x ))         --X distance from center
    local ydist = math.abs(ycenter -( y ))         --Y distance rom center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv_air = R/(M[1]*(gamma[1]-1))
    local Cv_hot = R/(M[1]*(gamma[1]-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 300 *Cv_air*1.0

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 3000*Cv_hot*1000.0

    if rdist < 0.2 then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Ehot   end
        if v==3 then return 1000.0 end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Eair   end
        if v==3 then return 1.0000 end
    end
end
