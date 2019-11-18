--
-- 2D (Thin 3D) Ideal Expansion

--Restart and Output Options
out_freq = 0                          --Screen Output Interval
restart_freq = 10000                  --Restart Write Interval
write_freq = 200                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
gamma = {1.40}                        --Array of Species Ratio of Specific Heats
M = {0.02897}                         --Array of Species Molar Masses [kg/mol]

--Time
nt = 4000                             --Time Step at which to end simulation
dt = 5e-7                             --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0

--Number of cells
ndim = 2
ni = 400            --Simulation size in x direction
nj = 400            --Simulation Size in y direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj

--MPI Processors
procsx = 2
procsy = 2

--Boundary Conditions (0: Freeflow, 1: Reflective)
bcXmin = 1
bcXmax = 0
bcYmin = 0
bcYmax = 0

--Periodic Boundary Conditions
xPer = 0
yPer = 0

-- C-Equation Coefficients
ceq = 0             --Enable/Disable Cequation
kappa = 10.0        --Smoothness Factor
epsilon = 1.0       --Support Factor
alpha = 10.0        --Anisotropic Coefficient
beta = 15.0         --Isotropic Coefficient
betae = 2.0         --Energy Equation Isotropic Coefficient

function f(i,j,k,v)
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
    local xdist = math.abs(xcenter -( (dx/2)+dx*i ))         --X distance from center
    local ydist = math.abs(ycenter -( (dy/2)+dy*j ))         --Y distance rom center
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
