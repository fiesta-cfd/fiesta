--
-- 3D Ideal Expansion
--

--Restart and Output Options
out_freq = 10                         --Screen Output Interval
restart_freq = 10000                  --Restart Write Interval
write_freq = 100                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
gamma = {1.400}                  --Array of Species Ratio of Specific Heats
M = {0.02897}                 --Array of Species Molar Masses [kg/mol]
--M = {0.02897,0.14606}                 --Array of Species Molar Masses [kg/mol]

--Time
nt = 3000                             --Time Step at which to end simulation
dt = 1e-6                             --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0
Lz = 10.0

--Number of cells
ndim = 3
ni = 75             --Simulation size in x direction
nj = 75             --Simulation Size in y direction
nk = 75             --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 2
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
    if rdist < 0.5 then
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
