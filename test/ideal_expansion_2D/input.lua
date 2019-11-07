-- Input File

--Restart and Output Options
out_freq = 10                         --Screen Output Interval
restart_freq = 10000                  --Restart Write Interval
write_freq = 200                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
gamma = {1.40, 1.40}                  --Array of Species Ratio of Specific Heats
M = {0.02897,0.02897}                 --Array of Species Molar Masses [kg/mol]
--M = {0.02897,0.14606}                 --Array of Species Molar Masses [kg/mol]

--Time
nt = 10000                            --Time Step at which to end simulation
dt = 0.0000001                        --Time Step Size [s]

--User Parameters
Lx = 10.0
Ly = 10.0
Lz = 0.05

--Number of cells
ni = 1000           --Simulation size in x direction
nj = 1000           --Simulation Size in y direction
nk = 5              --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 2
procsy = 1
procsz = 1

--Boundary Conditions
bcXmin = 0
bcXmax = 0
bcYmin = 0
bcYmax = 0
bcZmin = 0
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
    

    ---Quasi 2D ideal expansion---
    
    --find position and distance from expansion center
    local xdist = math.abs(0   -( (dx/2)+dx*i ))
    local ydist = math.abs(Ly/2-( (dy/2)+dy*j ))
    local rdist = math.sqrt(xdist*xdist + ydist*ydist)

    local Cv_air = R/(M[1]*(gamma[1]-1))
    local Cv_exp = R/(M[1]*(gamma[1]-1))
    Eair = 300 *Cv_air*1.0
    Eexp = 3000*Cv_exp*1000.0

    --setup energy pill
    if rdist < 0.2 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eexp   end
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