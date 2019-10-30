-- Input File

--Restart and Output Options
out_freq = 10                         --Screen Output Interval
restart_freq = 2000                   --Restart Write Interval
write_freq = 100                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-002000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 2                                --Number of Gas Species
gamma = {1.402, 1.402}                  --Array of Species Ratio of Specific Heats
M = {0.028966,0.028966}                 --Array of Species Molar Masses [kg/mol]
--gamma = {1.402, 1.092}                  --Array of Species Ratio of Specific Heats
--M = {0.028966,0.14606}                 --Array of Species Molar Masses [kg/mol]

--Time
nt = 2000                            --Time Step at which to end simulation
dt = 0.0000001                        --Time Step Size [s]

--User Parameters
Lx = 0.40
Ly = 0.16
Lz = 0.004

--Number of cells
ni = 500            --Simulation size in x direction
nj = 200            --Simulation Size in y direction
nk = 5              --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 1
procsy = 1
procsz = 1

--Boundary Conditions (ignored if periodic conditions specified in that direction)
-- 0 Free Flow
-- 1 Reflective
bcXmin = 0
bcXmax = 0
bcYmin = 0
bcYmax = 0
bcZmin = 0
bcZmax = 0

xPer = 0
yPer = 0
zPer = 1

-- C-Equation Coefficients
ceq = 1             --Enable/Disable Cequation
kappa = 1.0         --Smoothness Factor
epsilon = 1.0      --Support Factor
alpha = 0.0        --Anisotropic Coefficient
beta = 0.0         --Isotropic Coefficient
betae = 0.0        --Energy Equation Isotropic Coefficient

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
    

    ---Quasi 2D inclined cylinder---
    
    --find position and distance from expansion center
    local angle = 20
    local topdist = (Ly - ((dy/2)+dy*j))
    local xcenter = 0.04 + topdist*math.tan(angle*math.pi/180)
    local xdist = math.abs(xcenter   -( (dx/2)+dx*i ))
    local xabs = (dx/2)+dx*i

    --Shocked
    local rho1 = 3.1334
    local temp1 = 507.19
    local vel1 = 433.97

    --Driven
    local rho2 = 1.1766
    local temp2 = 300

    --Test
    local rhoSF6 = 5.9332

    local Cv_air = R/(M[1]*(gamma[1]-1))
    local Cv_sf6 = R/(M[2]*(gamma[2]-1))
    Eshocked = temp1 *Cv_air*rho1 + 0.5*rho1*vel1*vel1
    Edriven = temp2 *Cv_air*rho2
    Esf6 = temp2 *Cv_sf6*rhoSF6

    if xabs < 0.02 then
        if v==0 then return rho1*vel1 end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eshocked   end
        if v==4 then return rho1   end
        if v==5 then return 0      end
    end

    --setup energy pill
    if xdist < 0.005 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Esf6   end
        if v==4 then return 0      end
        if v==5 then return rhoSF6 end
    --ambient conditions
    else
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Edriven end
        if v==4 then return rho2   end
        if v==5 then return 0      end
    end
end
