-- Input File

--Restart and Output Options
out_freq = 5                          --Screen Output Interval
restart_freq = 1000                   --Restart Write Interval
write_freq = 10                       --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 2                                --Number of Gas Species
gamma = {1.40, 1.60}                  --Array of Species Ration of Specifi Heats
M = {0.02897,0.14606}                 --Array of Species Molar Masses [kg/mol]

--Time
nt = 100                              --Time Step at which to end simulation
dt = 0.001                            --Time Step Size [s]

--Number of Ghost Points (Don't Change)
ng = 3

--User Parameters
Lx = 3.0
Ly = 1.5
Lz = 0.7

--Number of cells
ni = 300            --Simulation size in x direction
nj = 150            --Simulation Size in y direction
nk = 70             --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 2
procsy = 1
procsz = 1

--Boundary Conditions
xPer = 0
yPer = 0
zPer = 0

-- C-Equation Coefficients
ceq = 1             --Enable/Disable Cequation
kappa = 10.0        --Smoothness Factor
epsilon = 1.0       --Support Factor
alpha = 1.0         --Anisotropic Coefficient
beta = 15.0         --Isotropic Coefficient
betae = 1.0         --Energy Equation Isotropic Coefficient

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
    

    --Inclined gass cylinder example
    
    --find position and distance from cylinder axis
    local angle = 20
    local topdist = (Ly - ((dy/2)+dy*j))
    local xcenter = 0.6 + topdist*math.tan(angle*math.pi/180)
    local xdist = math.abs(xcenter-((dx/2)+dx*i))
    local zdist = math.abs(Lz/2-((dz/2)+dz*k))
    local rdist = math.sqrt(xdist*xdist + zdist*zdist)
    local xabs = (dx/2)+dx*i

    --setup shock wave
    if xabs < 0.2 then
        if v==0 then return 1      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 3.5329 end
        if v==4 then return 1      end
        if v==5 then return 0      end
    end

    --setup cylinder
    if rdist < 0.3 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 0.6319 end
        if v==4 then return 0      end
        if v==5 then return 1.8784 end
    --ambient conditions
    else
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 0.9478 end
        if v==4 then return 0.4555 end
        if v==5 then return 0      end
    end
end
