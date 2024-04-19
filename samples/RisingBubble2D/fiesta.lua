frq = 1000
tt = 200000

title = "2D Plume Test"

--Restart and Output Options
out_freq = 0                          --Screen Output Interval
restart_freq = 0                      --Restart Write Interval
write_freq = frq                       --Solution Write Interval
stat_freq = frq
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
species_names = {'Air'}
gamma = {1.40}                        --Array of Species Ratio of Specific Heats
M = {0.02897}                         --Array of Species Molar Masses [kg/mol]
mu = {2.928e-5}
visc=0
scheme="weno5"
grid="cartesian"

buoyancy = 1

--Time
nt = tt                                 --Time Step at which to end simulation
dt = 1e-5                             --Time Step Size [s]

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
bcXmin = 0
bcXmax = 0
bcYmin = 1
bcYmax = 0

--Periodic Boundary Conditions
xPer = 0
yPer = 0

ceq = 0
noise = 0
n_dh = 0.1
n_eta = 0.1
n_coff = 0.0
n_nt = 1
particle = 0

function g(i,j,k,v)
    if v==0 then return dx*i end
    if v==1 then return dy*j end
end

function f(i,j,k,v)
    
    local xcenter = Lx/2    -- x coordinate of expansion center
    local ycenter = Ly/2 -- y coordinate of expansion center
    
    --find position and distance from expansion center
    local xdist = math.abs(xcenter -( (dx/2)+dx*i ))         --X distance from center
    local ydist = math.abs(ycenter -( (dy/2)+dy*j ))         --Y distance rom center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv = R/(M[1]*(gamma[1]-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eamb = 300 *Cv

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 350*Cv

    Ramb = 1.0
    Rhot = Eamb/Ehot

    if rdist < 1.0  then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Rhot*Ehot   end
        if v==3 then return Rhot   end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Ramb*Eamb   end
        if v==3 then return Ramb   end
    end
end
