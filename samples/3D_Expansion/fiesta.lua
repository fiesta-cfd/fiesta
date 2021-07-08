title = "3D Idealized Expansion"

--Restart and Output Options
restart_path = "./restarts"  --Path where hdf, xmf and vtk solution files will be written
progress_frequency = 10                          --Screen Output Interval
restart_frequency = 0                     --Restart Write Interval
status_frequency = 10

--restart = 1                           --Whether or not to use restart file
--time = 0.0                            --Start time of simulation
--tstart = 0                            --Start time index of simulation
--restart_name = "./restarts/restart-20.h5"   --Restart File Name

--MPI Options
mpi = "gpu-aware"		      --Use CUDA-aware MPI

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
gamma = {1.400}                  --Array of Species Ratio of Specific Heats
viscosity="on"
advection_scheme="weno5"
grid_type="cartesian"
buoyancy="off"

species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

blocks =  {
    {name="sol",path="./solution",frequency=10},
    {name="center", path="./centerline", frequency=10,
    start={51,51,74},limit={149,149,74},stride={2,2,1}}
}

--Time
nt = 1000
dt = 1e-8

ndim = 3
Lx,Ly,Lz = 10.0, 10.0, 10.0
ni       = {150, 150, 150}
dx       = {Lx/ni[1], Ly/ni[2], Lz/ni[3]}
procs    = {2, 1, 2}

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

cequations = 0             --Enable/Disable Cequation
noise = 0

function grid(i,j,k,v)
    if v==0 then return dx[1]*i end
    if v==1 then return dx[2]*j end
    if v==2 then return dx[3]*k end
end

function initial_conditions(x,y,z,v)
    --  --[[
    --  
    --  Apply initial conditions:  this function is called for every cell (but not ghost cells) and for each flow variable.
    --  Flow variables are indexed as follows:

    --  x-momentum:        0
    --  y-momentum:        1
    --  z-momentum:        2
    --  energy:            3
    --  species density 1: 4
    --  species density 2: 5
    --  species density n: n+3

    --  --]]
    --  

    --  ---3D ideal expansion---
    --  
    --  --find position and distance from expansion center (0,0,0)
    --  local xdist = math.abs(0 -( (dx/2)+dx*i ))
    --  local ydist = math.abs(0 -( (dy/2)+dy*j ))
    --  local zdist = math.abs(0 -( (dz/2)+dz*k ))
    --  --find radial distance from center
    --  local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)

    --  --calcualte heat capacity of ambient and hot regions
    --  local Cv_air = R/(M[1]*(gamma[1]-1))
    --  local Cv_hot = R/(M[1]*(gamma[1]-1))
    --  
    --  --calculate internal energies
    --  --Air at 300[k] 1.0[kg/m^3]
    --  Eair = 300 *Cv_air*1.0
    --  --Air at 3000[k] 1000[kg/m^3]
    --  Ehot = 3000*Cv_hot*1000.0

    --  --setup hot region
    --  -- if rdist < 1.0 then
    --  --     if v==0 then return 0      end
    --  --     if v==1 then return 0      end
    --  --     if v==2 then return 0      end
    --  --     if v==3 then return Ehot   end
    --  --     if v==4 then return 1000.0 end
    --  --     if v==5 then return 0      end
    --  -- --ambient conditions
    --  -- else
    --      if v==0 then return i      end
    --      if v==1 then return j      end
    --      if v==2 then return k      end
    --      if v==3 then return Eair   end
    --      if v==4 then return 1.0000 end
    --      if v==5 then return 0      end
    --  -- end
    local xc = Lx/2 -- x coordinate of expansion center
    local yc = Ly/2 -- y coordinate of expansion center
    local zc = Lz/2 -- y coordinate of expansion center
    
    --find position and distance from expansion center
    local xdist = math.abs(xc -( x ))         --X distance from center
    local ydist = math.abs(yc -( y ))         --Y distance rom center
    local zdist = math.abs(zc -( z ))         --Y distance rom center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv = R/(species[1].M*(species[1].gamma-1))

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
