title = "3D Idealized Expansion"

--Restart and Output Options
restart_path = "./restarts"  --Path where hdf, xmf and vtk solution files will be written
progress_frequency = 0                          --Screen Output Interval
restart_frequency = 0                     --Restart Write Interval
status_frequency = 100

--restart = 1                           --Whether or not to use restart file
--time = 0.0                            --Start time of simulation
--tstart = 0                            --Start time index of simulation
--restart_name = "./restarts/restart-20.h5"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
gamma = {1.400}                  --Array of Species Ratio of Specific Heats
viscosity="on"
advection_scheme="weno5"
grid_type="cartesian"
buoyancy="on"

species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

blocks =  {
    {name="sol",path="./solution",frequency=100},
    {name="center", path="./centerline", frequency=100,
     start={1,50,1},limit={99,50,99},stride={1,1,1}}
}

--Time
nt = 10000
dt = 1e-5

ndim = 3
Lx,Ly,Lz = 2.0, 2.0, 2.0
ni       = {100, 100, 100}
dx       = {Lx/ni[1], Ly/ni[2], Lz/ni[3]}
procs    = {1, 1, 1}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
bcXmin = 1
bcXmax = 1
bcYmin = 1
bcYmax = 1
bcZmin = 0
bcZmax = 0

xPer = 0
yPer = 0
zPer = 0

cequations = 1             --Enable/Disable Cequation
kappa=10
epsilon=0.4
alpha=20
beta=0
betae=0

noise = 0

function grid(i,j,k,v)
    if v==0 then return dx[1]*i end
    if v==1 then return dx[2]*j end
    if v==2 then return dx[3]*k end
end

function initial_conditions(x,y,z,v)

    --Calcualte heat capacity of ambient air and hot air
    local Cv = R/(species[1].M*(species[1].gamma-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 30*Cv*10.0

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 300*Cv*1.0

    if z < zheight(x-1,y-1)+1 then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Ehot   end
        if v==4 then return 1.0    end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eair   end
        if v==4 then return 10.0    end
    end
end

function zheight(x,y)
    local A = 0.15*math.exp(-16*((x-0.25)^2 + y^2))
    local B = 0.15*math.exp(-16*((x+0.25)^2 + y^2))
    return A+B
end
