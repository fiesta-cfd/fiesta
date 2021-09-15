fiesta.title = "3D Idealized Expansion"

--Restart and Output Options
fiesta.progress_frequency = 10
fiesta.status_frequency = 100

fiesta.restart.path = "./restarts"
fiesta.restart.frequency = 100
fiesta.restart.name = "restarts/restart-200.h5"
fiesta.restart.enabled = true
--fiesta.restart.reset = true
--fiesta.time.start_index = 0
--fiesta.time.start_time = 0e-6

--MPI Options
mpi = "gpu-aware"		      --Use CUDA-aware MPI

--Gas Properties
fiesta.eos.R = 8.314462
fiesta.advection_scheme="weno5"

fiesta.species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=100},
    {name="center", path="./centerline", frequency=100,
    start={11,11,24},limit={14,14,24},stride={2,2,1}}
    --start={3,4,9},limit={17,17,9},stride={2,2,1}}
}

--Time
fiesta.time.nt = 200
fiesta.time.dt = 1e-8

Lx,Ly,Lz = 10.0, 10.0, 10.0
fiesta.grid.type="cartesian"
fiesta.grid.ndim = 3
fiesta.grid.ni       = {50, 50, 50}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2], Lz/fiesta.grid.ni[3]}
fiesta.mpi.procs    = {2, 1, 1}
fiesta.mpi.type = "gpu-aware-ordered"

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xmin = "reflective"
fiesta.bc.xmax = "outflow"
fiesta.bc.ymin = "reflective"
fiesta.bc.ymax = "outflow"
fiesta.bc.zmin = "reflective"
fiesta.bc.zmax = "outflow"

function fiesta.initial_conditions(x,y,z,v)
    local xc = Lx/2 -- x coordinate of expansion center
    local yc = Ly/2 -- y coordinate of expansion center
    local zc = Lz/2 -- y coordinate of expansion center
    
    --find position and distance from expansion center
    local xdist = math.abs(xc -( x ))         --X distance from center
    local ydist = math.abs(yc -( y ))         --Y distance rom center
    local zdist = math.abs(zc -( z ))         --Y distance rom center
    local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv = fiesta.R/(fiesta.species[1].M*(fiesta.species[1].gamma-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 300 *Cv*1.0

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 3000*Cv*1000.0

    if rdist < 1.2 then
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
