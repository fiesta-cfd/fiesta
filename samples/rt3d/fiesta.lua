fiesta.title = "3D Idealized Expansion"

--Restart and Output Options
fiesta.status_frequency = 100

--restart = 1                           --Whether or not to use restart file
--time = 0.0                            --Start time of simulation
--tstart = 0                            --Start time index of simulation
--restart_name = "./restarts/restart-20.h5"   --Restart File Name

--Gas Properties
fiesta.R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
fiesta.viscosity.enabled=false
fiesta.advection_scheme="weno5"
fiesta.grid.type="cartesian"
fiesta.buoyancy.enabled=true

fiesta.species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=500},
    {name="center", path="./centerline", frequency=100,
     start={1,50,1},limit={99,50,99},stride={1,1,1}}
}

--Time
local time = 0.2
fiesta.time.dt = 0.5e-5
fiesta.time.nt = time//fiesta.time.dt

fiesta.grid.ndim = 3
Lx,Ly,Lz = 2.0, 2.0, 2.0
fiesta.grid.ni       = {100, 100, 100}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2], Lz/fiesta.grid.ni[3]}
fiesta.mpi.procs    = {1, 1, 1}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xmin = "reflective"
fiesta.bc.xmax = "reflective"
fiesta.bc.ymin = "reflective"
fiesta.bc.ymax = "reflective"
fiesta.bc.zmin = "hydrostatic"
fiesta.bc.zmax = "hydrostatic"

fiesta.ceq.enabled=false
fiesta.ceq.kappa=10
fiesta.ceq.epsilon=0.4
fiesta.ceq.alpha=20
fiesta.ceq.beta=0
fiesta.ceq.betae=0

fiesta.noise.enabled=false
fiesta.noise.dh=2e-5
fiesta.noise.eta=5e-4
fiesta.noise.coff=0.05
fiesta.noise.nt=1

function fiesta.initial_conditions(x,y,z,v)

    --Calcualte heat capacity of ambient air and hot air
    local Cv = fiesta.R/(fiesta.species[1].M*(fiesta.species[1].gamma-1))

    --Calcualte Internal energies
    --Air at 300[k] 1.0[kg/m^3]
    Eair = 300*Cv*1.0

    --Air at 3000[k] 1000[kg/m^3]
    Ehot = 600*Cv*0.5

    if z < zheight(x-1,y-1)+1 then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Ehot   end
        if v==4 then return 0.5    end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return Eair   end
        if v==4 then return 1.0    end
    end
end

function zheight(x,y)
    local A = 0.15*math.exp(-16*((x-0.25)^2 + y^2))
    local B = 0.15*math.exp(-16*((x+0.25)^2 + y^2))
    return A+B
end
