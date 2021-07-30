fiesta.title = "2D Idealized Double Bubble"

--Restart and Output Options
fiesta.status_frequency = 4000
--fiesta.progress_frequency = 1

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
    {name="sol",path="./euler_nc_c_test",frequency=4000}
}

--Time
local time = 0.2
fiesta.time.dt = 0.5e-5
fiesta.time.nt = 200000

fiesta.grid.ndim = 2
Lx,Ly = 2.0, 2.0
fiesta.grid.ni       = {128, 128}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2]}
fiesta.mpi.procs    = {1, 1}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xmin = "reflective"
fiesta.bc.xmax = "reflective"
fiesta.bc.ymin = "hydrostatic"
fiesta.bc.ymax = "hydrostatic"

fiesta.ceq.enabled=true
fiesta.ceq.kappa=10
fiesta.ceq.epsilon=0.4
fiesta.ceq.alpha=20
fiesta.ceq.beta=0
fiesta.ceq.betae=0

fiesta.noise.enabled=true
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

    if y < zheight(x-1)+1 then
        --setup hot region
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Ehot   end
        if v==3 then return 0.5    end
    else
        --ambient conditions
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return Eair   end
        if v==3 then return 1.0    end
    end
end

function zheight(x)
    local A = 0.15*math.exp(-16*((x-0.25)^2))
    local B = 0.15*math.exp(-16*((x+0.25)^2))
    return A+B
end
