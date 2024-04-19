fiesta.title = "3D Idealized Expansion"

fiesta.status.frequency = 100
fiesta.progress.frequency = 0

fiesta.restart.auto = false
-- fiesta.restart.frequency = 100
-- fiesta.restart.path = './restarts'

fiesta.eos.R = 8.314462
fiesta.advection_scheme="weno5"
fiesta.viscosity.enabled = true

fiesta.species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=100,
     stride={1,1,1}}
     -- start={1,1,1},limit={50,50,24},stride={1,1,1}}
}

fiesta.time.tend = 50000
fiesta.time.dt = 1e-6

Lx,Ly,Lz = 0.1, 0.05, 0.01
fiesta.grid.type="cartesian"
fiesta.grid.ndim = 3
fiesta.grid.ni   = {100, 50, 10}
fiesta.grid.dx   = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2], Lz/fiesta.grid.ni[3]}

fiesta.mpi.procs    = {2,2,1}
fiesta.mpi.type = "gpu-aware"

fiesta.bc.xperiodic = true
-- fiesta.bc.xmin = "outflow"
-- fiesta.bc.xmax = "outflow"
fiesta.bc.ymin = "noslip"
fiesta.bc.ymax = "reflective"
-- fiesta.bc.zmin = "outflow"
-- fiesta.bc.zmax = "outflow"
fiesta.bc.zperiodic = true

function fiesta.initial_conditions(x,y,z,v)
    local rho = 1.0
    local vel = 10.0
    local temp = 300
    local Cv = fiesta.eos.R/(fiesta.species[1].M*(fiesta.species[1].gamma-1))

    Eair = rho*Cv*temp + 0.5*rho*vel*vel

    if v==0 then return rho*vel end
    if v==1 then return 0       end
    if v==2 then return 0       end
    if v==3 then return Eair    end
    if v==4 then return rho     end
    return 0
end
