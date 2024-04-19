fiesta.title = "3D Idealized Expansion"

fiesta.status.frequency = 1000
fiesta.restart.auto = true

fiesta.eos.R = 8.314462
fiesta.advection_scheme="weno5"

fiesta.species={
    {name="Air",gamma=1.4,M=0.02897,mu=2.928e-5}
}

fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=1000},
    {name="center", path="./centerline", frequency=100,
     start={11,11,24},limit={34,34,24},stride={2,2,1}}
}

fiesta.time.tend = 5000
fiesta.time.dt = 1e-8

Lx,Ly,Lz = 10.0, 10.0, 10.0
fiesta.grid.type="cartesian"
fiesta.grid.ndim = 3
fiesta.grid.ni       = {50, 50, 50}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2], Lz/fiesta.grid.ni[3]}

fiesta.mpi.procs    = {1, 1, 1}
fiesta.mpi.type = "gpu-aware"

fiesta.bc.xmin = "reflective"
fiesta.bc.xmax = "outflow"
fiesta.bc.ymin = "reflective"
fiesta.bc.ymax = "outflow"
fiesta.bc.zmin = "reflective"
fiesta.bc.zmax = "outflow"

function fiesta.initial_conditions(x,y,z,v)
    --find position and distance from center
    local xdist = math.abs(Lx/2 -( x ))
    local ydist = math.abs(Ly/2 -( y ))
    local zdist = math.abs(Lz/2 -( z ))
    local rdist = math.sqrt(xdist*xdist + ydist*ydist + zdist*zdist)   --Radial distance from center

    --Calcualte heat capacity of ambient air and hot air
    local Cv = fiesta.eos.R/(fiesta.species[1].M*(fiesta.species[1].gamma-1))

    --Calcualte Internal energies
    Eair = 300 *Cv*1.0
    Ehot = 3000*Cv*1000.0

    if rdist < 0.5 then
        --setup hot region
        if v==3 then return Ehot   end
        if v==4 then return 1000.0 end
    else
        --ambient conditions
        if v==3 then return Eair   end
        if v==4 then return 1.0000 end
    end
    return 0
end
