fiesta.title = "Quasi-3D Rayleigh-Taylor Problem"

--Restart and Output Options
fiesta.status_frequency = 100

--Gas Properties
fiesta.R = 1.0                          --Universal Gas Constant [J/(K*mol)]
fiesta.viscosity.enabled=false
fiesta.advection_scheme="weno5"
fiesta.grid.type="cartesian"
fiesta.buoyancy.enabled=true
fiesta.buoyancy.acceleration=0.1

fiesta.species={
    {name="Air",gamma=1.4,M=1.0,mu=2.928e-5}
}

fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=100}
}

--Time
fiesta.time.dt = 0.0004
fiesta.time.nt = 10000

fiesta.grid.ndim = 3
Lx,Ly,Lz = 1.0/20.0, 0.5, 1.0/3.0
fiesta.grid.ni       = {10, 100, 100}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2], Lz/fiesta.grid.ni[3]}
fiesta.mpi.procs    = {1,1,1}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xperiodic=true
--fiesta.bc.xmin = "reflective"
--fiesta.bc.xmax = "reflective"
fiesta.bc.ymin = "hydrostatic"
fiesta.bc.ymax = "hydrostatic"
--fiesta.bc.zperiodic=true
fiesta.bc.zmin = "reflective"
fiesta.bc.zmax = "reflective"

fiesta.ceq.enabled=true
fiesta.ceq.kappa=10
fiesta.ceq.epsilon=5
fiesta.ceq.alpha=20
fiesta.ceq.beta=0
fiesta.ceq.betae=0

fiesta.noise.enabled=true
fiesta.noise.dh=1e-5
fiesta.noise.eta=1e-3
fiesta.noise.coff=0.05
fiesta.noise.nt=1
fiesta.noise.mode=1

function fiesta.initial_conditions(x,y,z,v)
    h0=height(z-1.0/6.0)

    gamma = fiesta.species[1].gamma
    g = fiesta.buoyancy.acceleration

    rU = 2.0 --upper density
    rL = 1.0 --lower density

    P0=2.4 --"background" pressure
    pU = P0 + rU*g*(1-y)  --upper pressure
    pL = P0 + rL*g*(h0-y) + rU*g*(1-h0)  --lower pressure

    eU = pU/(gamma-1)  --upper energy
    eL = pL/(gamma-1)  --lower energy

    if y > h0 then
        if v==0 then return 0  end
        if v==1 then return 0  end
        if v==2 then return 0  end
        if v==3 then return eU end
        if v==4 then return rU end
    else
        if v==0 then return 0  end
        if v==1 then return 0  end
        if v==2 then return 0  end
        if v==3 then return eL end
        if v==4 then return rL end
    end
end

function height(x)
    local A = 0.25+0.01*math.cos(6.0*math.pi*x)
    return A
end
