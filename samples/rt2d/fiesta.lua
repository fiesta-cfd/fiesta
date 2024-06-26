fiesta.title = "Two Dimensional Rayleigh-Taylor"

--Output Options
fiesta.status.frequency = 5000

--Gas Properties
fiesta.R = 1.0                          --Universal Gas Constant [J/(K*mol)]

--Solvers and Physics
fiesta.advection_scheme="weno5"
fiesta.viscosity.enabled=false
fiesta.buoyancy.enabled=true
fiesta.buoyancy.acceleration=0.1

--Gas Species Data
fiesta.species={
    {name="Air",gamma=1.4,M=1.0,mu=2.928e-5}
}

--Solution IO Views
fiesta.ioviews =  {
    {name="sol",path="./solution",frequency=5000}
}

--Time
fiesta.time.dt = 1e-4
fiesta.time.nt = 7e4

--Grid
fiesta.grid.type="cartesian"
fiesta.grid.ndim = 2
Lx,Ly = 1.0/3.0, 1.0
fiesta.grid.ni       = {100, 200}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2]}

--mpi
fiesta.mpi.procs    = {1,1}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xperiodic=true
fiesta.bc.ymin = "hydrostatic"
fiesta.bc.ymax = "hydrostatic"

--Cequations
fiesta.ceq.enabled=true
fiesta.ceq.kappa=10.0
fiesta.ceq.epsilon=0.4
fiesta.ceq.alpha=40.0
fiesta.ceq.beta=0.0
fiesta.ceq.betae=0

--Noise Filter
fiesta.noise.enabled=true
fiesta.noise.dh=2e-5
fiesta.noise.eta=5e-4
fiesta.noise.coff=0.05
fiesta.noise.nt=1
fiesta.noise.mode=1

function fiesta.initial_conditions(x,y,z,v)
    h0=height(x-1.0/6.0)

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
        if v==2 then return eU end
        if v==3 then return rU end
    else
        if v==0 then return 0  end
        if v==1 then return 0  end
        if v==2 then return eL end
        if v==3 then return rL end
    end
end

function height(x)
    local A = 0.5+0.01*math.cos(6.0*math.pi*x)
    return A
end
