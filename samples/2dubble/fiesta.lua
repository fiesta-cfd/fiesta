fiesta.title = "2D Double Bubble Problem"

--Restart and Output Options
fiesta.status_frequency = 1250

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
    {name="sol",path="./solution",frequency=500}
}

--Time
fiesta.time.dt = 0.001
fiesta.time.nt = 12000

fiesta.grid.ndim = 2
Lx,Ly = 2.0, 2.0
fiesta.grid.ni       = {128, 128}
fiesta.grid.dx       = {Lx/fiesta.grid.ni[1], Ly/fiesta.grid.ni[2]}
fiesta.mpi.procs    = {2, 2}

--Boundary Conditions (0 - Freeflow, 1 - Reflective)
fiesta.bc.xmin = "reflective"
fiesta.bc.xmax = "reflective"
fiesta.bc.ymin = "hydrostatic"
fiesta.bc.ymax = "hydrostatic"

fiesta.ceq.enabled=true
fiesta.ceq.kappa=10
fiesta.ceq.epsilon=0.4
fiesta.ceq.alpha=10
fiesta.ceq.beta=0
fiesta.ceq.betae=0

fiesta.noise.enabled=true
fiesta.noise.dh=2e-5
fiesta.noise.eta=5e-4
fiesta.noise.coff=0.05
fiesta.noise.nt=1
fiesta.noise.mode=1

function fiesta.initial_conditions(x,y,z,v)
    h0=height(x-1.0)

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
    local A = 0.5+0.15*math.exp(-16*(x-0.25)^2)
    local B = 0.5+0.15*math.exp(-16*(x+0.25)^2)
    return A+B
end
