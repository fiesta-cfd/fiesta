-- Input File

--Restart and Output Options
out_freq = 5                          --Screen Output Interval
restart_freq = 1000                   --Restart Write Interval
write_freq = 10                       --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                            --Start time of simulation
tstart = 0                            --Start time index of simulation
restartName = "restart-000000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 2                                --Number of Gas Species
gamma = {1.40, 1.60}                  --Array of Species Ration of Specifi Heats
M = {0.02897,0.14606}                 --Array of Species Molar Masses [kg/mol]

--Time
nt = 100                              --Time Step at which to end simulation
dt = 0.001                            --Time Step Size [s]

--Number of Ghost Points (Don't Change)
ng = 3

--User Parameters
Lx = 3.0
Ly = 1.5
Lz = 0.7

--Number of cells
ni = 300            --Simulation size in x direction
nj = 150            --Simulation Size in y direction
nk = 70             --Simulation Size in z direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj
dz = Lz/nk

--MPI Processors
procsx = 2
procsy = 1
procsz = 1

--Boundary Conditions
xPer = 0
yPer = 0
zPer = 0

-- C-Equation Coefficients
ceq = 1             --Enable/Disable Cequation
kappa = 10.0        --Smoothness Factor
epsilon = 1.0       --Support Factor
alpha = 1.0         --Anisotropic Coefficient
beta = 15.0         --Isotropic Coefficient
betae = 1.0         --Energy Equation Isotropic Coefficient

function f(i,j,k,v)
    local angle = 20
    local topdist = (Ly - ((dy/2)+dy*j))
    local xcenter = 0.6 + topdist*math.tan(angle*math.pi/180)
    local xdist = math.abs(xcenter-((dx/2)+dx*i))
    local zdist = math.abs(Lz/2-((dz/2)+dz*k))
    local rdist = math.sqrt(xdist*xdist + zdist*zdist)

    local xabs = (dx/2)+dx*i

    Eref = 0.9478
    Rref = 0.4555
    if xabs < 0.2 then
        if v==0 then return 1      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 3.5329 end
        if v==4 then return 1      end
        if v==5 then return 0      end
    end

    if rdist < 0.3 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return energy(rdist,Eref)   end
        if v==4 then return density1(rdist,Rref) end
        if v==5 then return density2(rdist,Rref) end
    else
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 0.9478 end
        if v==4 then return 0.4555 end
        if v==5 then return 0      end
    end
end

function density1(rad,Rref)
    mfs = mf(rad)
    M_mix = mfs*M[2]+(1-mfs)*M[1];
    return (1-mfs)*Rref*(M_mix/M[1])
end

function density2(rad,Rref)
    mfs = mf(rad)
    M_mix = mfs*M[2]+(1-mfs)*M[1];
    return mfs*Rref*(M_mix/M[1])
end

function energy(rad,Eref)
    mfs = mf(rad)
    Cps = (gamma[2]*R)/(M[2]*(gamma[2]-1))
    Cpa = (gamma[1]*R)/(M[1]*(gamma[1]-1))
    Cvs = R/(M[2]*(gamma[2]-1))
    Cva = R/(M[1]*(gamma[1]-1))
    gamma_mix = (mfs*Cps + (1-mfs)*Cpa)/(mfs*Cvs+(1-mfs)*Cva)
    return ((gamma[1]-1)/(gamma_mix-1))*Eref
end

function mf(rad)
    x = rad/0.150
    return math.exp(-math.pi*x*x)
end
