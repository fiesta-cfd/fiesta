--
-- 2D Ideal Expansion
d_fac = 4

--Restart and Output Options
out_freq = 0                          --Screen Output Interval
restart_freq = 0                      --Restart Write Interval
write_freq = 100                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                               --Start time of simulation
tstart = 0                           --Start time index of simulation
restartName = "restart-020000.cgns"   --Restart File Name

--Gas Properties
--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 2                                --Number of Gas Species
gamma = {1.402, 1.092}                  --Array of Species Ratio of Specific Heats
M = {0.028966,0.14606}                 --Array of Species Molar Masses [kg/mol]
visc = 0

--User Parameters
Ly = 0.1500
Lx = 0.40

--Number of cells
ndim = 2
ni = 2000           --Simulation size in x direction
nj = 750            --Simulation Size in y direction

--Cell Sizes
dx = Lx/ni
dy = Ly/nj

--Time
--dt = 0.6/(1000/dx+1000/dy)            --Time Step Size [s]
dt = 1e-8
--nt = math.ceil(0.0010/dt)             --Time Step at which to end simulation
nt = 1000

--MPI Processors
procsx = 1
procsy = 1

--Boundary Conditions (0: Freeflow, 1: Reflective)
bcXmin = 0
bcXmax = 0
bcYmin = 1
bcYmax = 1

--Periodic Boundary Conditions
xPer = 0
yPer = 0

-- C-Equation Coefficients
ceq = 0             --Enable/Disable Cequation
kappa = 10.0        --Smoothness Factor
epsilon = 1.0       --Support Factor
alpha = 10.0        --Anisotropic Coefficient
beta = 15.0         --Isotropic Coefficient
betae = 2.0         --Energy Equation Isotropic Coefficient

-- Problem Setup
gam = gamma[1]
M1 = 2.0
T1 = 300
P1 = 1e5

Rg = R/M[1]
Cv = Rg/(gam-1)
Cp = Rg + Cv
r1 = P1/(Rg*T1)

M2 = math.sqrt( (M1^2+(2/(gam-1))) / (2*gam/(gam-1)*M1^2-1) )
r_ratio = ((gam+1)*M1^2) / (((gam-1)*M1^2) + 2)
T_ratio = (1+((gam-1)/2)*M1^2) / (1+((gam-1)/2)*M2^2)

T2 = T_ratio*T1
r2 = r_ratio*r1

a1 = math.sqrt(gam*Rg*T1)
a2 = math.sqrt(gam*Rg*T2)

u1 = 0
us = M1*a1
u2p = M2*a2
u2 = us - u2p

e1 = r1*Cv*T1
e2 = r2*Cv*T2 + (r2/2)*u2^2

Cv_sf6 = R/(M[2]*(gamma[2]-1))

--print(M2)
--print(r1,r2,r_ratio)
--print(T1,T2,T_ratio)
--print(u1,us,u2p,u2,r2*u2)
--print(e1,e2)

--initial conditions
function f(i,j,k,v)
    local angle = 30
    local topdist = (Ly - ((dy/2)+dy*j))
    local xcenter = 0.04 + topdist*math.tan(angle*math.pi/180)
    local xdist = math.abs(xcenter   -( (dx/2)+dx*i ))
    local xabs = (dx/2)+dx*i

    --setup shock
    if xabs < 0.02 then
        if v==0 then return r2*u2 end
        if v==1 then return 0      end
        if v==2 then return e2   end
        if v==3 then return r2   end
        if v==4 then return 0      end
    end

    --setup cylinder
    if xdist < 3*0.0065 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return energy(xdist,e1) end
        if v==3 then return density1(xdist,r1) end
        if v==4 then return density2(xdist,r1) end
    --ambient conditions
    else
        if v==0 then return r1*u1 end
        if v==1 then return 0     end
        if v==2 then return e1    end
        if v==3 then return r1    end
        if v==4 then return 0     end
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
    x = rad/(2*0.0065/d_fac)
    return math.exp(-math.pi*x*x)
end
