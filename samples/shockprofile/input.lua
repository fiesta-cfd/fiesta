--
-- 2D Ideal Expansion

--Restart and Output Options
out_freq = 50                         --Screen Output Interval
restart_freq = 0                      --Restart Write Interval
write_freq = 500                      --Solution Write Interval
restart = 0                           --Whether or not to use restart file
time = 0.0                               --Start time of simulation
tstart = 0                           --Start time index of simulation
restartName = "restart-020000.cgns"   --Restart File Name

--Gas Properties
R = 8.314462                          --Universal Gas Constant [J/(K*mol)]
ns = 1                                --Number of Gas Species
gamma = {1.6667}                        --Array of Species Ratio of Specific Heats
M = {0.02897}                         --Array of Species Molar Masses [kg/mol]
visc = 1

--Time
nt = 36500
dt = 2e-11                            --Time Step Size [s]

--Number of cells
ndim = 2
ni = 10000          --Simulation size in x direction
nj = 10             --Simulation Size in y direction

--Cell Sizes
dx = 1e-5
dy = 1e-5

--MPI Processors
procsx = 1
procsy = 1

--Boundary Conditions (0: Freeflow, 1: Reflective)
bcXmin = 0
bcXmax = 0
bcYmin = 0
bcYmax = 0

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


-- problem statement

--gam = 5.0/3.0
--cp = 520
--rmachinit = 3.38
--tempinit = 300
--rhoinit = 0.0001068
--
--cv = cp/gam
--rg = cp-cv
--
--term = rmachinit*rmachinit
--
--densitytmp = (gam+1)/(gam-1+2/term)
--utmp = 2*(rmachinit^2-1)/((gam+1)*rmachinit)
--prtmp = 1+2*gam/(gam+1)*(rmachinit^2-1)
--ertmp = tempinit*cv*rhoinit
--
--cr = math.sqrt(gam*(gam-1)*ertmp)
--prb = rg*rhoinit*tempinit
--sound = math.sqrt(gam*prb/rhoinit)
--
--rke = 0.5*sound*utmp*densitytmp*rhoinit*sound*utmp
--temp = prb*prtmp/(rg*densitytmp*rhoinit)

gam = 5.0/3.0
Cp = 520
M1 = 3.38
T1 = 300
r1 = 0.0001068

Cv = Cp/gam
R  = Cp-Cv

M2 = math.sqrt( (M1^2+(2/(gam-1))) / (2*gam/(gam-1)*M1^2-1) )
p_ratio = 2*gam/(gam+1)*M1^2 - (gam-1)/(gam+1)
r_ratio = ((gam+1)*M1^2) / (((gam-1)*M1^2) + 2)
T_ratio = (1+((gam-1)/2)*M1^2) / (1+((gam-1)/2)*M2^2)

T2 = T_ratio*T1
r2 = r_ratio*r1

a1 = math.sqrt(gam*R*T1)
a2 = math.sqrt(gam*R*T2)

u1 = 0
us = M1*a1
u2p = M2*a2
u2 = us - u2p

e1 = r1*Cv*T1
e2 = r2*Cv*T2 + (r2/2)*u2^2

--initial conditions
function f(i,j,k,v)
    --find x and y position
    xdist = math.abs((dx/2)+dx*i)
    ydist = math.abs((dy/2)+dy*j)
    
    if xdist < 0.02 then
        --setup hot region
        if v==0 then return r2*u2 end
        if v==1 then return 0      end
        if v==2 then return e2 end
        if v==3 then return r2   end
    else
        --ambient conditions
        if v==0 then return r2*u1     end
        if v==1 then return 0      end
        if v==2 then return e1 end
        if v==3 then return r1 end
    end
end
