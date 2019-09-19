-- input file
out_freq = 1

R = 8.314462 --J/(K*mol)
ns = 2
gamma = {1.41, 1.67}
M = {0.02897,0.14606} --kg/mol
L = 1
nt = 10
ng = 3
ni = 100
nj = 100
nk = 100
nv = 6
dt = 0.001
dx = L/ni
dy = L/nj
dz = L/nk
procsx = 1
procsy = 1
procsz = 1
xPer = 0
yPer = 0
zPer = 0

function f(i,j,k,v)
    local xdist = math.abs(0.5-((dx/2)+dx*i))
    local zdist = math.abs(0.5-((dz/2)+dz*k))
    local rdist = math.sqrt(xdist*xdist + zdist*zdist)

    local xabs = (dx/2)+dx*i

    if xabs < 0.2 then
        if v==0 then return 1      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 3.5329 end
        if v==4 then return 1      end
        if v==5 then return 0      end
    end

    if rdist < 0.1 then
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 0.9478 end
        if v==4 then return 0      end
        if v==5 then return 1.8784 end
    else
        if v==0 then return 0      end
        if v==1 then return 0      end
        if v==2 then return 0      end
        if v==3 then return 0.9478 end
        if v==4 then return 0.4555 end
        if v==5 then return 0      end
    end
end

--function f(i,j,k,v)
--    local slant = 0.5 + (L-j*dy)/4
--    local xdist = math.abs(slant-((dx/2)+dx*i))
--    local zdist = math.abs(0.5-((dz/2)+dz*k))
--
--    local rdist = math.sqrt(xdist*xdist + zdist*zdist)
--
--    if rdist < 0.1 then
--        return 1+v
--    else
--        return 0
--    end
--    
--end
