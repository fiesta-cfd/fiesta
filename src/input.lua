-- input file
out_freq = 1

gamma = 1.42
L = 1
nt = 10
ng = 3
ni = 100
nj = 100
nk = 100
nv = 6
dt = 0.1
dx = L/ni
dy = L/nj
dz = L/nk
procsx = 1
procsy = 1
procsz = 1

function f(i,j,k,v)
    local slant = 0.5 + (L-j*dy)/4
    local xdist = math.abs(slant-((dx/2)+dx*i))
    local zdist = math.abs(0.5-((dz/2)+dz*k))

    local rdist = math.sqrt(xdist*xdist + zdist*zdist)

    if rdist < 0.1 then
        return 1+v
    else
        return 0
    end
    
end
