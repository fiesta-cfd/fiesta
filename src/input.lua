-- input file
gamma = 1.42
L = 1
ni = 10
nj = 10
nk = 10
dx = L/ni
dy = L/nj
dz = L/nk
procsx = 1
procsy = 1
procsz = 1

function f(x,y)
    --return(x^2 * math.sin(y))/(1-x)
    local rtn = (math.pi*x + math.pi*y)
    return rtn
    --return (x/y)
end
