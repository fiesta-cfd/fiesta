-- input file
gamma = 1.42
L = 1
ni = 10
nj = 10
nk = 10
dx = L/ni
dy = L/nj
dz = L/nk
procsx = 2
procsy = 2
procsz = 2

function f(x,y)
    return(x^2 * math.sin(y))/(1-x)
end
