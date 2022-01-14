local fmt = string.format
local function error(msg)
    print(fmt("ERROR:   %s",msg))
end
local function warning(msg)
    print(fmt("WARNING: %s",msg))
end
local function info(msg)
    print(fmt("INFO:    %s",msg))
end
fiesta = {}
fiesta.ceq = {}
fiesta.noise = {}
fiesta.restart = {}
fiesta.bc = {}
fiesta.time = {}
fiesta.viscosity = {}
fiesta.buoyancy = {}
fiesta.grid = {}
fiesta.mpi = {}
fiesta.hdf5 = {}
fiesta.eos = {}
fiesta.progress = {}
fiesta.status = {}
require("fiesta")

-- Check ioviews
for idx,value in ipairs(fiesta.ioviews) do
    local limit = {}
    local start = {}
    local stride = {}
    local name = value.name

    if (value.limit) then
        if(#value.limit ~= fiesta.grid.ndim) then
            error(fmt("Limits for '%s' needs %d elements for a %dd problem.",name,fiesta.grid.ndim,fiesta.grid.ndim))
        end
        limit = value.limit
    else
        warning("Default values will be used for '%s' limit indices.",name)
        limit = fiesta.grid.ni
    end

    if(value.start) then
        if(#value.start ~= fiesta.grid.ndim) then
            error(fmt("Limits for '%s' needs %d elements for a %dd problem.",name,fiesta.grid.ndim,fiesta.grid.ndim))
        end
        start = value.start
    else
        warning("Default values will be used for '%s' start indices.",name)
        start = {1,1,1}
    end

    if (value.stride) then
        if(#value.stride ~= fiesta.grid.ndim) then
            error(fmt("Limits for '%s' needs %d elements for a %dd problem.",name,fiesta.grid.ndim,fiesta.grid.ndim))
        end
        stride = value.stride
    else
        warning("Default values will be used for '%s' stride indices.",name)
        stride = {1,1,1}
    end

    for i=1,3 do
        if (start[i] > fiesta.grid.ni[i]) then
            error(fmt("Start index %d of '%s' is %d which is greater than the problem dimension of %d",i,name,start[i],fiesta.grid.ni[i]))
        end
        if (start[i] < 1) then
            error(fmt("Start index %d of '%s' is %d but must not be less than 1.",i,name,start[i]))
        end
        if (limit[i] > fiesta.grid.ni[i]) then
            error(fmt("limit index %d of '%s' is %d which is greater than the problem dimension of %d",i,name,limit[i],fiesta.grid.ni[i]))
        end
        if (limit[i] < 1) then
            error(fmt("Limit index %d of '%s' is %d but must not be less than 1.",i,name,limit[i]))
        end
        if (start[i] > limit[i]) then
            error(fmt("Start index %d of '%s' is %d which is greater than the limit %d",i,name,start[i],limit[i]))
        end
        if (stride[i] > (limit[i] - start[i] + 1)) then
            error(fmt("Stride index %d of '%s' is %d which is greater than limit-start for index %d",i,name,stride[i],i))
        end
    end
end

fiesta.initial_conditions(0,0,0,3)
if(string.lower(fiesta.grid.type)=="generalized") then
    fiesta.initialize_grid(0,0,0,0)
end

