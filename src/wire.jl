#export Wire1d, zones, npoints, nsegments
export wire_discr, wire_discr_dx, join_wires


struct Wire1d
    x::Vector{Float64}
    zones::Vector{Int}
    zones_idx::Vector{UnitRange{Int}}
end

zones(w::Wire1d) = w.zones
npoints(w::Wire1d) = length(w.x)
nsegments(w::Wire1d) = npoints(w)-1

import Base.getindex
Base.getindex(w::Wire1d, idx) = w.x[idx]

#=
function join_wires(w...; zones=nothing)
    nw = length(w)
    if nw == 0
        error("At least a single wire should be here")
    elseif nw == 1
        return w[1]
    end
    n = npoints(w[1])
    zones_idx = [1:n]
    x = w[1].x[1:end] 
    xlast = x[end]
    idx_start = 1
    idx_last = n
    for wi in w[2:end]
        # The first point will be the same as the last point
        n = npoints(wi)
        x0 = wi[1]
        for i in 2:n
            xi = wi[i]
            push!(x, xi-x0+xlast)
        end
        idx_start = idx_last
        idx_last = idx_start + n - 1
        push!(zones_idx, idx_start:idx_last)
    end
   
        
end


function wire_discr_nd_nr(npts, r)
    if npts < 2
        error("A wire should have at least 2 points!")
    end
    if r <= 0
        error("r should be a positive number!")
    end
    
    n = npts-1
    if isapprox(r, 1.0)
        Δx = L / n
        r = 1.0
    else
        Δx = L * (r-1) / (r^n -  1)
    end

    return npts, r, Δx
end

=#
function wire_discr(x1, x2; npts=11, r=1.0, zone=1)
    
    L = x2 - x1

    n =  npts-1 # Number of segments
    if isapprox(r,1.0)
	Δx = L / n
        r = 1.0
    else
	Δx = L * (r-1) / (r^n - 1)
    end
        
    # Generate the points
    x = zeros(npts)
    x[1] = x1
    dx = Δx
    for i in 1:n
	x[i+1] = x[i] + dx
	dx *= r
    end
              
    x[end] = x2  # Just eliminating any floating points errors...
    return x
    #return Wire1d(x, [zone], [1:npts])
end

wire_discr(L; npts=11, r=1.0, zone=1) = wire_discr(0.0, L; npts=npts, r=r, zone=zone)
              

function wire_discr_dx(x1, x2; dx0, dx1, zone=1)

    L = x2 - x1
    r = (L-dx0) / (L-dx1)
    if r ≈ 1.0
	n = ceil(Int, L/dx0)
	Δx = L / n
    else
	nx = log(dx1/dx0) / log(r) + 1
	if nx < 1
	    n = 1
	else
	    n = ceil(Int, nx)
	end
	Δx = L*(r-1) / (r^n - 1)
    end
    
    x = zeros(n+1)
    x[1] = x1
    dx = Δx
    for i in 1:n
	x[i+1] = x[i] + dx
	dx *= r
    end
    x[end] = x2

    #return Wire1d(x, [zone], [1:(n+1)])
    return x
end    

wire_discr_dx(L; dx0, dx1, zone=1) = wire_discr_dx(0.0, L; dx0=dx0, dx1=dx1, zone=zone)


function join_wires(x1, x2)

    xout = x1[:]
    xlast = x1[end]
    xi = x2[begin]

    for x in x2[2:end]
        push!(xout, x-xi+xlast)
    end

    return xout
    
end
