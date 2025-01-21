export blayer_mesh

function blayer_mesh(fname, H, L1, L2, dy0, dy1, dx0, dx1)

    
    y1 = wire_discr_dx(0.0, H; dx0=dy0, dx1=dy1)

    ny = length(y1)-1

    
    nosymm =  (L1 == 0.0) ? true : false

    if nosymm
        x1a = Float64[]
        x1b = wire_discr_dx(0.0, L2; dx0=dx0, dx1=dx1)
        nxa = 0
        nxb = length(x1b)-1
        x1 = x1b
    else
        x1a = wire_discr_dx(-L1, 0.0; dx0=dx1, dx1=dx0)
        x1b = wire_discr_dx(0.0, L2; dx0=dx0, dx1=dx1)
        nxa = length(x1a)-1
        nxb = length(x1b)-1
        x1 = join_wires(x1a, x1b)
    end
    nx = nxa+nxb
    npts = (nx+1)*(ny+1)
    open(fname, "w") do io
	println(io, "NDIME= 2")
	println(io, "NPOIN= $npts")
	cnt = 0
	for y in y1
	    for x in x1
		println(io, "$x\t$y\t$cnt")
		cnt += 1
	    end
	end
	nel = nx*ny
	cnt = 0
	println(io, "NELEM= $nel")
	for iy in 0:ny-1
	    for ix in 0:nx-1
		e1 = iy*(nx+1) + ix
		e2 = e1 + 1
		e3 = e2 + (nx+1)
		e4 = e3-1
		println(io, "9\t$e1\t$e2\t$e3\t$e4\t$cnt")
	    end
	end
	nmark = nosymm ? 4 : 5
	# Write boundary groups
	println(io, "NMARK= $nmark")
	println(io, "MARKER_TAG= farfield")
	println(io, "MARKER_ELEMS= $nx")
	for i in 0:(nx-1)
	    i1 = ny*(nx+1) + i
	    i2 = i1 + 1
	    println(io, "3\t$i1\t$i2")
	end
	println(io, "MARKER_TAG= inlet")
	println(io, "MARKER_ELEMS= $ny")
	for i in 0:ny-1
	    i1 = i*(nx+1)
	    i2 = i1 + nx + 1
	    println(io, "3\t$i1\t$i2")
	end
	println(io, "MARKER_TAG= outlet")
	println(io, "MARKER_ELEMS= $ny")
	for i in 0:ny-1
	    i1 = (i+1)*(nx+1) - 1
	    i2 = i1 + nx + 1
	    println(io, "3\t$i1\t$i2")
	end

        if !nosymm
	    println(io, "MARKER_TAG= symmetry")
	    println(io, "MARKER_ELEMS= $nxa")
	    for i in 0:nxa-1
	        i1 = i
	        i2 = i1 + 1
	        println(io, "3\t$i1\t$i2")
	    end
        end
	println(io, "MARKER_TAG= wall")
	println(io, "MARKER_ELEMS= $nxb")
	for i in 0:nxb-1
	    i1 = i+nxa
	    i2 = i1 + 1
	    println(io, "3\t$i1\t$i2")
	end
	
	
        
    end
    
    # Let's return the mesh
    symmidx = nosymm ? (1:1) : (1:nxa+1)
    wallidx = nosymm ? (1:nx+1) : (nxa:nx+1)

    # tensor
    idx = collect(reshape(1:npts, nx+1, ny+1))
    
        
    return (x=x1, y=y1, idx=idx, bcs=(inlet=(1,1:ny+1), outlet=(nx+1,1:ny+1), farfield=(1:nx,ny+1),
                                      wall=(wallidx, 1),symm=(symmidx, 1)))

end
