

export read_su2_mesh, msh2meshes, SU2Mesh

import DataStructures: OrderedDict
using Meshes




function find_tag(io::IO, tag)
    while !eof(io)
        s = readline(io)
        idx = findfirst("$(tag)=", uppercase(s))
        if !isnothing(idx)
            ntag = parse(Int, strip(s[last(idx)+1:end]))
            return ntag
        end
    end

    error("Could not find tag $tag in file")
end

find_tag(fname, tag) = open(fname, "r") do io
    find_tag(io, tag)
end

    
function read_ndime(fname)

    open(fname, "r") do io
        s = readline(io)
        idx = findfirst("NDIME=", s)
        if isnothing(idx)
            error("Substring 'NDIME=' not found in the first lines")
        end
        
        ndime = parse(Int, strip(s[last(idx)+1:end]))
        ndime
    end
    
end


function read_points(io::IO, ndime=2)
    npoin = find_tag(io, "NPOIN")
    # Read the points - x and y
    pts1 = zeros(npoin, ndime)
    for i in 1:npoin
        s = split(readline(io))
        for k in 1:ndime
            pts1[i,k] = parse(Float64, s[k])
        end
    end
    pts1
    
end
read_points(fname, ndime=2) = open(fname, "r") do io
    read_points(io, ndime)
end

struct ElemType
    id::Int
    name::Symbol
    nverts::Int
    ndim::Int
end

const elem_types = OrderedDict(3=>ElemType(3, :Line, 2, 1),
                               5=>ElemType(5, :Tri, 3, 2),
                               9=>ElemType(9, :Qua, 4, 2),
                               10=>ElemType(10, :Tet, 4, 3),
                               12=>ElemType(12, :Hex, 8, 3),
                               13=>ElemType(13, :Prism, 6, 3),
                               14=>ElemType(14, :Pyr, 5, 3))

const ConnType = Tuple{Int, Vararg{Int}}


function read_elements(io::IO; tag="NELEM", maxdim=3)
    
    nelem = find_tag(io, tag)
    conn = Vector{Int}[]
    etype = zeros(Int, nelem)
    kk = keys(elem_types)
    for e in 1:nelem
        ii = parse.(Int, split(readline(io)))
        nt = ii[1]
        etype[e] = nt
        
        if nt ∉ kk
            error("Element $e has element type $nt not implemented!")
        end
        ndim = elem_types[nt].ndim
        if ndim > maxdim
            error("Element $(e-1) has incompatible dimension $(ndim)!")
        end
        
        nverts = elem_types[nt].nverts
        
        # Remember: SU2 vertex ids start from 0. In julia, arrays are 1 based:
        push!(conn, ii[2:nverts+1] .+ 1)
        
    end
    
    etype, conn
end

read_elements(fname; tag="NELEM", maxdim=3) = open(fname, "r") do io
    read_elements(io, tag=tag, maxdim=maxdim)
end

struct SU2Mesh
    ndim::Int
    pts::Matrix{Float64}
    etype::Vector{Int}
    conn::Vector{Vector{Int}}
    bcs::OrderedDict{String,Tuple{Vector{Int},Vector{Vector{Int}}}}
end


    

function read_su2_mesh(fname)
    ndime = read_ndime(fname)

    points = read_points(fname, ndime)
    elems = read_elements(fname; tag="NELEM", maxdim=ndime)

    # Now lets read the BCs
    ConnType = Tuple{Int, Vararg{Int}}
    
    BCs = OrderedDict{String,Tuple{Vector{Int},Vector{Vector{Int}}}}()
    open(fname, "r") do io
        nmark = find_tag(io, "NMARK")
        for i in 1:nmark
            s = readline(io)
            idx = findfirst("MARKER_TAG=", s)
            tag = strip(s[last(idx)+1:end])
            # Read the elements
            etype, conn = read_elements(io; tag="MARKER_ELEMS", maxdim=ndime-1)
            BCs[tag] = (etype, conn)
        end
    end
    return SU2Mesh(ndime, points, elems[1], elems[2], BCs)
                                                                     
    #    return (ndim=ndime, points=points, elems=elems, bcs=BCs)
end

const vtkid_meshes = OrderedDict(3=>Line,
                                 5=>Triangle, 9=>Quadrangle,
                                 10=>Tetrahedron, 12=>Hexahedron, 13=>Wedge, 14=>Pyramid)

function meshes_connect(etype, conn)
    
    [connect((conn[i]...,), vtkid_meshes[etype[i]]) for i in eachindex(etype)]
    
end

#msh2meshes(pts, elems) = SimpleMesh(pts, meshes_connect(elems[1], elems[2]))

function msh2meshes(msh::SU2Mesh)

    ndim = msh.ndim

    pts = [Point(p...) for p in eachrow(msh.pts)]

    mesh = SimpleMesh(pts, meshes_connect(msh.etype, msh.conn))

    # Let's get the meshes of each BC
    if ndim==2
        bcs = OrderedDict(k=>[(e[1], e[2]) for e in v[2]] for (k,v) in msh.bcs)
    else
        bcs = OrderedDict(k=>SimpleMesh(pts, meshes_connect(v[1], v[2])) for (k,v) in msh.bcs)
    end
    
    return (mesh=mesh, bcs=bcs)
        
end
