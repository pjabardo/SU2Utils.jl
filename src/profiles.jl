export read_profile


function read_profile(io::IO)
    nmark = find_tag(io, "NMARK")

    prof = OrderedDict{String,Matrix{Float64}}()
    for i in 1:nmark
        s = readline(io)
        idx = findfirst("MARKER_TAG=", s)
        if isnothing(idx)
            error("Expected 'MARKER_TAG=' but not found.")
        end
        marker = strip(s[last(idx)+1:end])
        nrow = find_tag(io, "NROW")
        ncol = find_tag(io, "NCOL")
        tab = zeros(nrow, ncol)
        i = 1
        while true
            s = readline(io)
            idx = findfirst("#", s)
            if !isnothing(idx)
                s = s[begin:first(idx)-1]
            end
            if strip(s)==""
                continue
            end
            xrow = parse.(Float64, split(s))
            tab[i,:] .= xrow
            i = i + 1
            i > nrow && break
        end
        prof[marker] = tab
    end
    return prof
          
end

read_profile(fname) = open(read_profile, fname, "r")

