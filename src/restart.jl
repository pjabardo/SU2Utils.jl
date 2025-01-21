export read_restart_file

const CGNS_STRING_SIZE = 33


function read_restart_file(fname)
    open(fname, "r") do io
        # Read header
        hdr = zeros(Int32,5)
        read!(io, hdr)

        if hdr[1] != 535532
            error("File $fname is not a binary SU2 file!")
        end
        nfields = hdr[2]
        npoints = hdr[3]

        fields = [strip(read(IOBuffer(read(io,CGNS_STRING_SIZE)), String), '\0') for i in 1:nfields]
        

        # Let's read the data!
        x = zeros(Float64, nfields, npoints)
        read!(io, x)
        
        fields, x
            
        
    end
    
end
