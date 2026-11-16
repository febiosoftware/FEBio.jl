module FEBio

using XML, Dates

export febiojl_dir, feb_doc_initialize, aen, read_logfile, run_febio
export showlog

# Export imported modules for later possible use
export XML

include("functions.jl")

end # module FEBio
