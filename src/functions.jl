# Call required packages
using XML

# Functions 

function febiojl_dir()
    pkgdir(@__MODULE__)
end

function feb_doc_initialize()
    # Define document level
    doc = XML.Document()

    # Define declaration: <?xml version="1.0" encoding="UTF-8"?>
    declarationNode = XML.Declaration(version = "1.0", encoding = "UTF-8")
    push!(doc, declarationNode)

    # Define febio_spec node: <febio_spec version="4.0">
    febio_spec_node = XML.Element("febio_spec"; version = "4.0")
    push!(doc, febio_spec_node)

    # Add comment: <!--Created using Julia-->
    comment_node = XML.Comment("Created using Julia")
    push!(febio_spec_node, comment_node)
    
    return doc,febio_spec_node
end

function aen(main_node::Node,sub_node_name::String, args...; kwargs...)
    if isnothing(args)
        sub_node = XML.Element(sub_node_name; kwargs...)    
    else
        sub_node = XML.Element(sub_node_name, args...; kwargs...)    
    end
    push!(main_node, sub_node)
    return sub_node
end

struct stepData
    data::Vector{Vector{Float64}}
    dataFlag::String
    indices::Vector{Int64}
    time::Float64
end

function read_logfile(filename)
    
    f = open(filename,"r")  # Open file in read mode
    nParse = 1    

    dataFlagNow = ""
    stepIndexNow = -1
    timeNow = -1.0
    D = Vector{Vector{Float64}}()
    I = Vector{Int64}()
    DD = Dict{Int64,stepData}() 
       
    # *Step  = 0
    # *Time  = 0
    # *Data  = ux;uy;uz
    # 1,0,0,0
    for l ∈ eachline(f) # Loop over each line in the file
        if contains(l,"*") # This means 
            sv = strip(split(l,'=')[2]) # Split the string using = and get string after             
            if contains(l,"*Step")
                if stepIndexNow>=0
                    DD[stepIndexNow] = stepData(D,dataFlagNow,I,timeNow)
                    D = Vector{Vector{Float64}}()
                    I = Vector{Int64}()
                end
                stepIndexNow=parse(Int64,sv)                            
            elseif contains(l,"*Time")
                timeNow=parse(Float64,sv)            
            elseif contains(l,"*Data")
                dataFlagNow=sv  
            end
        else
            ss = strip.(split(l,',')) # Split based on , as delimiter
            push!(I,parse(Int64,ss[1])) # Add index
            push!(D,parse.(Float64,ss[2:end])) # Add data
        end
    end
    DD[stepIndexNow]=stepData(D,dataFlagNow,I,timeNow) # Add last data

    close(f)

    return DD
end

function run_febio(filename_FEB,FEBIO_PATH="febio")
    run_filename = filename_FEB
    
    @static if Sys.islinux()
        runCommand = `nice "$FEBIO_PATH" "$run_filename"`    
    elseif Sys.isapple()
        runCommand = `nice "$FEBIO_PATH" "$run_filename"`
    else
        runCommand = `"$FEBIO_PATH" "$run_filename"` # runCommand = `start /min "$FEBIO_PATH" "$run_filename"`
    end

    
    run(runCommand)
end