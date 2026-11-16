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
    date_time_now = replace(string(now()), "T" => " at ")
    
    comment_node = XML.Comment("Created using FEBio.jl on "*date_time_now)
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
    step::Int64
end

function read_logfile(filename)    
    f = open(filename,"r")  # Open file in read mode
    
    dataFlagNow = ""
    stepIndexNow = -1
    timeNow = -1.0
    stepNow = 0
    D = Vector{Vector{Float64}}()
    I = Vector{Int64}()
    DD = Dict{Int64,stepData}() 
    
    # *Step  = 0
    # *Time  = 0
    # *Data  = ux;uy;uz
    # 1,0,0,0    
    for l ∈ eachline(f) # Loop over each line in the file
        if contains(l,"*") # Header line
            sv = strip(split(l,'=')[2]) # Split the string using = and get string after             
            if contains(l,"*Step")
                if stepIndexNow>=0                    
                    DD[stepNow] = stepData(D, dataFlagNow, I, timeNow, stepIndexNow)
                    D = Vector{Vector{Float64}}()
                    I = Vector{Int64}()
                    stepNow += 1
                end
                stepIndexNow=parse(Int64,sv)                            
            elseif contains(l,"*Time")
                timeNow=parse(Float64,sv)            
            elseif contains(l,"*Data")
                dataFlagNow=sv  
            end
        else # Data line
            ss = strip.(split(l,',')) # Split based on , as delimiter
            push!(I,parse(Int64,ss[1])) # Add index
            push!(D,parse.(Float64,ss[2:end])) # Add data
        end
    end
    DD[stepNow]=stepData(D, dataFlagNow, I, timeNow, stepIndexNow) # Add last data

    close(f)

    return DD
end

function showlog(logFileName)
    log_IO = open(logFileName)
    for line in eachline(log_IO)
        println(line)
    end
end

function _wait4logfile(process_id, logFileName; max_wait_time=10.0, sleep_time=0.01)
    log_file_flag = 0 
    t1 = time()
    while process_running(process_id)
        if isfile(logFileName)
            return 1            
        elseif time() - t1 > max_wait_time 
            return -1            
        end
        sleep(sleep_time)
    end
    return log_file_flag
end

function _kill_febio(process_id)
    try 
        kill(process_id)
    catch
        @warn "Killing process unsuccesful. Perhaps due to insufficient permissions. FEBio may still be running."
    end
end

function check_febio_log(process_id::Base.Process, logFileName; max_wait_time=20.0, sleep_time=0.01, quiet=false)    
    if isfile(logFileName) # If the log file exists 
        termination_type = :failed_process # Initiate as failed for now
        t1 = time()
        open(logFileName) do io
            if !quiet
                println("--------------- FEBio ---------------")
            end
            while true
                line = readline(io)
                if !quiet && contains(line, "number of iterations")
                    println(line)                    
                elseif !quiet && contains(line, "number of reformations")
                    println(line)                    
                elseif !quiet && contains(line, "------- converged at time :")
                    s = replace(line, "------- converged at time : " => "-> Converged for t=")
                    println(s)                    
                elseif !quiet && contains(line, "------- failed to converge at time :")
                    s = replace(line, "------- failed to converge at time : " => "-> Failed to converge for t=")
                    println(s)                    
                elseif !quiet && contains(line, "Total elapsed time")
                    s = replace(lstrip(line), "Total elapsed time .............. : " => "Elapsed time: ") 
                    println(s)                     
                elseif contains(line, "N O R M A L   T E R M I N A T I O N")
                    if !quiet
                        println(lstrip(line))
                        println("--------------- DONE  ---------------")
                        println("")
                    end
                    return termination_type = :normal                
                elseif contains(line, "E R R O R   T E R M I N A T I O N")
                    if !quiet
                        println(lstrip(line))
                        println("--------------- DONE  ---------------")
                        println("")
                    end
                    @warn """FEBio terminated with an error! Study the log file for more information, e.g. use showlog("$logFileName")."""
                    return termination_type = :error                
                end
                if eof(io) # Currently at end of file so wait to see if it will grow 
                    if process_running(process_id)               
                        if time() - t1 > max_wait_time
                            @warn "FEBio process time out. Log file has not changed for more than $max_wait_time seconds (max_wait_time). To trouble shoot, try increasing max_wait_time. Alternatively use run_mode=:internal and quiet=false and check terminal for output and potential error information."
                            _kill_febio(process_id)
                            return termination_type = :failed_process_timeout                        
                        end
                    else
                        @warn """FEBio process appears to have terminated prematurely. Study terminal output and the log file for more information, e.g. use showlog("$logFileName"). In addition, study the .feb file and importing it into FEBioStudio may help to diagnose issues."""
                        return termination_type = :failed_process                    
                    end
                    sleep(sleep_time)
                else
                    t1 = time() # Not end of file yet so reset timer       
                end
            end
        end
        return termination_type
    else
        @warn "No log file found. Perhaps the FEBio process terminated prematurely. Use run_mode=:internal and quiet=false and check terminal for output and potential error information. In addition, study the .feb file and importing it into FEBioStudio may help to diagnose issues."
        return termination_type = :failed_process
    end    
end

function run_febio(run_filename::String, logFileName::String; febio_path="febio4", run_mode=:external, max_wait_time=30.0, quiet=false, sleep_time=0.01)   
    # Set up run command for FEBio process
    @static if Sys.islinux()
        runCommand = `nice "$febio_path" "$run_filename"`    
    elseif Sys.isapple()
        runCommand = `nice "$febio_path" "$run_filename"`
    else
        runCommand = `"$febio_path" "$run_filename"` # runCommand = `start /min "$FEBIO_PATH" "$run_filename"`
    end 

    # Remove preexisting log files 
    if isfile(logFileName)
        try
            rm(logFileName)
        catch
            @warn "Unable to delete pre-existing FEBio log file: $logFileName, via rm function. Check if Julia has permissions to delete files."
        end
    end

    # Trigger FEBio
    if run_mode == :internal
        if quiet
            process_id = run(pipeline(runCommand, stdout=devnull); wait=true)
        else 
            process_id = run(pipeline(runCommand); wait=true)            
        end
        if isfile(logFileName)            
            # note quiet is alway true in check_febio_log below since either it is true or we've just shown all output in the terminal already
            run_flag = check_febio_log(process_id, logFileName; max_wait_time=max_wait_time, sleep_time=sleep_time, quiet=true)
            return run_flag
        else
            @warn "No log file created. Perhaps FEBio instantly failed without output. Check terminal for output, use run_mode=:internal and quiet=false to see potential error information."
            return :failed_process
        end
return run_flag    
    elseif run_mode == :external
        process_id = run(pipeline(runCommand); wait=false) # Trigger FEBio run and do not wait
        log_file_flag = _wait4logfile(process_id, logFileName; max_wait_time=max_wait_time, sleep_time=sleep_time)        
        if log_file_flag == 1 
            run_flag = check_febio_log(process_id, logFileName; max_wait_time=max_wait_time, sleep_time=sleep_time, quiet=quiet)
            return run_flag
        elseif log_file_flag == 0
            @warn "No log file created and FEBio process not running. Perhaps it instantly failed without output. Check terminal for output, use run_mode=:internal and quiet=false to see potential error information."
            return :failed_process
        elseif log_file_flag == -1        
            @warn "Log file not created within $max_wait_time seconds (max_wait_time). To trouble shoot, try increasing max_wait_time. Alternatively use run_mode=:internal and quiet=false and check terminal for output and potential error information."
            return :failed_logfile_timeout
        end
    end
    

end