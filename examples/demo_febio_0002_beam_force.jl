using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.LinearAlgebra
using FEBio
using FEBio.XML
using Printf 

######
# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

###### 
# Control parameters 
pointSpacing = 3.0
appliedForce = (0.0, 0.0, -2e-3) 

# Material parameters
c1 = 1e-3 #Shear-modulus-like parameter
m1 = 8  #Material parameter setting degree of non-linearity
k_factor = 1e2  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus

###### 
# Creating a hexahedral mesh for a cube 
boxDim = [10.0,40.0,10.0] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

# Create face sets to define node sets later 
Fb_top = Fb[CFb_type.==1]
Fb_bottom = Fb[CFb_type.==2]
Fb_s1 = Fb[CFb_type.==6]
Fb_s2 = Fb[CFb_type.==3]

indNodesBack = elements2indices(Fb[CFb_type.==4])
indNodesFront = elements2indices(Fb[CFb_type.==3])

######
# Define file names
saveDir = joinpath(febiojl_dir(),"assets","temp") # Main directory to save FEBio input and output files
if !isdir(saveDir)
    mkdir(saveDir)      
end

filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file
filename_xplt = joinpath(saveDir,"febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
filename_log = joinpath(saveDir,"febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
filename_stress = "febioInputFile_01_STRESS.txt"

######
# Define febio input file XML
doc,febio_spec_node = feb_doc_initialize()

aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>

control_node = aen(febio_spec_node,"Control") # Define Control node: <Control>
    aen(control_node,"analysis","STATIC")               
    aen(control_node,"time_steps",10)
    aen(control_node,"step_size",0.1)
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)

time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",5)
    aen(time_stepper_node,"opt_iter",10)
    aen(time_stepper_node,"dtmin",1e-3)
    aen(time_stepper_node,"dtmax",0.1)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",1)
    aen(solver_node,"equation_scheme",1)
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    aen(solver_node,"max_refs",70)
    aen(solver_node,"check_zero_diagonal",0)
    aen(solver_node,"zero_diagonal_tol",0)
    aen(solver_node,"force_partition",0)
    aen(solver_node,"reform_each_time_step",1)
    aen(solver_node,"reform_augment",0)
    aen(solver_node,"diverge_reform",1)
    aen(solver_node,"min_residual",1e-20)
    aen(solver_node,"max_residual",0)
    aen(solver_node,"dtol",1e-3)
    aen(solver_node,"etol",1e-2)
    aen(solver_node,"rtol",0)
    aen(solver_node,"rhoi",0)
    aen(solver_node,"alpha",1)
    aen(solver_node,"beta",2.5e-01)
    aen(solver_node,"gamma",5e-01)
    aen(solver_node,"logSolve",0)
    aen(solver_node,"arc_length",0)
    aen(solver_node,"arc_length_scale",0)
qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
    aen(qn_method_node,"max_ups",0)
    aen(qn_method_node,"max_buffer_size",0)
    aen(qn_method_node,"cycle_buffer",0)
    aen(qn_method_node,"cmax",0)

Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Material1", type="Ogden")
    aen(material_node,"c1",c1)
    aen(material_node,"m1",m1)
    aen(material_node,"c2",c1)
    aen(material_node,"m2",-m1)
    aen(material_node,"k",k)

Mesh_node = aen(febio_spec_node,"Mesh")

# Nodes
Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for q ∈ eachindex(V)
        # aen(Nodes_node,"node",@sprintf("%.2f, %.2f, %.2f",V[q][1],V[q][2],V[q][3]); id = q)
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ V[q]],','); id = q)     
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for q ∈ eachindex(E)
        # aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q][1],E[q][2],E[q][3],E[q][4],E[q][5],E[q][6],E[q][7],E[q][8]); id = q)
        aen(Elements_node,"elem",join(map(string, E[q]), ','); id = q)
    end
    
# Node sets
bcSupportList = "bcSupportList"
bcPrescribe = "bcPrescribe"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indNodesBack],','); name=bcSupportList)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indNodesFront],','); name=bcPrescribe)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xyz", node_set=bcSupportList, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

Loads_node = aen(febio_spec_node, "Loads")
load_node = aen(Loads_node,"nodal_load"; name="PrescribedForce", node_set=bcPrescribe, type="nodal_force")
    aen(load_node,"value",join([@sprintf("%.16e",x) for x ∈ appliedForce./length(bcPrescribe)],','); lc=@sprintf("%i",1))

LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1,1))

Output_node = aen(febio_spec_node,"Output")

plotfile_node = aen(Output_node,"plotfile"; type="febio")
    aen(plotfile_node,"var"; type="displacement")
    aen(plotfile_node,"var"; type="stress")
    aen(plotfile_node,"var"; type="relative volume")
    aen(plotfile_node,"var"; type="reaction forces")
    aen(plotfile_node,"var"; type="contact pressure")
    aen(plotfile_node,"compression",@sprintf("%i",0))

logfile_node = aen(Output_node,"logfile"; file=filename_log)
    aen(logfile_node,"node_data"; data="ux;uy;uz", delim=",", file=filename_disp)
    aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
# <logfile file="tempModel.txt">
#   <node_data data="ux;uy;uz" delim="," file="tempModel_disp_out.txt">1, 2, 3, 4, 5, 6, 7, 8, 

#######
# Write FEB file
XML.write(filename_FEB, doc)

#######
# Run FEBio
run_febio(filename_FEB,FEBIO_EXEC)

#######
# Import results
DD_disp = read_logfile(joinpath(saveDir,filename_disp))
DD_stress = read_logfile(joinpath(saveDir,filename_stress))
numInc = length(DD_disp)
incRange = 0:1:numInc-1

# Create time varying vectors
UT = fill(V,numInc) 
VT = fill(V,numInc)
UT_mag = fill(zeros(length(V)),numInc)
ut_mag_max = zeros(numInc)
@inbounds for i in 0:1:numInc-1    
    UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
    VT[i+1] += UT[i+1]
    UT_mag[i+1] = norm.(UT[i+1])
    ut_mag_max[i+1] = maximum(UT_mag[i+1]) 
end

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

titleString = lift(hSlider.value) do stepIndex
    "Step: "*string(stepIndex)
end
#######
# Visualization
fig = Figure(size=(800,800))

ax = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString,
limits=(min_p[1],max_p[1], min_p[2],max_p[2], min_p[3],max_p[3]))
hp = poly!(GeometryBasics.Mesh(VT[end],Fb), strokewidth=2,color=UT_mag[end], transparency=false, overdraw=false,
colormap = Reverse(:Spectral),colorrange=(0,maximum(ut_mag_max)))
Colorbar(fig[1, 2],hp.plots[1],label = "Displacement magnitude [mm]") 


hSlider = Slider(fig[2, 1], range = incRange, startvalue = incRange[end],linewidth=30)
on(hSlider.value) do stepIndex 
    hp[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp.color = UT_mag[stepIndex+1]
    ax.title = "Step: "*string(stepIndex)
end

slidercontrol(hSlider,ax)

fig