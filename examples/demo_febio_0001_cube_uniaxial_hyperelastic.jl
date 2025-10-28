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
sampleSize = 10.0
pointSpacing = 2.0
strainApplied = 0.5 # Equivalent linear strain
loadingOption ="compression" # "tension" or "compression"

E_youngs = 1
ν = 0.4

###### 
# Creating a hexahedral mesh for a cube 
boxDim = sampleSize.*[1,1,1] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

# Create face sets to define node sets later 
Fb_top = Fb[CFb_type.==1]
Fb_bottom = Fb[CFb_type.==2]
Fb_s1 = Fb[CFb_type.==6]
Fb_s2 = Fb[CFb_type.==3]

# Defining displacement of the top surface in terms of x, y, and z components
if loadingOption=="tension"
    displacement_prescribed = strainApplied*sampleSize
elseif loadingOption=="compression"
    displacement_prescribed = -strainApplied*sampleSize
end

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
filename_force = "febioInputFile_01_REACTION_FORCE.txt"

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

material_node = aen(Material_node,"material"; id = 1, name="Material1", type="neo-Hookean")
    aen(material_node,"E",E_youngs)
    aen(material_node,"v",ν)

Mesh_node = aen(febio_spec_node,"Mesh")

# Nodes
Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for (i,v) in enumerate(V)        
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ v],','); id = i)     
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for (i,e) in enumerate(E)     
        aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = i)
    end
    
# Node sets
bcPrescribeList_z = "bcPrescribeList_z"
bcSupportList_x = "bcSupportList_x"
bcSupportList_y = "bcSupportList_y"
bcSupportList_z = "bcSupportList_z"
indNodes_bcPrescribeList_z = elements2indices(Fb_bottom)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indNodes_bcPrescribeList_z],','); name=bcPrescribeList_z)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_s1)],','); name=bcSupportList_x)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_s2)],','); name=bcSupportList_y)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_top)],','); name=bcSupportList_z)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_x", node_set=bcSupportList_x, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",0)
    aen(bc_node,"z_dof",0)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_y", node_set=bcSupportList_y, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",0)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_z", node_set=bcSupportList_z, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",0)
    aen(bc_node,"z_dof",1)

bc_node4 = aen(Boundary_node,"bc"; name="prescribed_disp_z", node_set=bcPrescribeList_z, type="prescribed displacement")
    aen(bc_node4,"dof","z")
    aen(bc_node4,"value",displacement_prescribed; lc=@sprintf("%i",1))
    aen(bc_node4,"relative",@sprintf("%i",0))

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
    aen(logfile_node,"element_data"; data="sz", delim=",", file=filename_stress)
    aen(logfile_node,"node_data"; data="Rx;Ry;Rz", delim=",", file=filename_force)


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
DD_force = read_logfile(joinpath(saveDir,filename_force))
numInc = length(DD_disp)
incRange = 0:1:numInc-1

F = element2faces(E)
indBoundaryFaces = boundaryfaceindices(F)
Fb = F[indBoundaryFaces]
# Fbs, Vs = separate_vertices(Fb,V)
Sz = [d[1] for d in DD_stress[numInc-1].data] 
Sz_F = repeat(Sz,inner=6)
Sz_Fb = Sz_F[indBoundaryFaces]
Sz_V = simplex2vertexdata(Fb, Sz_Fb,V)

# Create z-force curve data 

Fz_curve = zeros(numInc)
time_curve = zeros(numInc)
for i in 0:1:numInc-1    
    Fz = [d[3] for d in DD_force[i].data]
    Fz_curve[i+1] = sum(Fz[indNodes_bcPrescribeList_z])
    time_curve[i+1] = DD_force[i].time
end

#######
# Visualization

# Create time varying vectors
UT = fill(V,numInc) 
VT = fill(V,numInc)
SzT = fill(Sz_V,numInc)
UT_mag = fill(zeros(length(V)),numInc)
ut_mag_max = zeros(numInc)
Sz_max = zeros(numInc)
Sz_min = zeros(numInc)
@inbounds for i in 0:1:numInc-1    
    UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
    VT[i+1] += UT[i+1]
    UT_mag[i+1] = norm.(UT[i+1])
    ut_mag_max[i+1] = maximum(UT_mag[i+1]) 

    Sz = [d[1] for d in DD_stress[i].data]
    Sz_F = repeat(Sz,inner=6)
    Sz_Fb = Sz_F[indBoundaryFaces]
    Sz_V = simplex2vertexdata(Fb, Sz_Fb, V)
    SzT[i+1] = Sz_V
    Sz_min[i+1] = minimum(Sz)
    Sz_max[i+1] = maximum(Sz)
end

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

GLMakie.closeall()

fig = Figure(size=(1200, 800))
stepStart = incRange[end]
ax1 = AxisGeom(fig[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=2, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral),colorrange=(0.0, maximum(ut_mag_max)))
Colorbar(fig[1, 2], hp1.plots[1], label = "Displacement magnitude [mm]") 

ax2 = AxisGeom(fig[1, 3], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp2 = meshplot!(ax2, Fb, VT[end]; strokewidth=2, color=SzT[end], transparency=false, colormap = :viridis,colorrange=(minimum(Sz_min), maximum(Sz_max)))
Colorbar(fig[1, 4], hp2.plots[1], label = "σ [MPa]") 

ax3 = Axis(fig[1, 5], title = "Step: $stepStart", aspect = AxisAspect(1), xlabel="Time [s]", ylabel="Force [N]")
lines!(ax3, time_curve, Fz_curve, color=:red, linewidth=3)
hp3 = scatter!(ax3, Point{2,Float64}(time_curve[stepStart+1], Fz_curve[stepStart+1]), markersize=15, color=:red)
hSlider = Slider(fig[2, :], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider.value) do stepIndex 
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp1.color = UT_mag[stepIndex+1]

    hp2[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp2.color = SzT[stepIndex+1]
        
    hp3[1] = Point{2,Float64}(time_curve[stepIndex+1], Fz_curve[stepIndex+1])
    ax1.title = "Step: $stepIndex"
    ax2.title = "Step: $stepIndex"
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider,ax1)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")