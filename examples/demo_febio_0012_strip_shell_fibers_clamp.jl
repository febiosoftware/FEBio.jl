using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.LinearAlgebra
using FEBio
using FEBio.XML
using Printf 

######
GLMakie.closeall()

# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

###### 
# Control parameters 
pointSpacing = 2.0
strainApplied = 1.0 # Equivalent linear strain
strainClamp = 0.3 # Equivalent linear strain
loadingOption ="tension" # "tension" or "compression"

c = 1.0
m = 2.0
κp = c*100
d = 1e-9
shellThickness = 0.25 

ξ₁  = 10.0
α₁  = 0.0
β₁  = 2.0
λ0₁ = 1.0
θ₁  = 90.0
ϕ₁  = 45.0

ξ₂  = 10.0
α₂  = 0.0
β₂  = 2.0
λ0₂ = 1.0
θ₂  = θ₁
ϕ₂  = ϕ₁ + 90.0

testHeight = 45
testWidth = 35
sampleThickness = 6.0
clampHeight = 35.0

# Defining displacement of the top surface in terms of x, y, and z components
displacement_prescribed_z = strainApplied*testHeight
displacement_prescribed_x = strainClamp*sampleThickness

###### 
# Creating a hexahedral mesh for a cube 
boxDim1 = [sampleThickness, testWidth, testHeight] # Dimensionsions for the box in each direction
boxEl1 = ceil.(Int64, boxDim1./pointSpacing) # Number of elements to use in each direction 
E1, V1, F1, Fb1, Cb1 = hexbox(boxDim1,boxEl1)

boxDim2 = [sampleThickness, testWidth, clampHeight] # Dimensionsions for the box in each direction
boxEl2 = [boxEl1[1], boxEl1[2], ceil.(Int64, boxDim2[3]/pointSpacing)]# Number of elements to use in each direction 
E2, V2, F2, Fb2, Cb2 = hexbox(boxDim2,boxEl2)
E3 = deepcopy(E2)
Fb3 = deepcopy(Fb2)
V3 = deepcopy(V2)
V2 = [Point{3,Float64}(v[1], v[2], v[3]+testHeight/2.0+clampHeight/2.0) for v in V2]
V3 = [Point{3,Float64}(v[1], v[2], v[3]-testHeight/2.0-clampHeight/2.0) for v in V3]
F3 = deepcopy(F2) 

indexOffset_E2 = length(V1)
E2 = [Hex8{Int64}(e.+indexOffset_E2) for e in E2]
Fb2 = [QuadFace{Int64}(f.+indexOffset_E2) for f in Fb2]
F2 = [QuadFace{Int64}(f.+indexOffset_E2) for f in F2]

indexOffset_E3 = indexOffset_E2 + length(V2)
E3 = [Hex8{Int64}(e.+indexOffset_E3) for e in E3]
Fb3 = [QuadFace{Int64}(f.+indexOffset_E3) for f in Fb3]
F3 = [QuadFace{Int64}(f.+indexOffset_E3) for f in F3]

E = [E1; E2; E3]
V = [V1; V2; V3]
F = [F1; F2; F3]

E, V, indUnique, indMap = mergevertices(E, V; pointSpacing=pointSpacing)
indexmap!(Fb1, indMap)
indexmap!(Fb2, indMap)
indexmap!(Fb3, indMap)
indexmap!(F, indMap)

Fb = boundaryfaces(F)

# # Create face sets to define node sets later 
# Fb_bottom = Fb[CFb_type.==1]
# Fb_top = Fb[CFb_type.==2]

Fb_middle_front = Fb1[Cb1.==6]
Fb_top_front = Fb2[Cb2.==6]
Fb_top_back = Fb2[Cb2.==5]

Fb_bottom_front = Fb3[Cb2.==6]
Fb_bottom_back = Fb3[Cb2.==5]

F_shell = [Fb_top_front; Fb_middle_front; Fb_bottom_front]

indicesNodes_Ux_pre_pos = [elements2indices(Fb_top_front); elements2indices(Fb_bottom_front)]
indicesNodes_Ux_pre_neg = [elements2indices(Fb_top_back); elements2indices(Fb_bottom_back)] 

indicesNodes_Uz_pre = [elements2indices(Fb_top_front); elements2indices(Fb_top_back)] 
indicesNodes_Uz_fix = [elements2indices(Fb_bottom_front); elements2indices(Fb_bottom_back)]

indicesNodes_Uy_fix = [indicesNodes_Uz_pre; indicesNodes_Uz_fix]; 

# indicesNodes_bottom = elements2indices(Fb_bottom)

# # Visualisation
Fp, Vp = separate_vertices(Fb,V)
# Fp_bottom = Fp[CFb_type.==1]
# Fp_top = Fp[CFb_type.==2]
# Fp_s1 = Fp[CFb_type.==6]
# Fp_s2 = Fp[CFb_type.==3]

fig = Figure(size=(1200,1000))
ax1 = AxisGeom(fig[1, 1], title="Material regions")
hp1 = meshplot!(ax1, Fp, Vp; color=:white, strokewidth=1.0)
hp2 = meshplot!(ax1, F_shell, V; color=:lightgreen, strokewidth=1.0)

Legend(fig[1, 2],
    [hp1, hp2],
    ["Silicone hexahedral elements", "Fabric shell quadrilateral elements"])


ax2 = AxisGeom(fig[1, 3], title="Boundary conditions step 1")
meshplot!(ax2, Fp, Vp; color=(:white, 0.5), strokewidth=0.0, transparency=true)

hp3 = scatter!(ax2, V[indicesNodes_Ux_pre_pos], color=:red, markersize=10, depth_shift=-0.01f0)
hp4 = scatter!(ax2, V[indicesNodes_Ux_pre_neg], color=:black, markersize=10, depth_shift=-0.01f0)

Legend(fig[1, 4],
    [hp3, hp4],
    ["Front clamp nodes", "Back clamp nodes"])

ax3 = AxisGeom(fig[1, 5], title="Boundary conditions step 2")
meshplot!(ax3, Fp, Vp; color=(:white, 0.5), strokewidth=0.0, transparency=true)

hp5 = scatter!(ax3, V[indicesNodes_Uz_pre], color=:red, markersize=10, depth_shift=-0.01f0)
hp6 = scatter!(ax3, V[indicesNodes_Uz_fix], color=:black, markersize=10, depth_shift=-0.01f0)

Legend(fig[1, 6],
    [hp5, hp6],
    ["Top clamp nodes", "Bottom clamp nodes"])

screen = display(GLMakie.Screen(), fig)


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

Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Material1", type="Ogden unconstrained")
    aen(material_node,"c1",  c)
    aen(material_node,"m1",  m)
    aen(material_node,"c2",  c)
    aen(material_node,"m2", -m)
    aen(material_node,"cp", κp)
    aen(material_node,"density", d)

    material_node = aen(Material_node,"material"; id = 2, name="Material2", type="solid mixture")

    solid_node_01 = aen(material_node,"solid"; type="Ogden unconstrained")
        aen(solid_node_01,"c1",  c)
        aen(solid_node_01,"m1",  m)
        aen(solid_node_01,"c2",  c)
        aen(solid_node_01,"m2", -m)
        aen(solid_node_01,"cp", κp)
        aen(solid_node_01,"density", d)

    solid_node_02 = aen(material_node,"solid"; type="fiber-exp-pow")
        aen(solid_node_02,"ksi", ξ₁)
        aen(solid_node_02,"alpha", α₁)
        aen(solid_node_02,"beta", β₁)
        aen(solid_node_02,"lam0", λ0₁)
        mat_axis_node = aen(solid_node_02,"fiber"; type="angles")
            aen(mat_axis_node,"theta", θ₁)
            aen(mat_axis_node,"phi", ϕ₁)

    solid_node_03 = aen(material_node,"solid"; type="fiber-exp-pow")
        aen(solid_node_03,"ksi", ξ₂)
        aen(solid_node_03,"alpha", α₂)
        aen(solid_node_03,"beta", β₂)
        aen(solid_node_03,"lam0", λ0₂)
        mat_axis_node = aen(solid_node_03,"fiber"; type="angles")
            aen(mat_axis_node,"theta", θ₂)
            aen(mat_axis_node,"phi", ϕ₂)

            # mat_axis_node = aen(solid_node_02,"mat_axis"; type="vector")
            # aen(mat_axis_node,"a", join([@sprintf("%.16e",x) for x ∈ [1.0, 0.0, 0.0]]))
            # aen(mat_axis_node,"d", join([@sprintf("%.16e",x) for x ∈ [0.0, 1.0, 0.0]]))

    # <material id="1" type="solid mixture">
    #     <mat_axis type="local">0,0,0</mat_axis>
    #     <solid type="neo-Hookean">
    #         <E>1000.0</E>
    #         <v>0.45</v>
    #     </solid>
    #     <solid type="fiber-exp-pow">
    #         <ksi>5</ksi>
    #         <alpha>20</alpha>
    #         <beta>3</beta>
    #         <mat_axis type="angles">
                # <theta>0</theta>
                # <phi>90</phi>
    #         </mat_axis>
    #     </solid>
    # </material>

Mesh_node = aen(febio_spec_node,"Mesh")

# Nodes
Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for (i,v) in enumerate(V)        
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ v],','); id = i)     
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for (iElem,e) in enumerate(E)     
        aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = iElem)
    end

Elements_node = aen(Mesh_node,"Elements"; name="Part2", type="quad4")
    for (i,e) in enumerate(F_shell)
        iElem = i + length(E)     
        aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = iElem)
    end

# Node sets
# indicesNodes_Ux_pos = [elements2indices(Fb_top_front); elements2indices(Fb_bottom_front)]
# indicesNodes_Ux_neg = [elements2indices(Fb_top_back); elements2indices(Fb_bottom_back)] 
# indicesNodes_Uz_pos = [elements2indices(Fb_top_front); elements2indices(Fb_top_back)] 
# indicesNodes_Uz_neg = [elements2indices(Fb_bottom_front); elements2indices(Fb_bottom_back)]
# indicesNodes_Uy_fix = [indicesNodes_Uz_pos; indicesNodes_Uz_neg]; 

bcPrescribe_z_pre = "bcPrescribe_z_pre"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesNodes_Uz_pre],','); name=bcPrescribe_z_pre)

bcSupport_z_fix = "bcSupport_z_fix"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesNodes_Uz_fix],','); name=bcSupport_z_fix)

bcPrescribe_x_pre_pos = "bcPrescribe_x_pre_pos"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesNodes_Ux_pre_pos],','); name=bcPrescribe_x_pre_pos)

bcPrescribe_x_pre_neg = "bcPrescribe_x_pre_neg"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesNodes_Ux_pre_neg],','); name=bcPrescribe_x_pre_neg)

bcSupport_y_fix = "bcSupport_y_fix"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesNodes_Uy_fix],','); name=bcSupport_y_fix)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    SolidDomain_node = aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")
    ShellDomain_node = aen(MeshDomains_node,"ShellDomain"; mat = "Material2", name="Part2")
                            aen(ShellDomain_node,"shell_thickness", shellThickness)



Step_node = aen(febio_spec_node,"Step") # Define Step node: <Step>

# STEP 1
step_node_1 = aen(Step_node,"step"; id=1) # Define step node: <step>
control_node = aen(step_node_1,"Control") # Define Control node: <Control>
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

# STEP 2    
step_node_2 = aen(Step_node,"step"; id=2) # Define step node: <step>
control_node = aen(step_node_2,"Control") # Define Control node: <Control>
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
    

Boundary_node = aen(step_node_1, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_y", node_set=bcSupport_y_fix, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",0)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_z", node_set=bcPrescribe_z_pre, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",0)
    aen(bc_node,"z_dof",1)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_z", node_set=bcSupport_z_fix, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",0)
    aen(bc_node,"z_dof",1)

bc_node = aen(Boundary_node,"bc"; name="prescribed_displacement_x", node_set=bcPrescribe_x_pre_pos, type="prescribed displacement")
    aen(bc_node,"dof","x")
    aen(bc_node,"value", displacement_prescribed_x/2.0; lc=@sprintf("%i",1))
    aen(bc_node,"relative",@sprintf("%i",0))

bc_node = aen(Boundary_node,"bc"; name="prescribed_displacement_x", node_set=bcPrescribe_x_pre_neg, type="prescribed displacement")
    aen(bc_node,"dof","x")
    aen(bc_node,"value", -displacement_prescribed_x/2.0; lc=@sprintf("%i",1))
    aen(bc_node,"relative",@sprintf("%i",0))

Boundary_node = aen(step_node_2, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_y", node_set=bcSupport_y_fix, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",0)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_z", node_set=bcSupport_z_fix, type="zero displacement")
    aen(bc_node,"x_dof",0)
    aen(bc_node,"y_dof",0)
    aen(bc_node,"z_dof",1)

bc_node = aen(Boundary_node,"bc"; name="prescribed_displacement_x", node_set=bcPrescribe_x_pre_pos, type="prescribed displacement")
    aen(bc_node,"dof","x")
    aen(bc_node,"value", displacement_prescribed_x/2.0; lc=@sprintf("%i",1))
    aen(bc_node,"relative",@sprintf("%i",0))

bc_node = aen(Boundary_node,"bc"; name="prescribed_displacement_x", node_set=bcPrescribe_x_pre_neg, type="prescribed displacement")
    aen(bc_node,"dof","x")
    aen(bc_node,"value", -displacement_prescribed_x/2.0; lc=@sprintf("%i",1))
    aen(bc_node,"relative",@sprintf("%i",0))

bc_node = aen(Boundary_node,"bc"; name="prescribed_displacement_z", node_set=bcPrescribe_z_pre, type="prescribed displacement")
    aen(bc_node,"dof","z")
    aen(bc_node,"value", displacement_prescribed_z; lc=@sprintf("%i",2))
    aen(bc_node,"relative",@sprintf("%i",1))

LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0.0, 0.0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1.0, 1.0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",2.0, 1.0))

load_controller_node = aen(LoadData_node,"load_controller"; id=2, name="LC_2", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0.0, 0.0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1.0, 0.0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",2.0, 1.0))

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
time_vec = zeros(numInc)
@inbounds for i in 0:1:numInc-1    
    UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
    VT[i+1] += UT[i+1]
    UT_mag[i+1] = norm.(UT[i+1])
    ut_mag_max[i+1] = maximum(UT_mag[i+1]) 
    time_vec[i+1] = DD_disp[i].time
end

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

#######
# Visualization
fig = Figure(size=(800,800))
stepStart = incRange[end]
timeStart = time_vec[stepStart+1]
ax = AxisGeom(fig[1, 1], title = "Step: $stepStart, Time: $timeStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp = meshplot!(ax, Fb, VT[end]; strokewidth=2, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral),colorrange=(0,maximum(ut_mag_max)))
Colorbar(fig[1, 2],hp.plots[1],label = "Displacement magnitude [mm]") 

hSlider = Slider(fig[2, 1], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider.value) do stepIndex     
    hp[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp.color = UT_mag[stepIndex+1]
    timeNow = time_vec[stepIndex+1]
    ax.title = "Step: $stepIndex, Time: $timeNow"
end

slidercontrol(hSlider,ax)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")