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

strainApplied = 0.6 # Equivalent linear strain
loadingOption = "tension" # "tension" or "compression"
elementType = "hex20" # "hex8"
if elementType == "hex8"
    pointSpacing = 3.0/2.0
elseif elementType == "hex20"
    pointSpacing = 3.0
end

c = 1.0
m = 2.0
κp = c*100
d = 1e-9
shellThickness = 0.5

ξ₁ = 50.0
α₁ = 0.0
β₁ = 2.0
λ0₁ = 1.0
θ₁ = 90.0
ϕ₁ = 45.0

ξ₂ = 50.0
α₂ = 0.0
β₂ = 2.0
λ0₂ = 1.0
θ₂ = θ₁
ϕ₂ = ϕ₁ + 90.0

# FEA control settings
numTimeSteps = 20 # Number of time steps desired
max_refs = 50 # Max reforms
max_ups = 0 # Set to zero to use full-Newton iterations
opt_iter = 10 # Optimum number of iterations
max_retries = 5 # Maximum number of retires
dtmin = (1.0/numTimeSteps)/100.0 # Minimum time step size
dtmax = 1.0/numTimeSteps  # Maximum time step size
symmetric_stiffness = 1
min_residual = 1e-30

###### 
# Creating a hexahedral mesh for a cube 
boxDim = [5.0, 40.0, 40.0] # Dimensionsions for the box in each direction
boxEl1 = ceil(Int64, boxDim[1]/pointSpacing)
boxEl2 = ceil(Int64, boxDim[2]/pointSpacing)
boxEl3 = ceil(Int64, (boxDim[3]*(1.0+strainApplied))/pointSpacing)
if !iseven(boxEl3)
    boxEl3+=1
end

boxEl = [boxEl1, boxEl2, boxEl3]

E, V, F, Fb, CFb_type = hexbox(boxDim, boxEl)

if elementType == "hex20"
    E, V = hex8_hex20(E,V)
    F = element2faces(E)
    indBoundary = boundaryfaceindices(F)
    Fb = F[indBoundary]
    faceType = "quad8"
else 
    faceType = "quad4"
end

# Create face sets to define node sets later 
Fb_bottom = Fb[CFb_type .== 1]
Fb_top = Fb[CFb_type .== 2]
Fb_s1 = Fb[CFb_type .== 6]
Fb_s2 = Fb[CFb_type .== 3]

if elementType == "hex20"
    E2 = ngon8_quad8(Fb_s1)
else
    E2 = Fb_s1 
end
# Visualisation

cmap_cat = Makie.Categorical(:Spectral) 
Fbs, Vs = separate_vertices(Fb, V)
Cb_Vs = simplex2vertexdata(Fbs, CFb_type, Vs)

fig = Figure(size=(1200, 1000))

ax1 = AxisGeom(fig[1, 1], title="Boundary labels")
hp1 = meshplot!(ax1, Fbs, Vs; color=Cb_Vs, strokewidth=1.0, colormap=cmap_cat)
Colorbar(fig[1, 2], hp1)

ax2 = AxisGeom(fig[1, 3], title="Boundary conditions and shell layer")

hp2 = meshplot!(ax2, Fb, V; color=(:white, 0.25), strokewidth=0.0, transparency=true)
hp3 = meshplot!(ax2, Fb_s1, V; color=:green, strokewidth=1.0)
hp4 = meshplot!(ax2, Fb_top, V; color=:red, strokewidth=1.0)
hp5 = meshplot!(ax2, Fb_bottom, V; color=:blue, strokewidth=1.0)

# scatter!(ax1, V, color=:black, markersize=10, depth_shift=-0.01f0)
screen = display(GLMakie.Screen(), fig)


# Defining displacement of the top surface in terms of x, y, and z components
if loadingOption=="tension"
    displacement_prescribed = strainApplied*boxDim[3]
elseif loadingOption=="compression"
    displacement_prescribed = -strainApplied*boxDim[3]
end

######
# Define file names
saveDir = joinpath(febiojl_dir(), "assets", "temp") # Main directory to save FEBio input and output files
if !isdir(saveDir)
    mkdir(saveDir)
end

filename_FEB = joinpath(saveDir, "febioInputFile_01.feb")   # The main FEBio input file
filename_xplt = joinpath(saveDir, "febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
filename_log = joinpath(saveDir, "febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
filename_stress = "febioInputFile_01_STRESS.txt"
######
# Define febio input file XML
doc, febio_spec_node = feb_doc_initialize()

aen(febio_spec_node, "Module"; type="solid") # Define Module node: <Module type="solid"/>

control_node = aen(febio_spec_node, "Control") # Define Control node: <Control>
aen(control_node, "analysis", "STATIC")
aen(control_node, "time_steps", 10)
aen(control_node, "step_size", 0.1)
aen(control_node, "plot_zero_state", 1)
aen(control_node, "plot_range", @sprintf("%.2f, %.2f", 0, -1))
aen(control_node, "plot_level", "PLOT_MAJOR_ITRS")
aen(control_node, "plot_stride", 1)
aen(control_node, "output_level", "OUTPUT_MAJOR_ITRS")
aen(control_node, "adaptor_re_solve", 1)

time_stepper_node = aen(control_node, "time_stepper"; type="default")
aen(time_stepper_node, "max_retries", 5)
aen(time_stepper_node, "opt_iter", 10)
aen(time_stepper_node, "dtmin", 1e-3)
aen(time_stepper_node, "dtmax", 0.1)
aen(time_stepper_node, "aggressiveness", 0)
aen(time_stepper_node, "cutback", 5e-1)
aen(time_stepper_node, "dtforce", 0)

solver_node = aen(control_node, "solver"; type="solid")
aen(solver_node, "symmetric_stiffness", symmetric_stiffness)
aen(solver_node, "equation_scheme", 1)
aen(solver_node, "equation_order", "default")
aen(solver_node, "optimize_bw", 0)
aen(solver_node, "lstol", 9e-1)
aen(solver_node, "lsmin", 1e-2)
aen(solver_node, "lsiter", 5)
aen(solver_node, "max_refs", max_refs)
aen(solver_node, "check_zero_diagonal", 0)
aen(solver_node, "zero_diagonal_tol", 0)
aen(solver_node, "force_partition", 0)
aen(solver_node, "reform_each_time_step", 1)
aen(solver_node, "reform_augment", 0)
aen(solver_node, "diverge_reform", 1)
aen(solver_node, "min_residual", min_residual)
aen(solver_node, "max_residual", 0)
aen(solver_node, "dtol", 1e-3)
aen(solver_node, "etol", 1e-2)
aen(solver_node, "rtol", 0)
aen(solver_node, "rhoi", 0)
aen(solver_node, "alpha", 1)
aen(solver_node, "beta", 2.5e-01)
aen(solver_node, "gamma", 5e-01)
aen(solver_node, "logSolve", 0)
aen(solver_node, "arc_length", 0)
aen(solver_node, "arc_length_scale", 0)
qn_method_node = aen(solver_node, "qn_method"; type="BFGS")
aen(qn_method_node, "max_ups", max_ups)
aen(qn_method_node, "max_buffer_size", 0)
aen(qn_method_node, "cycle_buffer", 0)
aen(qn_method_node, "cmax", 0)

Globals_node = aen(febio_spec_node, "Globals")

Constants_node = aen(Globals_node, "Constants")
aen(Constants_node, "R", 8.3140000e-06)
aen(Constants_node, "T", 298)
aen(Constants_node, "F", 9.6485000e-05)

Material_node = aen(febio_spec_node, "Material")

material_node = aen(Material_node, "material"; id="1", name="Material1", type="Ogden unconstrained")
aen(material_node, "c1", c)
aen(material_node, "m1", m)
aen(material_node, "c2", c)
aen(material_node, "m2", -m)
aen(material_node, "cp", κp)
aen(material_node, "density", d)

material_node = aen(Material_node, "material"; id="2", name="Material2", type="solid mixture")

solid_node_01 = aen(material_node, "solid"; type="Ogden unconstrained")
aen(solid_node_01, "c1", c)
aen(solid_node_01, "m1", m)
aen(solid_node_01, "c2", c)
aen(solid_node_01, "m2", -m)
aen(solid_node_01, "cp", κp)
aen(solid_node_01, "density", d)

solid_node_02 = aen(material_node, "solid"; type="fiber-exp-pow")
aen(solid_node_02, "ksi", ξ₁)
aen(solid_node_02, "alpha", α₁)
aen(solid_node_02, "beta", β₁)
aen(solid_node_02, "lam0", λ0₁)
mat_axis_node = aen(solid_node_02, "fiber"; type="angles")
aen(mat_axis_node, "theta", θ₁)
aen(mat_axis_node, "phi", ϕ₁)

solid_node_03 = aen(material_node, "solid"; type="fiber-exp-pow")
aen(solid_node_03, "ksi", ξ₂)
aen(solid_node_03, "alpha", α₂)
aen(solid_node_03, "beta", β₂)
aen(solid_node_03, "lam0", λ0₂)
mat_axis_node = aen(solid_node_03, "fiber"; type="angles")
aen(mat_axis_node, "theta", θ₂)
aen(mat_axis_node, "phi", ϕ₂)

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

Mesh_node = aen(febio_spec_node, "Mesh")

# Nodes
Nodes_node = aen(Mesh_node, "Nodes"; name="nodeSet_all")
for (i, v) in enumerate(V)
    aen(Nodes_node, "node", join([@sprintf("%.16e", x) for x ∈ v], ','); id=@sprintf("%i", i))
end

# Elements
Elements_node = aen(Mesh_node, "Elements"; name="Part1", type=elementType)
for (iElem, e) in enumerate(E)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", iElem))
end

Elements_node = aen(Mesh_node, "Elements"; name="Part2", type=faceType)
for (i, e) in enumerate(E2)
    iElem = i + length(E)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", iElem))
end

# Node sets
bcPrescribeList_z = "bcPrescribeList_z"
bcSupportList_z = "bcSupportList_z"
aen(Mesh_node, "NodeSet", join([@sprintf("%i", x) for x ∈ elements2indices(Fb_top)], ','); name=bcPrescribeList_z)
aen(Mesh_node, "NodeSet", join([@sprintf("%i", x) for x ∈ elements2indices(Fb_bottom)], ','); name=bcSupportList_z)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
SolidDomain_node = aen(MeshDomains_node, "SolidDomain"; mat="Material1", name="Part1")
ShellDomain_node = aen(MeshDomains_node, "ShellDomain"; mat="Material2", name="Part2")
aen(ShellDomain_node, "shell_thickness", shellThickness)

Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node, "bc"; name="zero_displacement_z_bottom", node_set=bcSupportList_z, type="zero displacement")
aen(bc_node, "x_dof", 1)
aen(bc_node, "y_dof", 1)
aen(bc_node, "z_dof", 1)

bc_node = aen(Boundary_node, "bc"; name="zero_displacement_xy_top", node_set=bcPrescribeList_z, type="zero displacement")
aen(bc_node, "x_dof", 1)
aen(bc_node, "y_dof", 1)
aen(bc_node, "z_dof", 0)

bc_node4 = aen(Boundary_node, "bc"; name="prescribed_disp_z", node_set=bcPrescribeList_z, type="prescribed displacement")
aen(bc_node4, "dof", "z")
aen(bc_node4, "value", displacement_prescribed; lc=@sprintf("%i", 1))
aen(bc_node4, "relative", @sprintf("%i", 0))

LoadData_node = aen(febio_spec_node, "LoadData")

load_controller_node = aen(LoadData_node, "load_controller"; id="1", name="LC_1", type="loadcurve")
aen(load_controller_node, "interpolate", "LINEAR")

points_node = aen(load_controller_node, "points")
aen(points_node, "pt", @sprintf("%.2f, %.2f", 0.0, 0.0))
aen(points_node, "pt", @sprintf("%.2f, %.2f", 1.0, 1.0))

Output_node = aen(febio_spec_node, "Output")

plotfile_node = aen(Output_node, "plotfile"; type="febio")
aen(plotfile_node, "var"; type="displacement")
aen(plotfile_node, "var"; type="stress")
aen(plotfile_node, "var"; type="relative volume")
aen(plotfile_node, "var"; type="reaction forces")
aen(plotfile_node, "var"; type="contact pressure")
aen(plotfile_node, "compression", @sprintf("%i", 0))

logfile_node = aen(Output_node, "logfile"; file=filename_log)
aen(logfile_node, "node_data"; data="ux;uy;uz", delim=",", file=filename_disp)
aen(logfile_node, "element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
# <logfile file="tempModel.txt">
#   <node_data data="ux;uy;uz" delim="," file="tempModel_disp_out.txt">1, 2, 3, 4, 5, 6, 7, 8, 

#######
# Write FEB file
XML.write(filename_FEB, doc)

#######
# Run FEBio
run_febio(filename_FEB, FEBIO_EXEC)

#######
# Import results
DD_disp = read_logfile(joinpath(saveDir, filename_disp))
DD_stress = read_logfile(joinpath(saveDir, filename_stress))
numInc = length(DD_disp)
incRange = 0:1:(numInc-1)

# Create time varying vectors
UT = fill(V, numInc)
VT = fill(V, numInc)
UT_mag = fill(zeros(length(V)), numInc)
ut_mag_max = zeros(numInc)
@inbounds for i in 0:1:(numInc-1)
    UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
    VT[i+1] += UT[i+1]
    UT_mag[i+1] = norm.(UT[i+1])
    ut_mag_max[i+1] = maximum(UT_mag[i+1])
end

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

#######
# Visualization
fig = Figure(size=(800, 800))
stepStart = incRange[end]
ax = AxisGeom(fig[1, 1], title="Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp = meshplot!(ax, Fb, VT[end]; strokewidth=2, color=UT_mag[end], transparency=false, colormap=Reverse(:Spectral), colorrange=(0, maximum(ut_mag_max)))
Colorbar(fig[1, 2], hp.plots[1], label="Displacement magnitude [mm]")

hSlider = Slider(fig[2, 1], range=incRange, startvalue=stepStart, linewidth=30)
on(hSlider.value) do stepIndex
    hp[1] = GeometryBasics.Mesh(VT[stepIndex+1], Fb)
    hp.color = UT_mag[stepIndex+1]
    ax.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")