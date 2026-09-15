using Comodo
using Comodo.Rotations
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.LinearAlgebra
using Comodo.Statistics
using FEBio
using FEBio.XML
using Printf
using FileIO

GLMakie.closeall()

######
# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

## Control parameters 

# Material parameters
E_youngs = 1
ν = 0.45

# Contact parameters
contactPenalty = 100.0;
laugon = 0;
minaug = 1;
maxaug = 10;
fric_coeff = 0.01;

# FEA control settings
numTimeSteps = 10 # Number of time steps desired
max_refs = 50 # Max reforms
max_ups = 0 # Set to zero to use full-Newton iterations
opt_iter = 20 # Optimum number of iterations
max_retries = 5 # Maximum number of retires
dtmin = (1.0/numTimeSteps)/100.0 # Minimum time step size
dtmax = 1.0/numTimeSteps  # Maximum time step size
symmetric_stiffness = 0
min_residual = 1e-20

cylDisplacement_Z = -20.0

fileName_stl = joinpath(febiojl_dir(), "assets", "stl", "femur_iso.stl")

## 
M = load(fileName_stl)
F1 = tofaces(faces(M))
V1 = [Point{3,Float64}(v) for v in topoints(coordinates(M))]
F1, V1, _, _ = mergevertices(F1, V1)
C1 = ones(Int, length(F1))

p = surface_centroid(F1, V1)
SVD = surface_svd(F1, V1)
V1 = [Point{3,Float64}(SVD.Vt*1000.0*(v-p)) for v in V1]

pointSpacing = pointspacingmean(F1, V1) # Point spacing of current mesh 

initialContactSpacing = pointSpacing/10.0

vol1 = (pointSpacing^3.0) / (6.0*sqrt(2.0))

V_regions = [faceinteriorpoint(F1, V1, 1)]
V_holes = Vector{Point{3,Float64}}()
region_vol = [vol1]

stringOpt = "paAqYQ"
E, V, CE, Fb, Cb = tetgenmesh(F1, V1; facetmarkerlist=C1, V_regions=V_regions, V_holes=V_holes, region_vol=region_vol, stringOpt=stringOpt, element_type=Tet4{Int})
F = element2faces(E)

## Create an position cylinders 

r = 12.0 # Cylinder radius
h = 120.0 # Cylinder height
nr = ceil(Int, (2π*r) ./ pointSpacing) # Derived number of radial points to approximate point spacing
nh = 1 + ceil(Int, h ./ pointSpacing) # Derived number of height direction points to approximate point spacing
F2, V2 = cylinder(r, h, nr, nh; direction=:both, face_type=:tri_even, face_orientation=:outward)

# Reorient cylinder 
Q = RotXYZ(0.5*π, 0.0, 0.0)
V2 = [GeometryBasics.Point{3,Float64}(Q*v) for v ∈ V2]

# Copy to create others 
V3 = deepcopy(V2)
pShift = Point{3,Float64}(-120.0, 0.0, 0.0)
V3 .+= pShift

V4 = deepcopy(V2)
pShift = Point{3,Float64}(120.0, 0.0, 0.0)
V4 .+= pShift

function getLocalMinMax_Z(Vb, Vc)
    Xc = [v[1] for v in Vc] # X coordinates 
    Xc_min = minimum(Xc)
    Xc_max = maximum(Xc)
    zMax = -Inf
    zMin = Inf
    for (i, v) in enumerate(Vb)
        if v[1]>=Xc_min && v[1]<=Xc_max
            z = v[3]
            zMax=max(zMax, z)
            zMin=min(zMin, z)
        end
    end
    return zMax, zMin
end

function get_z_shift(Fb, Vb, Vc, shiftDir)
    t = Inf
    for v in Vc
        _, _, T, _, _ = ray_triangle_intersect(Fb, Vb, v, shiftDir; rayType=:ray, triSide=0)
        if !isempty(T)
            t = min(t, minimum(T))
        end
    end
    return t
end

zMax, zMin = getLocalMinMax_Z(V1, V2)
pShift = Point{3,Float64}(0.0, 0.0, r+zMax)
V2 .+= pShift
t = get_z_shift(F1, V1, V2, Vec{3,Float64}(0.0, 0.0, -1.0))
pShift = Point{3,Float64}(0.0, 0.0, -t)
V2 .+= pShift

zMax, zMin = getLocalMinMax_Z(V1, V3)
pShift = Point{3,Float64}(0.0, 0.0, -r+zMin)
V3 .+= pShift
t = get_z_shift(F1, V1, V3, Vec{3,Float64}(0.0, 0.0, 1.0))
pShift = Point{3,Float64}(0.0, 0.0, t)
V3 .+= pShift

zMax, zMin = getLocalMinMax_Z(V1, V4)
pShift = Point{3,Float64}(0.0, 0.0, -r+zMin)
V4 .+= pShift
t = get_z_shift(F1, V1, V4, Vec{3,Float64}(0.0, 0.0, 1.0))
pShift = Point{3,Float64}(0.0, 0.0, t)
V4 .+= pShift

## Add cylinder models

# Shift sphere indices as they are appended after solid mesh
F2 = [eltype(F2)(f .+ length(V)) for f in F2]
F3 = deepcopy(F2)
F3 = [eltype(F3)(f .+ length(V2)) for f in F3]
F4 = deepcopy(F3)
F4 = [eltype(F4)(f .+ length(V3)) for f in F4]

# Append cylinder nodes
append!(V, V2)
append!(V, V3)
append!(V, V4)


## Visualization
cmap = cgrad(:Spectral, 5, categorical=true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE, inner=4)

Fbs, Vbs = separate_vertices(Fb, V)
Cbs_V = simplex2vertexdata(Fbs, Cb)

Fs, Vs = separate_vertices(F, V)
CE_Vs = simplex2vertexdata(Fs, CE_F)
M = GeometryBasics.Mesh(Vs, Fs)

strokewidth = 0.5

fig1 = Figure(size=(800, 800))

ax1 = AxisGeom(fig1[1, 1][1, 1], title="Boundary surfaces")
hp1_fig1 = meshplot!(ax1, Fbs, Vbs, color=(:white, 0.5), transparency=true, colorrange=(1, 3), colormap=cmap, strokewidth=0.0)

scatter!(ax1, V_regions, color=:black, markersize=25)
scatter!(ax1, V_holes, color=:black, markersize=25)

ax2 = AxisGeom(fig1[1, 1][1, 2], title="Cut mesh")
meshplot!(ax2, Fbs, Vbs, color=(:white, 0.15), transparency=true, colorrange=(1, 3), colormap=cmap, strokewidth=0.0)

hp2_fig1 = meshplot!(ax2, Fs, Vs, color=CE_Vs, strokewidth=strokewidth, strokecolor=:black, colorrange=(1, 2), colormap=cmap)
# hp3 = scatter!(ax2, V, color=:black, markersize=10)

VE = simplexcenter(E, V)
XE = [v[1] for v in VE]
X = [v[1] for v in V]
xMax = maximum(X)
xMin = minimum(X)
numSlicerSteps = 2*ceil(Int, (xMax-xMin)/pointSpacing)

stepRange = range(xMin, xMax, numSlicerSteps)
hSlider1 = Slider(fig1[2, 1], range=stepRange, startvalue=mean(stepRange), linewidth=30)

on(hSlider1.value) do x
    B = XE .>= x
    indShow = findall(B)
    if isempty(indShow)
        hp2_fig1.visible=false
    else
        hp2_fig1.visible=true
        Fs = element2faces(E[indShow])
        Cs = repeat(CE[indShow], inner=4)

        indB = boundaryfaceindices(Fs)
        Fs = Fs[indB]
        Cs = Cs[indB]
        Fs, Vs = separate_vertices(Fs, V)
        CE_Vs = simplex2vertexdata(Fs, Cs)
        Ms = GeometryBasics.Mesh(Vs, Fs)
        hp2_fig1[1] = Ms
        hp2_fig1.color = CE_Vs
        # hp3[1] = V[elements2indices(E[indShow])]
    end
end
# hSlider.selected_index[]+=1
slidercontrol(hSlider1, ax2)

screen = display(GLMakie.Screen(), fig1)

##

fig = Figure(size=(800, 800))

ax1 = AxisGeom(fig[1, 1], title="Boundary surfaces")
hp1 = meshplot!(ax1, Fb, V, color=:white, strokewidth=1.0)
hp2 = meshplot!(ax1, F2, V, color=:red, strokewidth=1.0)
hp3 = meshplot!(ax1, F3, V, color=:green, strokewidth=1.0)
hp4 = meshplot!(ax1, F4, V, color=:blue, strokewidth=1.0)

screen = display(GLMakie.Screen(), fig)

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
aen(control_node, "time_steps", numTimeSteps)
aen(control_node, "step_size", 1.0/numTimeSteps)
aen(control_node, "plot_zero_state", 1)
aen(control_node, "plot_range", @sprintf("%.2f, %.2f", 0, -1))
aen(control_node, "plot_level", "PLOT_MAJOR_ITRS")
aen(control_node, "plot_stride", 1)
aen(control_node, "output_level", "OUTPUT_MAJOR_ITRS")
aen(control_node, "adaptor_re_solve", 1)

time_stepper_node = aen(control_node, "time_stepper"; type="default")
aen(time_stepper_node, "max_retries", max_retries)
aen(time_stepper_node, "opt_iter", opt_iter)
aen(time_stepper_node, "dtmin", dtmin)
aen(time_stepper_node, "dtmax", dtmax)
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

material_node = aen(Material_node, "material"; id="1", name="Material1", type="neo-Hookean")
aen(material_node, "E", E_youngs)
aen(material_node, "v", ν)

material_node = aen(Material_node, "material"; id="2", name="Material2", type="rigid body")
aen(material_node, "density", 1.0)
aen(material_node, "center_of_mass", join([@sprintf("%.16e", x) for x ∈ mean(V2)], ','))

material_node = aen(Material_node, "material"; id="3", name="Material3", type="rigid body")
aen(material_node, "density", 1.0)
aen(material_node, "center_of_mass", join([@sprintf("%.16e", x) for x ∈ mean(V3)], ','))

material_node = aen(Material_node, "material"; id="4", name="Material4", type="rigid body")
aen(material_node, "density", 1.0)
aen(material_node, "center_of_mass", join([@sprintf("%.16e", x) for x ∈ mean(V4)], ','))

# Mesh     
Mesh_node = aen(febio_spec_node, "Mesh")

# Nodes
Nodes_node = aen(Mesh_node, "Nodes"; name="nodeSet_all")
for (i, v) in enumerate(V)
    aen(Nodes_node, "node", join([@sprintf("%.16e", x) for x ∈ v], ','); id=@sprintf("%i", i))
end

# Elements
Elements_node = aen(Mesh_node, "Elements"; name="Part1", type="tet4")
for (i, e) in enumerate(E)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", i))
end

Elements_node = aen(Mesh_node, "Elements"; name="Part2", type="tri3")
for (i, e) in enumerate(F2)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", i+length(E)))
end

Elements_node = aen(Mesh_node, "Elements"; name="Part3", type="tri3")
for (i, e) in enumerate(F3)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", i+length(E)+length(F2)))
end

Elements_node = aen(Mesh_node, "Elements"; name="Part4", type="tri3")
for (i, e) in enumerate(F4)
    aen(Elements_node, "elem", join([@sprintf("%i", i) for i ∈ e], ", "); id=@sprintf("%i", i+length(E)+length(F2)+length(F3)))
end

surfaceName1 = "Surface_Bone"
Surface_node = aen(Mesh_node, "Surface"; name=surfaceName1)
for (i, e) in enumerate(Fb)
    aen(Surface_node, "tri3", join([@sprintf("%i", j) for j in e], ','); id=@sprintf("%i", i))
end

surfaceName2 = "Surface_cyl1"
Surface_node = aen(Mesh_node, "Surface"; name=surfaceName2)
for (i, e) in enumerate(F2)
    aen(Surface_node, "tri3", join([@sprintf("%i", j) for j in e], ','); id=@sprintf("%i", i))
end

surfaceName3 = "Surface_cyl2"
Surface_node = aen(Mesh_node, "Surface"; name=surfaceName3)
for (i, e) in enumerate(F3)
    aen(Surface_node, "tri3", join([@sprintf("%i", j) for j in e], ','); id=@sprintf("%i", i))
end

surfaceName4 = "Surface_cyl4"
Surface_node = aen(Mesh_node, "Surface"; name=surfaceName4)
for (i, e) in enumerate(F4)
    aen(Surface_node, "tri3", join([@sprintf("%i", j) for j in e], ','); id=@sprintf("%i", i))
end

surfacePairName1 = "SurfacePair_Cyl1_Bone"
SurfacePair_node = aen(Mesh_node, "SurfacePair"; name=surfacePairName1)
aen(SurfacePair_node, "primary", surfaceName1)
aen(SurfacePair_node, "secondary", surfaceName2)

surfacePairName2 = "SurfacePair_Cyl2_Bone"
SurfacePair_node = aen(Mesh_node, "SurfacePair"; name=surfacePairName2)
aen(SurfacePair_node, "primary", surfaceName1)
aen(SurfacePair_node, "secondary", surfaceName3)

surfacePairName3 = "SurfacePair_Cyl3_Bone"
SurfacePair_node = aen(Mesh_node, "SurfacePair"; name=surfacePairName3)
aen(SurfacePair_node, "primary", surfaceName1)
aen(SurfacePair_node, "secondary", surfaceName4)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
aen(MeshDomains_node, "SolidDomain"; mat="Material1", name="Part1")

shelldomain_node = aen(MeshDomains_node, "ShellDomain"; mat="Material2", name="Part2")
# aen(shelldomain_node,"shell_thickness",0.1)

shelldomain_node = aen(MeshDomains_node, "ShellDomain"; mat="Material3", name="Part3")
# aen(shelldomain_node,"shell_thickness",0.1)

shelldomain_node = aen(MeshDomains_node, "ShellDomain"; mat="Material4", name="Part4")
# aen(shelldomain_node,"shell_thickness",0.1)

Rigid_node = aen(febio_spec_node, "Rigid")
rigid_bc = aen(Rigid_node, "rigid_bc"; name="RigidFixRot_cyl2", type="rigid_fixed")
aen(rigid_bc, "rb", 3)
aen(rigid_bc, "Rx_dof", 1)
aen(rigid_bc, "Ry_dof", 1)
aen(rigid_bc, "Rz_dof", 1)
aen(rigid_bc, "Ru_dof", 1)
aen(rigid_bc, "Rv_dof", 1)
aen(rigid_bc, "Rw_dof", 1)

Rigid_node = aen(febio_spec_node, "Rigid")
rigid_bc = aen(Rigid_node, "rigid_bc"; name="RigidFixRot_cyl3", type="rigid_fixed")
aen(rigid_bc, "rb", 4)
aen(rigid_bc, "Rx_dof", 1)
aen(rigid_bc, "Ry_dof", 1)
aen(rigid_bc, "Rz_dof", 1)
aen(rigid_bc, "Ru_dof", 1)
aen(rigid_bc, "Rv_dof", 1)
aen(rigid_bc, "Rw_dof", 1)

Rigid_node = aen(febio_spec_node, "Rigid")
rigid_bc = aen(Rigid_node, "rigid_bc"; name="RigidFixRot_cyl1", type="rigid_fixed")
aen(rigid_bc, "rb", 2)
aen(rigid_bc, "Rx_dof", 1)
aen(rigid_bc, "Ry_dof", 1)
aen(rigid_bc, "Rz_dof", 0)
aen(rigid_bc, "Ru_dof", 1)
aen(rigid_bc, "Rv_dof", 1)
aen(rigid_bc, "Rw_dof", 1)

rigid_bc = aen(Rigid_node, "rigid_bc"; name="RigidPrescribe_cyl1_Z", type="rigid_displacement")
aen(rigid_bc, "rb", 2)
aen(rigid_bc, "dof", "z")
aen(rigid_bc, "value", cylDisplacement_Z; lc="1")
aen(rigid_bc, "relative", 0)

for surfacePairName in (surfacePairName1, surfacePairName2, surfacePairName3)
    Contact_node = aen(febio_spec_node, "Contact")
    contact_node = aen(Contact_node, "contact"; type="sticky", surface_pair=surfacePairName)
    aen(contact_node, "penalty", contactPenalty)
    aen(contact_node, "laugon", laugon)
    aen(contact_node, "tolerance", 0.2)
    aen(contact_node, "minaug", minaug)
    aen(contact_node, "maxaug", maxaug)
    aen(contact_node, "snap_tol", 0.01)
    aen(contact_node, "max_traction", 0.01)
    aen(contact_node, "search_tolerance", 2.0*pointSpacing)
end

LoadData_node = aen(febio_spec_node, "LoadData")
load_controller_node = aen(LoadData_node, "load_controller"; id="1", name="LC_1", type="loadcurve")
aen(load_controller_node, "interpolate", "LINEAR")
points_node = aen(load_controller_node, "points")
aen(points_node, "pt", @sprintf("%.2f, %.2f", 0.0, 0.0))
aen(points_node, "pt", @sprintf("%.2f, %.2f", 0.5, 0.5))
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

# #######
# Write FEB file
XML.write(filename_FEB, doc)

# #######
# Run FEBio
run_febio(filename_FEB, FEBIO_EXEC)

# #######
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

function get_elementData_limits(D)
    s1_max = -Inf
    s1_min = Inf
    @inbounds for i in 0:1:(length(D)-1)
        S = [s[1] for s in D[i].data]
        s1_max = max(s1_max, maximum(S))
        s1_min = min(s1_min, minimum(S))
    end
    return s1_min, s1_max
end
s1_min, s1_max = get_elementData_limits(DD_stress)
min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

## Visualization

fig = Figure(size=(1200, 800))
stepStart = incRange[end]
ax1 = AxisGeom(fig[1, 1], title="Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=0.0, color=UT_mag[end], transparency=false, colormap=Reverse(:Spectral), colorrange=(0, maximum(ut_mag_max)))
hp2 = meshplot!(ax1, F2, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)
hp3 = meshplot!(ax1, F3, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)
hp4 = meshplot!(ax1, F4, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)

Colorbar(fig[1, 2], hp1, label="Displacement magnitude [mm]")

S_E = [s[1] for s in DD_stress[stepStart].data]
S_F = repeat(S_E[1:length(E)], inner=4)
indB = boundaryfaceindices(F)

# Fs,Vs = separate_vertices(F[indB],VT[stepStart+1])Cbs
S_Vs = simplex2vertexdata(F[indB], S_F[indB], V)

ax2 = AxisGeom(fig[1, 3], title="Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp5 = meshplot!(ax2, F[indB], VT[end]; strokewidth=0.0, color=S_Vs, colormap=:viridis, colorrange=(s1_min/30.0, s1_max/30.0))
hp6 = meshplot!(ax2, F2, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)
hp7 = meshplot!(ax2, F3, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)
hp8 = meshplot!(ax2, F4, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)

Colorbar(fig[1, 4], hp5, label="S1")

hSlider = Slider(fig[2, :], range=incRange, startvalue=stepStart, linewidth=30)
on(hSlider.value) do stepIndex
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex+1], Fb)
    hp1.color = UT_mag[stepIndex+1]

    hp2[1] = GeometryBasics.Mesh(VT[stepIndex+1], F2)
    hp3[1] = GeometryBasics.Mesh(VT[stepIndex+1], F3)
    hp4[1] = GeometryBasics.Mesh(VT[stepIndex+1], F4)

    ax1.title = "Step: $stepIndex"
    S_E = [s[1] for s in DD_stress[stepIndex].data]
    S_F = repeat(S_E[1:length(E)], inner=4)
    # Fs,Vs = separate_vertices(F[indB],VT[stepIndex+1])
    S_Vs = simplex2vertexdata(F[indB], S_F[indB], V)

    hp5[1] = GeometryBasics.Mesh(VT[stepIndex+1], F[indB])
    hp5.color = S_Vs
    hp6[1] = GeometryBasics.Mesh(VT[stepIndex+1], F2)
    hp7[1] = GeometryBasics.Mesh(VT[stepIndex+1], F3)
    hp8[1] = GeometryBasics.Mesh(VT[stepIndex+1], F4)
    ax2.title = "Step: $stepIndex"
end
slidercontrol(hSlider, ax1)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")