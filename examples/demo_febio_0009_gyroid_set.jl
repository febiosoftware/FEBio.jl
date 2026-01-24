using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.LinearAlgebra
using Comodo.Statistics
using Geogram # https://github.com/COMODO-research/Geogram.jl
using FEBio
using FEBio.XML
using Printf 

GLMakie.closeall()

######
# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

###### 
# Control parameters 
E_youngs = 1.0
ν = 0.3

displacement_prescribed = -2.0

###### 
pointSpacing = 0.5
level = 0.75
s1 = -0.1*π
s2 =  4.1*π
ns = spacing2numsteps(s2-s1, pointSpacing/3.0; close_loop=false)
x = range(s1, s2, ns)

type=:G
F1, V1 = tpms(type; x=x, level= -level, cap = true, padValue=1e8, side=:positive)
F2, V2 = tpms(type; x=x, level= level, cap = true, padValue=1e8, side=:negative)

np1 = spacing2numvertices(F1, V1, pointSpacing)
F1,V1 = ggremesh(F1,V1; nb_pts=np1)

np2 = spacing2numvertices(F2, V2, pointSpacing)
F2,V2 = ggremesh(F2,V2; nb_pts=np2)

V1_in = faceinteriorpoint(F1,V1,1)
V2_in = faceinteriorpoint(F2,V2,1)

Fb, Vb, Cb = joingeom(F1, V1, F2, V2)

V_regions = [V1_in, V2_in]

vol1 = pointSpacing^3 / (6.0*sqrt(2.0))
region_vol = [vol1, vol1]
stringOpt = "paAqYQ"
E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=V_regions, region_vol=region_vol, stringOpt)

## Visualization ------------------------------------------------------------
cmap = cgrad(:Spectral, 2, categorical = true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4)

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)

Fs,Vs = separate_vertices(F,V)
CE_Vs = simplex2vertexdata(Fs,CE_F)
M = GeometryBasics.Mesh(Vs,Fs)

strokewidth = 1 

fig1 = Figure(size=(800,800))

ax1 = AxisGeom(fig1[1, 1][1, 1], title = "Boundary surfaces")
hp1 = meshplot!(ax1, Fb[Cb.==1], V, color=:brown, strokewidth=0.5)
hp1 = meshplot!(ax1, Fb[Cb.==2], V, color=:white, strokewidth=0.5)

scatter!(ax1, V_regions, color=:black, markersize=25)

ax2 = AxisGeom(fig1[1, 1][1, 2], title = "Cut mesh")
hp2 = meshplot!(ax2, Fs, Vs, color=CE_Vs, strokewidth=strokewidth,colorrange = (1,2), colormap=cmap)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider1 = Slider(fig1[2, 1], range = stepRange, startvalue = mean(stepRange),linewidth=30)

on(hSlider1.value) do z 
    B = ZE .<= z
    indShow = findall(B)
    if isempty(indShow)
        hp2.visible=false        
    else        
        hp2.visible=true
        Fs = element2faces(E[indShow])
        Cs = repeat(CE[indShow],inner=4)
        
        indB = boundaryfaceindices(Fs)        
        Fs = Fs[indB]
        Cs = Cs[indB]
        Fs,Vs = separate_vertices(Fs,V)
        CE_Vs = simplex2vertexdata(Fs,Cs)
        Ms = GeometryBasics.Mesh(Vs,Fs)
        hp2[1] = Ms
        hp2.color = CE_Vs
    end
end

screen = display(GLMakie.Screen(), fig1)
GLMakie.set_title!(screen, "FEBio example")
## ----------------------------------------------------------------------------

# Find boundary condition faces 
indicesTopNodes = Vector{Int}()
indicesBottomNodes = Vector{Int}()
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
search_tol = pointSpacing/2.0
for f in Fb    
    bTop = Z[f].>(zMax-search_tol)
    bBottom = Z[f].<(zMin+search_tol)
    if all(bTop)
        append!(indicesTopNodes, f)        
    end
    if all(bBottom)
        append!(indicesBottomNodes, f)        
    end
end
indicesTopNodes = unique(indicesTopNodes)
indicesBottomNodes = unique(indicesBottomNodes)

## Visualization ------------------------------------------------------------
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title = "Boundary surfaces")
hp1 = meshplot!(ax1, Fb, V, color=:white)
scatter!(ax1, V[indicesTopNodes], color=:red, markersize=10)
scatter!(ax1, V[indicesBottomNodes], color=:blue, markersize=10)
screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")
## ----------------------------------------------------------------------------

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
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="tet4")
    for (i,e) in enumerate(E)     
        aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = i)
    end

# Node sets
bcPrescribeList = "bcPrescribeList_z"
bcSupportList = "bcSupportList"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesTopNodes   ],','); name=bcPrescribeList)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indicesBottomNodes],','); name=bcSupportList)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xyz", node_set=bcSupportList, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xy", node_set=bcPrescribeList, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

bc_node = aen(Boundary_node,"bc"; name="prescribed_disp_z", node_set=bcPrescribeList, type="prescribed displacement")
    aen(bc_node,"dof","z")
    aen(bc_node,"value", displacement_prescribed; lc=@sprintf("%i",1))
    aen(bc_node,"relative",@sprintf("%i",0))

LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f", 0.0, 0.0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f", 1.0, 1.0))

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
    aen(logfile_node,"element_data"; data="s1", delim=",", file=filename_stress)
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
S1 = [d[1] for d in DD_stress[numInc-1].data] 
S1_F = repeat(S1,inner=6)
S1_Fb = S1_F[indBoundaryFaces]
S1_V = simplex2vertexdata(Fb, S1_Fb,V)

# Create z-force curve data 
Fz_curve = zeros(numInc)
time_curve = zeros(numInc)
for i in 0:1:numInc-1    
    Fz = [d[3] for d in DD_force[i].data]
    Fz_curve[i+1] = sum(Fz[indicesTopNodes])
    time_curve[i+1] = DD_force[i].time
end

#######
# Visualization

# Create time varying vectors
UT = fill(V,numInc) 
VT = fill(V,numInc)
S1T = fill(S1_V,numInc)
UT_mag = fill(zeros(length(V)),numInc)
ut_mag_max = zeros(numInc)
S1_max = zeros(numInc)
S1_min = zeros(numInc)
@inbounds for i in 0:1:numInc-1    
    UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
    VT[i+1] += UT[i+1]
    UT_mag[i+1] = norm.(UT[i+1])
    ut_mag_max[i+1] = maximum(UT_mag[i+1]) 

    S1 = [d[1] for d in DD_stress[i].data]
    S1_F = repeat(S1,inner=6)
    S1_Fb = S1_F[indBoundaryFaces]
    S1_V = simplex2vertexdata(Fb, S1_Fb, V)
    S1T[i+1] = S1_V
    S1_min[i+1] = minimum(S1)
    S1_max[i+1] = maximum(S1)
end

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])


fig2 = Figure(size=(1200, 800))
stepStart = incRange[end]
ax1 = AxisGeom(fig2[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=1, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral), colorrange=(0.0, maximum(ut_mag_max)))
Colorbar(fig2[1, 2], hp1.plots[1], label = "Displacement magnitude [mm]") 

ax2 = AxisGeom(fig2[1, 3], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp2 = meshplot!(ax2, Fb, VT[end]; strokewidth=1, color=S1T[end], transparency=false, colormap = :viridis, colorrange=(minimum(S1_min), maximum(S1_max)/10.0))
Colorbar(fig2[1, 4], hp2.plots[1], label = "σ [MPa]") 

ax3 = Axis(fig2[1, 5], title = "Step: $stepStart", aspect = AxisAspect(1), xlabel="Time [s]", ylabel="Force [N]")
lines!(ax3, time_curve, Fz_curve, color=:red, linewidth=3)
hp3 = scatter!(ax3, Point{2,Float64}(time_curve[stepStart+1], Fz_curve[stepStart+1]), markersize=15, color=:red)
hSlider2 = Slider(fig2[2, :], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider2.value) do stepIndex 
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp1.color = UT_mag[stepIndex+1]

    hp2[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp2.color = S1T[stepIndex+1]
        
    hp3[1] = Point{2,Float64}(time_curve[stepIndex+1], Fz_curve[stepIndex+1])
    ax1.title = "Step: $stepIndex"
    ax2.title = "Step: $stepIndex"
    ax3.title = "Step: $stepIndex"
end

screen = display(GLMakie.Screen(), fig2)
GLMakie.set_title!(screen, "FEBio example")