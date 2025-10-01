using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.LinearAlgebra
using Comodo.Statistics
using FEBio
using FEBio.XML
using Printf 

GLMakie.closeall()

######
# Set FEBio exec path or name
const FEBIO_EXEC = "febio4" # FEBio executable

###### 
# Control parameters 
pointSpacing = 4.0

# Material parameters
c1 = 1e-3 #Shear-modulus-like parameter
m1 = 2  #Material parameter setting degree of non-linearity
k_factor = 1000.0  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus

# Contact parameters
contactPenalty = 100.0;
laugon = 0;
minaug = 1;
maxaug = 10;
fric_coeff = 0.01;
initialSpacing = pointSpacing/100.0

# FEA control settings
numTimeSteps=10 # Number of time steps desired
max_refs=50 # Max reforms
max_ups=0 # Set to zero to use full-Newton iterations
opt_iter=20 # Optimum number of iterations
max_retries=5 # Maximum number of retires
dtmin=(1.0/numTimeSteps)/100.0 # Minimum time step size
dtmax=1.0/numTimeSteps  # Maximum time step size
symmetric_stiffness=0 
min_residual=1e-18


###### 

hCylinder = 40.0
n = 2
rCylinder = 50.0
nh = 0
rSphere = rCylinder/2.0

E, V, F, Fb, Cb = hexcylinder(rCylinder, hCylinder, n; direction=:negative)  

F2, V2, C2 = quadsphere(rSphere, pointSpacing)

# Shift sphere indices as they are appended after solid mesh
vb = eltype(V2)(0.0, 0.0, rSphere+initialSpacing)
V2 .+= vb
F2 = [eltype(F2)(f.+length(V)) for f in F2]

# Append sphere nodes
append!(V, V2) # Add sphere nodes to set 

sphereDisplacement_XYZ = Point{3,Float64}(0.0, 0.0, -((0.5*hCylinder)-initialSpacing))
V2_centre = Point{3,Float64}(0.0, 0.0, 0.0) # mean(V2)

# Create face sets to define node sets later 
Fb_top = Fb[Cb.==2]
Fb_support = Fb[Cb.==1] # Bottom
# Fb_support = Fb[[in(c,[0,1]) for c in Cb]] # Bottom and sides

bcSupportList = elements2indices(Fb_support)

######
# ---------- Visualisation
Fbs, Vbs = separate_vertices(Fb,V)
Cbs = simplex2vertexdata(Fbs,Cb)

cmap = Makie.Categorical(:Spectral) 
fig = Figure(size=(1200,1000))
ax1 = AxisGeom(fig[1, 1], title="Hexahedral mesh of a cylinder, n=$n refinement steps")
hp1 = meshplot!(ax1, Fbs, Vbs; color=Cbs, colormap=cmap)
hp2 = meshplot!(ax1, F2, V; color=:white)
Colorbar(fig[1, 1][1, 2], hp1)

ax2 = AxisGeom(fig[1, 2], title="Boundary conditions")
hp3 = meshplot!(ax2, Fbs, Vbs; color=RGBA{Float64}(1.0, 1.0, 1.0, 0.5), transparency=true)
hp4 = meshplot!(ax2, F2, V; color=RGBA{Float64}(0.0, 1.0, 0.0, 0.5), transparency=true)
hp5 = meshplot!(ax2, F2, V.+sphereDisplacement_XYZ; color=RGBA{Float64}(1.0, 0.0, 0.0, 0.5), transparency=true)
hp6 = scatter!(ax2, V[bcSupportList], color=:black, markersize=10)
Legend(fig[1,3], [hp3, hp4, hp5, hp6],["Solid mesh", "Rigid body initial", "Rigid body final", "Supported nodes"])
screen = display(GLMakie.Screen(), fig)


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
    aen(control_node,"time_steps",numTimeSteps)
    aen(control_node,"step_size",1.0/numTimeSteps)
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)

time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",max_retries)
    aen(time_stepper_node,"opt_iter",opt_iter)
    aen(time_stepper_node,"dtmin",dtmin)
    aen(time_stepper_node,"dtmax",dtmax)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",symmetric_stiffness)
    aen(solver_node,"equation_scheme",1)
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    aen(solver_node,"max_refs",max_refs)
    aen(solver_node,"check_zero_diagonal",0)
    aen(solver_node,"zero_diagonal_tol",0)
    aen(solver_node,"force_partition",0)
    aen(solver_node,"reform_each_time_step",1)
    aen(solver_node,"reform_augment",0)
    aen(solver_node,"diverge_reform",1)
    aen(solver_node,"min_residual",min_residual)
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
    aen(qn_method_node,"max_ups",max_ups)
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

material_node = aen(Material_node,"material"; id = 2, name="Material2", type="rigid body")
    aen(material_node,"density",1.0)
    aen(material_node,"center_of_mass",join([@sprintf("%.16e",x) for x ∈ V2_centre],','))

# Mesh     
Mesh_node = aen(febio_spec_node,"Mesh")

    # Nodes
    Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
        for q ∈ eachindex(V)            
            aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ V[q]],','); id = q)     
        end
        
    # Elements
    Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
        for (i,e) in enumerate(E)     
            aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = i)
        end
        
    Elements_node = aen(Mesh_node,"Elements"; name="Part2", type="quad4")
        for (i,e) in enumerate(F2)     
            aen(Elements_node,"elem", join([@sprintf("%i", i) for i ∈ e], ", "); id = i)
        end    

    # Node sets
    bcSupportListName = "bcSupportList"
        aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ bcSupportList],','); name=bcSupportListName)

    surfaceName1 = "Surface_Top"
    Surface_node = aen(Mesh_node,"Surface"; name=surfaceName1)
    for (i,e) in enumerate(Fb_top)        
        aen(Surface_node,"quad4",join([@sprintf("%i",j) for j in e], ','); id = i)
    end

    surfaceName2 = "Surface_Sphere"
    Surface_node = aen(Mesh_node,"Surface"; name=surfaceName2)
    for (i,e) in enumerate(F2)        
        aen(Surface_node,"quad4",join([@sprintf("%i",j) for j in e], ','); id = i)
    end

    surfacePairName1 = "SurfacePair_Top_Sphere"
    SurfacePair_node = aen(Mesh_node,"SurfacePair"; name=surfacePairName1)
        aen(SurfacePair_node, "primary", surfaceName1)
        aen(SurfacePair_node, "secondary", surfaceName2)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

    shelldomain_node = aen(MeshDomains_node,"ShellDomain"; mat = "Material2", name="Part2")
    # aen(shelldomain_node,"shell_thickness",0.1)
    
Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xyz", node_set=bcSupportListName, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

Rigid_node = aen(febio_spec_node, "Rigid")
    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidFixRot_sphere1", type="rigid_fixed")
        aen(rigid_bc, "rb", 2)
        aen(rigid_bc, "Ru_dof", 1)
        aen(rigid_bc, "Rv_dof", 1)
        aen(rigid_bc, "Rw_dof", 0)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_sphere1_X", type="rigid_displacement")
        aen(rigid_bc, "rb", 2)
        aen(rigid_bc, "dof", "x")
        aen(rigid_bc, "value", sphereDisplacement_XYZ[1]; lc=1)
        aen(rigid_bc, "relative", 0)
    
    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_sphere1_Y", type="rigid_displacement")
        aen(rigid_bc, "rb", 2)
        aen(rigid_bc, "dof", "y")
        aen(rigid_bc, "value", sphereDisplacement_XYZ[2]; lc=1)
        aen(rigid_bc, "relative", 0)

    rigid_bc = aen(Rigid_node,"rigid_bc"; name="RigidPrescribe_sphere1_Z", type="rigid_displacement")
        aen(rigid_bc,"rb",2)
        aen(rigid_bc,"dof","z")
        aen(rigid_bc,"value", sphereDisplacement_XYZ[3]; lc=1)
        aen(rigid_bc,"relative",0)

Contact_node = aen(febio_spec_node,"Contact")
    contact_node = aen(Contact_node,"contact"; type="sliding-elastic", surface_pair=surfacePairName1)
        aen(contact_node,"two_pass",0)
        aen(contact_node,"laugon",laugon)
        aen(contact_node,"tolerance",0.2)
        aen(contact_node,"gaptol",0.0)
        aen(contact_node,"minaug",minaug)
        aen(contact_node,"maxaug",maxaug)
        aen(contact_node,"search_tol",0.01)
        aen(contact_node,"search_radius", 2.0*pointSpacing)
        aen(contact_node,"symmetric_stiffness",0)
        aen(contact_node,"auto_penalty",1)
        aen(contact_node,"penalty",contactPenalty)
        aen(contact_node,"fric_coeff",fric_coeff)
            
LoadData_node = aen(febio_spec_node,"LoadData")
    load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
        aen(load_controller_node,"interpolate","LINEAR")
        points_node = aen(load_controller_node,"points")
            aen(points_node,"pt",@sprintf("%.2f, %.2f", 0.0, 0.0))
            aen(points_node,"pt",@sprintf("%.2f, %.2f", 0.5, 0.5))
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
    aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)

# #######
# Write FEB file
XML.write(filename_FEB, doc)

# #######
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

function get_elementData_limits(D)
    s1_max = -Inf
    s1_min = Inf
    @inbounds for i in 0:1:length(D)-1            
        S = [s[1] for s in D[i].data]
        s1_max = max(s1_max, maximum(S))
        s1_min = min(s1_min, minimum(S))
    end
    return s1_min, s1_max
end
s1_min, s1_max = get_elementData_limits(DD_stress)
min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

#######
# Visualization

fig = Figure(size=(1200,800))
stepStart = incRange[end]
ax1 = AxisGeom(fig[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=1.0, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral), colorrange=(0,maximum(ut_mag_max)))
hp2 = meshplot!(ax1, F2, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)

Colorbar(fig[1, 2], hp1, label = "Displacement magnitude [mm]") 

S_E = [s[1] for s in DD_stress[stepStart].data]
S_F = repeat(S_E,inner=6)

Fs,Vs = separate_vertices(F,VT[stepStart+1])
S_Vs = simplex2vertexdata(Fs,S_F)

ax2 = AxisGeom(fig[1, 3], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp3 = meshplot!(ax2, Fs, Vs; strokewidth=1.0, color=S_Vs , colormap = :viridis, colorrange=(s1_min, s1_max))
hp4 = meshplot!(ax2, F2, VT[end]; strokewidth=0.0, color=RGBA{Float64}(1.0, 1.0, 1.0, 0.25), transparency=true)

Colorbar(fig[1, 4], hp3, label = "S1") 

hSlider = Slider(fig[2, :], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider.value) do stepIndex 
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp1.color = UT_mag[stepIndex+1]
    hp2[1] = GeometryBasics.Mesh(VT[stepIndex+1],F2)
    ax1.title = "Step: $stepIndex"

    S_E = [s[1] for s in DD_stress[stepIndex].data]
    S_F = repeat(S_E,inner=6)
    Fs,Vs = separate_vertices(F,VT[stepIndex+1])
    S_Vs = simplex2vertexdata(Fs,S_F)

    hp3[1] = GeometryBasics.Mesh(Vs, Fs)
    hp4[1] = GeometryBasics.Mesh(VT[stepIndex+1],F2)
    hp3.color = S_Vs    
    ax2.title = "Step: $stepIndex"
end

slidercontrol(hSlider,ax)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")