using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.LinearAlgebra
using Comodo.Statistics
using FEBio
using FEBio.XML
using Printf 

#=
This demo aims to show how to define the XML structure for the FEBio .feb 
input files. In particular the creation of the following is shown: 

<?xml version="1.0" encoding="UTF-8"?>
<febio_spec version="4.0">
  <!--Created using Julia-->
  <Material>
    <material id="1" name="Material1" type="Neo-Hookean">
      <E>1.0</E>
      <v>0.45</v>
    </material>
  </Material>
</febio_spec>

As can be seen the main XML "node" is called febio_spec. The term "node" should 
not be confused with a finite element node. Here "node" refers to a level in the 
XML file. It may be helpful to think of a file folder structure on a computer 
for instance `febio_spec` is a main folder containing a folder called `Material` 
which will contain each `material`. Besides these different levels we can see 
that some nodes have "attributes". Attributes say something about that node. For 
instance the attribute `id="1"` for the material. Attributes are different from 
"data" for instance the `1.0` in `<E>1.0</E>` is not an attribute but the data 
or value for the parameter `E`. 

Once one is familiar with the basic XML file structure and how values and 
attributes are different, one is ready to learn how to define the XML using 
code.
 =#

# First a file name is defined for the XML (.feb) file
saveDir = joinpath(febiojl_dir(),"assets","temp") # Folder for saving file
if !isdir(saveDir)
    mkdir(saveDir)      
end

# Define file name
filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file

#=
 Next we can initiate an XML object for the .feb file
=#

doc,febio_spec_node = feb_doc_initialize() # Start with FEBio feb file template

#=
Here `doc` is the XML file and `febio_spec_node` a handle to the `febio_spec` 
node or level. 
=#

# Next add a `Material` node 
Material_node = aen(febio_spec_node,"Material")

#=
Here `aen` is short for `add element node`, but one can just remember it as a 
way to add a new level. As can be seen this `Material` node has no attributes or 
"data", it is simply defined and next we simply move on to the next level. 
Now `Material_node` here is a handle to this level and can be used to add 
things to it. 
=#

#= 
Next we need a `material` node which does have attributes, like `id` and 
`name`. We once again use `aen` but this time the attribrutes are added as 
keyword arguments. Note that the first input to `aen` is the level where the 
new node is inserted, here `Material_node` is used. 
Note how the `material` level does not have data, just attributes. 
=#

# Add a material with "attributes" (e.g. id) to the Material section 
material_node = aen(Material_node,"material"; id = 1, name="Material1", type="Neo-Hookean")

#= 
Finally the parameters `E` and `v` are added to the `material` level. We again 
use `aen` but this time these parameters do have data, i.e. the material 
constant values. The data for a node is added not as a keyword argument but 
simply as the third normal input argument for `aen`.  
=#

# Add parameters for the material 
aen(material_node,"E", E)
aen(material_node,"v", v)

#=
The above shows an example of creating different levels, e.g. the main 
`febio_spec` level, next the `Material` level, and next it's `material` level. 
It is also shown how any level can have attributes and/or data added. Using the
above, users can code any XML file structure they like. Typically also the FEBio
user manual shows an example of what the XML file entry would look like so users
can use that for inspiration. 

The final step is .feb file creation. 
=#

# Write FEB file
XML.write(filename_FEB, doc)
