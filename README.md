![](assets/img/febio_jl_logo_banner.png)

# FEBio.jl: A julia wrapper for FEBio
This projects aims to provide a Julia interface to the [FEBio](https://febio.org/) finite element solver. 

## Installation
```julia
pkg> add https://github.com/febiosoftware/FEBio.jl
```

## Configuration
No special configuration is required at present. However, users do need to specify the path to the FEBio executable at the top of their code, e.g.: 
```julia
const FEBIO_PATH = "/home/kevin/FEBioStudio/bin/febio4" # Path to FEBio executable
```
Next febio can be called to run the analysis for the input file `filename_FEB` using: 
```julia
# Run FEBio
runMonitorFebio(filename_FEB,FEBIO_PATH)
```

## Getting started
The `docs` folder contains usage examples. See for instance: `demo_febio_input_file_01.jl`. This example should produce:

![](assets/img/febio_example_01.gif)

## Documentation 
Under construction

## Testing 
Under construction

## Roadmap
A detailed roadmap is under construction but the below is a list of major components which are currently been worked on: 

- [x] Julia based .feb file creation, import, and editing  
- [x] Julia calls/triggers FEBio executable
- [ ] Julia monitors/polices FEBio simulation progress
- [x] Julia based log-file importing 
- [ ] Julia based .xplt-file importing 
- [ ] Shipw with or download FEBio binaries upon package adding, and configure with FEBio.jl automatically  

## How to contribute? 
Your help would be greatly appreciated! If you can contribute please do so by posting a pull-request. 

## License 
This wrapper is released open source under the [Apache 2.0 license](https://github.com/febiosoftware/FEBio.jl/blob/main/LICENSE).