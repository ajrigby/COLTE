# COLTE: CO Local Thermodynamic Equilibrium
A program that produces LTE-based column density cubes from 12CO, 13CO (&amp; C18O),
largely following the procedure outlined in [Rigby et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...632A..58R/abstract), which uses COHRS 12CO and CHIMPS 13CO and C18O 3-2 data.

## Installation
The easiest way to install COLTE is via pip
    
    pip install git+https://github.com/ajrigby/COLTE.git


## Running COLTE
COLTE is very simple to run, and simply requires a parameter file as the input.

    from colte.LTE import make_cubes

    make_cubes('params.json')

It requires that you have 12CO, 13CO, and C18O cubes of the same transition in
the main beam brightness temperature scale that are spatially and spectrally
matched, with the same size.

The code produces three cubes: excitation temperature, 13CO optical depth, and
13CO column density.


### Parameters:
The performance is controlled by a bunch of parameters that are stored in a .json file.
The file format is fairly simple. Please note that all paths must be given relative to
the location in which you are running `colte` (not relative to the location of the 
parameter file itself). Some suggested default values are stored in the 
`defaultparams.json` file on the GitHub page, which should give you an idea of how it works.

`data_in_12CO`: string. Path to the input 12CO cube.                                       
`data_in_13CO`: string. Path to the input 13CO cube.                                          
`data_in_C18O`: string. Path to the input C18O cube.                                          

`data_out`: string. A prefix for the output files. 

`transition`: string. The transition being used.                                                  

`snrlim12`: float. Signal-to-noise ratio threshold for 12CO.                                                           
`snrlim13`: float. Signal-to-noise ratio threshold for 13CO.                                                         
`snrlim18`: float. Signal-to-noise ratio threshold for C18O.                                                          

`Tex_mode`: integer. A value of 1 will calculate Tex primarily from 12CO.
                     A value of 2 will adopt single_Tex for all pixels
`single_Tex`: float. The assumed single  excitation temperature.             
`Tau_mode`: integer. A value of 1 will calculate Tau primarly from 13CO/Tex
                     A value of 2 will use the low tau assumption with single_Tau
                     A value of 3 will use the optically thin approximation.  
`R1318`: integer. Assumeed abundance ratio of 13CO to C18O. We refer the reader to [Wilson & Rood 1994](https://ui.adsabs.harvard.edu/abs/1994ARA%26A..32..191W/abstract) for selecting suitable values.  
`single_Tau:` float. The assumed optical depth for use in Tau_mode 2.  
`fill_method`: string. Either "interpolate", "median", "rescale" or "none".     

#### Note on fill methods:
The fill_method parameter controls how tau values are filled, which often arise as a result of the 13CO brightness temperature approaching or exceed the value for 12CO, which may happen due to non-LTE conditions, such as self absoprtions, or line of sight temperature gradients.

`"interpolate"`: Uses a linear interpolation to fill missing values.  
`"median"`: Uses a median filter to fill missing values.  
`"rescale"`: First fills tau values based on the average 13CO emission to optical depth relationship, and then median filters the map to smooth out outlying values.  
`"none"`: will not perform any filling of tau values.
