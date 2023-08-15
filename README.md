# COLTE: CO Local Thermodynamic Equilibrium
A program that produces LTE-based column density cubes from 12CO, 13CO (&amp; C18O).

# Running COLTE
COLTE is very simple to run.

import colte

colte.make_LTE_cubes('default.params')

The performance is controlled by a bunch of parameters that are stored in a .json file.

# Parameters:
data_in_12CO: string. Path to the input 12CO cube. ["tmb_12CO.fits"]                                         
data_in_13CO: string. Path to the input 13CO cube. ["tmb_13CO.fits"]                                            
data_in_C18O: string. Path to the input C18O cube. ["tmb_C18O.fits"]                                            
single_Tex: float. Give it a value of 0.0 to calculate Tex from 12CO, 
                   or a non-zero value will provide an assumed single 
                   constant excitation temperature.                           
data_out: string. A prefix for the output files. ["colte_"]                                                    
transition: string. The transition being used with the format "Jup-Jlo". ["3-2"]                                                        
snrlim12: float. Signal-to-noise ratio threshold for 12CO. [3.0]                                                            
snrlim13: float. Signal-to-noise ratio threshold for 13CO. [3.0]                                                            
snrlim18: float. Signal-to-noise ratio threshold for C18O. [5.0]                                                            
fill_method: string. Either "interpolate" or "median". ["interpolate"]       

## Fill methods:
The fill_method parameter controls how tau values are filled, which often arise as a result of the 13CO brightness temperature approaching or exceed the value for 12CO, which may happen due to non-LTE conditions, such as self absoprtions, or line of sight temperature gradients.
interpolate: Uses a linear interpolation to fill missing values
median: Uses a median filter to fill missing values.
extrapolate: First fills tau values based on the average 13CO emission to optical depth relationship, and then median filters the map to smooth out outlying values.
