# AbsorptionTMM
We intend to use the transfer matrix representation to calculate the local absorbed energy density as a function of space. 
This way multiple reflections and interferances in thin films, different incident angles and polarizations, s & p, are considered.
The description of the formulas and the way we are implementing them can be found here: [Description](https://github.com/udcm-su/AbsorptionTMM/blob/master/Transfermatrix_Description.pdf)

  --- 
### Use
The [TMM_abs.py](https://github.com/udcm-su/AbsorptionTMM/blob/master/TMM_abs.py) file works as function liberary and the [Cobalt_example.py](https://github.com/udcm-su/AbsorptionTMM/blob/master/Cobalt_example.py) is a file which gives an example on how to use the code. 
### An example

  
| Platinum       	| 3nm layer    	| n = 1.0433 + i3.0855 	|
|----------------	|--------------	|----------------------	|
| Cobalt         	| 15nm layer   	| n = 1.0454 + i3.2169 	|
| Chromium       	| 5nm layer    	| n = 2.0150 + i2.8488 	|
| Magnesium Oxid 	| semi inf- layer 	| n = 1.7660           	|


  
  <img src="https://github.com/udcm-su/AbsorptionTMM/blob/master/5Layer.png" width="620" height="500" />



