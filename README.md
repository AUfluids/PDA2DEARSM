# PDA2DEARSM
Progressively data-augmented explicit algebraic Reynolds stress model (PDA-EARSM) enhanced for the prediction of secondary flows based on RANS $k-\omega$ SST turbulence model.
Developed by Fluid Mechanics & Turbulence research group at Aarhus University.

PDA2DEARSM - Implementation of the PDA2DEARSM RANS model
         as proposed by Rincón and Amarloo (2023) for OpenFOAM.

## License
    This program is free software: you can redistribute and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Description
Implementation of the progressively data-augmented explicit algebraic Reynolds stress model (PDAEARSM)
enhanced for the prediction of secondary flows based on standard RANS $k-\omega$ SST.
The model has been developed for 2D flows but its applicability has been tested on 3D
flows, yielding similar results.
Three coefficients can be modified by the user to change the model's performance.
Standard optimised values are given by default in the model.

## Target platform
The code is known to work with OpenFOAM v2112 and v2212.

## Authors
Mario Javier Rincón <mjrp@mpe.au.dk>
Ali Amarloo <amarloo@mpe.au.dk>

## Instructions

1. Download the source code using git:

         git clone https://github.com/AUfluids/PDA2DEARSM.git

2. Enter the directory where the source code has been extracted, and compile it by typing: 

         wmake

3. Add the following line to the _controlDict_ of your case:

         libs ( "libPDA2DEARSMIncompressibleTurbulenceModels" ) ;

4. Specify

         RASModel kOmegaSSTPDA2DEARSM;

in _turbulentProperties_.

5. Add the subdictionary

         PDA2DEARSM 2;  // 0: off | 1: ModelI | 2: ModelII

to _turbulentProperties_.

NOTE: You might have to define the bijDelta term in system/fvSchemes file, here you have an example:

         divSchemes
         {
                  div(dev(((2*k)*bijDelta)))          Gauss linear;
         }

## Test results

For more details, refer to the publication at: 
[Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation](https://doi.org/10.1016/j.ijheatfluidflow.2023.109242)
Or [ArXiv](https://arxiv.org/abs/2308.12720)
Qualitative results for a duct flow of aspect ratio 1 at bulk Reynolds number 3500 for Model **II**:
![alt text](https://github.com/AUfluids/PDA2DEARSM/blob/main/testCases/ductFlowAR1Reb3500/SD_u.png)
Quantitative results:
![alt text](https://github.com/AUfluids/PDA2DEARSM/blob/main/testCases/ductFlowAR1Reb3500/SD_profiles.png)

## How to cite
Please, cite this library using the following Publication:

Rincón and Amarloo (2023)

         @article{rincon2023progressive,
         title = {Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation},
         journal = {International Journal of Heat and Fluid Flow},
         volume = {104},
         pages = {109242},
         year = {2023},
         issn = {0142-727X},
         doi = {https://doi.org/10.1016/j.ijheatfluidflow.2023.109242},
         author = {Mario Javier Rincón and Ali Amarloo and Martino Reclari and Xiang I.A. Yang and Mahdi Abkar},
         keywords = {Turbulence modelling, RANS, Progressive augmentation, Surrogate modelling, Kriging, Secondary flows}
         }
                
For release-specific DOIs, click on the badge and find the DOI corresponding to the desired version in the version list.

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.

Detailed information on the OpenFOAM trademark can be found at

http://www.openfoam.com/legal/trademark-policy.php
http://www.openfoam.com/legal/trademark-guidelines.php
For further information on OpenCFD and OpenFOAM, please refer to

http://www.openfoam.com
