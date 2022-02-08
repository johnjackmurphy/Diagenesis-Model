# Diagenesis-model
## Numerical model on fluid-rock interaction in carbonate rocks (Ahm et al. 2018)
### Version JGM_v5
### Updated: 2021-02-08

Building on published modeling efforts (Banner & Hanson, 1990; Blattler & Higgins, 2015; Fantle & Higgins, 2014), we here present a MatLab based (R2016b) numerical model that simulates early marine diagenesis of carbonate sediments by stoichiometrically recrystallizing bulk carbonate sediment of mixed mineralogy (aragonite, high-Mg calcite, low-Mg calcite, dolomite) to low-Mg calcite (neomorphism) or to dolomite (dolomitization) assuming conservation of mass of carbon in the sediment. For a detailed description of the model setup we refer to Ahm et al. 2018 (GCA Vol. 236, pp 140-159, doi:10.1016/j.gca.2018.02.042).

To run the model you will need the following files in your MatLab path:

(1) **DiagModel_JGM_BASE_5.m**

This is the main script used to run the model. In this script you define the initial conditions such as primary mineralogy, elemental & isotopic composition of the primary minerals, diagenetic mineralogy (M), reaction rate (R), flow rate (u), length scale of flow path (number of boxes), and elemental & isotopic composition of the diagenetic fluid.  

Edit "MODEL VARIABLES" (~ lines 25 - 133) to set up your model run. 


(2) **DiagModel_JGM_constants_5.m**

A script that sets some basic constants e.g. volume of each box, sediment porosity, stoichiometry of the primary minerals, etc. 

You may want to adjust "Physical Features of Box" and Initial "Elemental Composition of Solid". 

(3) **DiagModel_JGM_dJdt_5.m**

The main ODE solver than calculates the change in mass and isotopic composition over time. You shouldn't need to change anything here. 


(4) **DiagModel_JGM_fluxes_5.m**

This script calculates the input and output fluxes into each box at every time step. These fluxes are then pulled into the ODE solver (DiagModel_JGM_dJdt_5.m) to solve dM/dt and dM_delta_/dt. Nothing to change here. 

============================================================================
*The following are optional: 
If you don't want to use them, comment out their function calls at the end of the BASE script (~lines 511 - 512)*


(5) **savediagmodelparaset_v5.m**

This saves the adjustable model parameter values to a spreadsheet so that you can review parameter sets for past runs. 

!! You must adjust the script to include a path to the directory on your computer where you want to save the spreadsheet !!


(6) **savediagmodelrun.m**

This saves the output of the current model run (your MatLab workspace) so you can pull up the results of past model runs without actually running the model again. 

!! You must adjust the script to include a path to the directory on your computer where you want to save the spreadsheet !!


(7) **odewbar.m**

This provides a model progress bar in a pop-up window (option can be turned off by commenting out line185 in DiagModel_JGM_BASE_5.m). 
