Folder 'triple-porosity' works as a module that must be inserted inside SINTEF produced software: mrst-202Xx/modules. 
The recent version of mrst-202Xx is avaiable at https://www.sintef.no/projectweb/mrst/.
Our code is developed based on the SINTEF modules, especially 'dual-continuum-mech' module written by Ashworth and Doster.
We decided to keep certain genearlisations of 'dual-continuum-mech' fluid model that are not discussed in Adamus et al (2023) 
for possible extended use of the triple-porosity simulations. 
For instance, our module supports extensions of the multi-phase flow, fluid compressibility factor, fluid well injection, and gravity force. 
The other parts of our code, like poroelastic parameterisation or fluid interflow follow directly the derivations in Adamus et al (2023). 
Apart from performing simulations to obtain results from Figure 5-6, user can verify time-dependent variables not illustrated 
in the paper, such as: local pressurization, fluid contant change, or volumetric strain.


Obtaining results from Figures 5-6 that make part of the numerical section in Adamus et al (2023):

0. Download MRST software and insert folder triple-porosity inside mrst-202Xx/modules
1. Open the script RUN.m in Matlab.
2. Run ''startup'' in the command line to activate MRST. 
3. Choose the relevant parameters listed in Table 2 and described in the numerical section.
4. Run the script.


Example: obtaining results from Figure 5a.

3a. Choose:  x=0.8 ,  y=0.5 ,  z=10 . 
3b. Keep other parameters unchanged.


References:
Adamus, F. P., Healy, D., Meredith, P. G., Mitchell, T. M., Stanton-Yonge, A. (2023). Multi-porous extension of anisotropic poroelasticity: consolidation and related coefficients. JGR: Solid Earth.