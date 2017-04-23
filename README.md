# DMSensitivity
# Detector response needs to be modified to poisson distribution
# wimp-xs needs to be modified to sigma0\*F^{2}
This is the program used to calculate the sensitivity of a direct dark matter detection experiment.  

Usage  
You need to configure the codes first with the help of Makefile: make  

Description  
DMCalc.h: head file for sensitivity calculation  
DMCalc.C: source file contains functions which are necessary to calculate sensitivity including cross section, form factor, event rate, detector response  
DMCalc::DMCalc(const char\* material, bool IsSI, int mode)  
material= "Ar": material is Ar40;  
material= "Xe": material is Xe132;  
mode=1: vearth= 0, Vescape= inf;  
mode=2: vearth= 0, Vescape= Vesc;  
mode=3: vearth= Ve, Vescape= inf;  
mode=4: vearth= Ve, Vescape= Vesc;  

calc.C: main function to do calculation  


Reference  
J.D. Lewin, P.F. Smith. Astroparticle Physics 6 (1996) 87-112  
arXiv:1502.0266v3 [hep-ph]  
