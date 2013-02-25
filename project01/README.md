Run the program using no arguments (it will print the default settings and use them), or use the following 10 arguments (in that order):
 - timestep dt (in MD units)
 - number of steps nSteps
 - initial "temperature" T (MD units)
 - thermostat (0 = no thermostat, 1 = berendsen, 2 = andersen)
 - temperature of heat bath Tbath (MD units)
 - tau_ for thermostat (tau = dt*tau_)
 - length of unit cell L (SI units)
 - number of unit cells N_C
 - bool calculate statistics? (1 = yes, 0 = no) - this is needed for the thermostat to work
 - bool save states in xyz-file each timestepd? (1 = yes, 0 = no)
 
