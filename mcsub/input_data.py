default = 0        #if set to 1 the default numbers are taken, 
                    #use this, if you're provided data throws Error

nt  = 1.33          #tissue index of refraction
Nphotons = 100000    #number of photons

environment = {
    "mcflag": 0,            	        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    "radius": 0,                        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    "tissueRefractive": 1.4,    #refractive index of tissue, former n1
    "tissueExtMedium": 1,      #refractive index of external medium, former n2
    "zfocus": 1,                        #depth of focus (if Gaussian (if mcflag = 1)
    "PI": 3.1415926,
    "xs": 0,                               #used if mcflag = 2, isotropic pt source
    "ys": 0,                               #used if mcflag = 2, isotropic pt source
    "zs": 0,                               #used if mcflag = 2, isotropic pt source
    "mua": 1,                              #mua: absorption coefficient [cm^-1],
    "mus": 100,                              #mus: scattering coefficient [cm^-1]
    "boundaryflag": 1,           #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    "dr": 0.002,                    #radial bin size [cm]
    "dz": 0.002,                     #depth bin size [cm]
    "NR": 100,                       #number of radial bins
    "NZ": 100,                        #number of depth bins
    "waist": 0.1,                         #1/e radius of Gaussian focus
    "excitAnisotropy": 0.9,                    #excitation anisotropy [dimensionless]
    "THRESHOLD": 0.0001,                  #used in roulette        
    "CHANCE": 0.1,                        #used in roulette
}