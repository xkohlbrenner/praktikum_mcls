Nphotons = 100000    #number of photons

environmentGeneral = {
    "mcflag": 0,            	        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    "radius": 0,                        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    "tissueExtMedium": 1,               #refractive index of external medium
    "zfocus": 1,                        #depth of focus (if Gaussian (if mcflag = 1)
    "PI": 3.1415926,
    "xs": 0,                               #isotropic pt source, used if mcflag = 2
    "ys": 0,                               #isotropic pt source, used if mcflag = 2
    "zs": 0,                               #isotropic pt source, used if mcflag = 2
    "boundaryflag": 1,           #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    "radialSize": 0.003,                #radial size in [cm]
    "depthSize": 0.003,                 #depth size in [cm]
    "NR": 256,                       #number of radial bins
    "NZ": 256,                        #number of depth bins
    "bins": 64,
    "waist": 0.1,                         #1/e radius of Gaussian focus
    "THRESHOLD": 0.0001,                  #used in roulette        
    "CHANCE": 0.1,                        #used in roulette
    "mua": 1,                              #mua: absorption coefficient [cm^-1],
    "mus": 100,                              #mus: scattering coefficient [cm^-1]
    "excitAnisotropy": 0.9,
}

envDetail = [
    {
        "name": "standard tissue",
        "mua": 1,                              #mua: absorption coefficient [cm^-1],
        "mus": 100,                              #mus: scattering coefficient [cm^-1]
        "excitAnisotropy": 0.9,                    #excitation anisotropy [dimensionless]
        "refractiveIndex": 1.4,
        "height": [0.001, 0.002],                           #only needed, if default is 0
        "radius": 0.001,
        "default": 1
    },
    {
        "name": "blood",
        "mua": 230.5427,                              #mua: absorption coefficient [cm^-1],
        "mus": 93.9850,                              #mus: scattering coefficient [cm^-1]
        "excitAnisotropy": 0.9000,                    #excitation anisotropy [dimensionless]
        "refractiveIndex": 1.351,
        "height": [0.001, 0.002],                           #in [cm]
        "radius": 0.001,                                #in [cm]
        "default": 0
    }
]