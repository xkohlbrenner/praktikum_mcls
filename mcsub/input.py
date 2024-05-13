default = 0        #if set to 1 the default numbers are taken, 
                    #use this, if you're provided data throws Error

nt  = 1.33          #tissue index of refraction
Nphotons = 100000    #number of photons

environment = {
    "mcflag": mcflag,            	        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    "radius": radius,                        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    "tissueRefractive": tissueRefractive,    #refractive index of tissue, former n1
    "tissueExtMedium": tissueExtMedium,      #refractive index of external medium, former n2
    "zfocus": zfocus,                        #depth of focus (if Gaussian (if mcflag = 1)
    "PI": 3.1415926,
    "xs": xs,                               #used if mcflag = 2, isotropic pt source
    "ys": ys,                               #used if mcflag = 2, isotropic pt source
    "zs": zs,                               #used if mcflag = 2, isotropic pt source
    "mua": mua,                              #mua: absorption coefficient [cm^-1],
    "mus": mus,                              #mus: scattering coefficient [cm^-1]
    "boundaryflag": boundaryflag,           #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    "dr": radialBinSize,                    #radial bin size [cm]
    "dz": depthBinSize,                     #depth bin size [cm]
    "NR": radialBins,                       #number of radial bins
    "NZ": depthBins,                        #number of depth bins
    "waist": waist,                         #1/e radius of Gaussian focus
    "excitAnisotropy": g,                    #excitation anisotropy [dimensionless]
    "THRESHOLD": THRESHOLD,                  #used in roulette        
    "CHANCE": CHANCE,                        #used in roulette
}