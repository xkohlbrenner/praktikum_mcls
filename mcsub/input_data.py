# This is the input file. 
# If you are unsure about your changes, I highly recommend to make a copy of the initial file.



#number of photons
Nphotons = 100000    

#all data about the tissue and the beams
environmentGeneral = {
#0 = collimated uniform, 1 = Gaussian, 2 = isotropic point    
    "mcflag": 0,            	        

#variables for uniform and gaussian beams (mcflag= 1, 2)
    "radius": 0,                        #radius of beam (1/e width if Gaussian)

#gaussian beam variables (mcflag=1)
    "waist": 0.1,                        #1/e radius of Gaussian focus
    "zfocus": 1,                        #depth of focus

#isotropic point source variables (mcflag=2)   
    "xs": 0,                               
    "ys": 0,                               
    "zs": 0,                               

#tissue variables
    "radialSize": 0.003,                #radial size in [cm]
    "depthSize": 0.003,                 #depth size in [cm]
    "bins": 128,                         #number of radial and depth bins, needs to be 2^n
    "tissueExtMedium": 1,               #refractive index of external medium

#roulette    
    "THRESHOLD": 0.0001,                        
    "CHANCE": 0.1,                      

#constant variables
    "PI": 3.1415926,
}

#Here all environment data, inclusive the position at the tissue is declared
envDetail = [
    {
        "name": "standard tissue",
        "mua": 1,                               #mua: absorption coefficient [cm^-1],
        "mus": 100,                             #mus: scattering coefficient [cm^-1]
        "excitAnisotropy": 0.9,                 #excitation anisotropy [dimensionless]
        "refractiveIndex": 1.4,                 #refractive index
        "height": [0.001, 0.002],               #only needed, if default is 0
        "radius": 0.003,
        "default": 1
    },
    {
        "name": "blood",
        "mua": 230.5427,                              #mua: absorption coefficient [cm^-1],
        "mus": 93.9850,                              #mus: scattering coefficient [cm^-1]
        "excitAnisotropy": 0.9000,                    #excitation anisotropy [dimensionless]
        "refractiveIndex": 1.351,
        "height": [0.001, 0.002],                           #in [cm]
        "radius": 0.003,                                #in [cm]
        "default": 0
    }
]