import random
import math









pi = 3.1415926
lightspeed = 2.997925*pow(10,10)    #in vacuo in [cm/s]
alive = 1           #photon not yet terminated
dead = 0            #photon is terminated
threshold = 0.01
chance = 0.1
cos90d = 1.0*pow(0.1,6)     #If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad.
one_minus_coszero = 1.0*pow(0.1,12)

def sign(x): 
    return 1 if x>=0 else -1

initRandomGen = float(1-random.random()) #random takes uniformly a number of [0,1), but 1 is needed and zero not
randomNum = float(1-random.random())


if __name__ == '__main__':

    x, y, z= 0     #photon position
    ux, uy, uz  = 0 #photon trajectory as cosines
    uxx, uyy, uzz  = 0	#temporary values used during SPIN 
    s  = 0          #step sizes. s = -log(RND)/mus [cm] 
    costheta  = 0   #cos(theta) 
    sintheta  = 0   #sin(theta) 
    cospsi  = 0     #cos(psi) 
    sinpsi  = 0     #sin(psi) 
    psi  = 0        #azimuthal angle 
    i_photon  = 0   #current photon 
    W  = 0          #photon weight 
    absorb  = 0     #weighted deposited in a step due to absorption 
    photon_status  = 0  #flag = ALIVE=1 or DEAD=0 

    #other variables 
    Csph = [0]*101  #spherical   photon concentration CC[ir=0..100] 
    Ccyl = [0]*101  #cylindrical photon concentration CC[ir=0..100] 
    Cpla = [0]*101  #planar      photon concentration CC[ir=0..100] 
    Fsph  = 0       #fluence in spherical shell 
    Fcyl  = 0       #fluence in cylindrical shell 
    Fpla  = 0       #fluence in planar shell 
    mua  = 0        #absorption coefficient [cm^-1] 
    mus  = 0        #scattering coefficient [cm^-1] 
    g  = 0          #anisotropy [-] 
    albedo  = 0     #albedo of tissue 
    nt  = 0         #tissue index of refraction 
    Nphotons  = 0   #number of photons in simulation 
    NR  = 0         #number of radial positions 
    radial_size  = 0  #maximum radial size 
    r  = 0          #radial position 
    dr  = 0         #radial bin size 
    ir  = 0         #index to radial position 
    shellvolume  = 0  #volume of shell at radial position r 

    #dummy variables 
    rnd  = 0        #assigned random value 0-1 
    temp  = 0    #dummy variables 
    FILE*	target  = 0     #point to output file 


#INPUT

    mua         = 1.0     #cm^-1 
    mus         = 0.0     #cm^-1 
    g           = 0.90  
    nt          = 1.33
    Nphotons    = 10000 #set number of photons in simulation 
    radial_size = 3.0   #cm, total range over which bins extend 
    NR          = 100	 #set number of bins.  
    #IF NR IS ALTERED, THEN USER MUST ALSO ALTER THE ARRAY DECLARATION TO A SIZE = NR + 1. 
    dr          = radial_size/NR  #cm 
    albedo      = mus/(mus + mua)


    while i_photon < Nphotons:
        #LAUNCh
        #Initialize photon position and trajectory.
        #Implements an isotropic point source.
        i_photon +=1
        W = 1.0
        photon_status = 1

        x=0
        y=0
        z=0

        #Randomly set photon trajectory to yield an isotropic source.
        costheta = 2.0*(1-random.random()) - 1.0   
        sintheta = pow(1.0 - costheta*costheta, 2)	#sintheta is always positive
        psi = 2.0*pi*random.random()
        ux = sintheta*math.cos(psi)
        uy = sintheta*math.sin(psi)
        uz = costheta
    
        while photon_status:
            rnd = 1 - random.random()
            s = -math.log(rnd)/(mua + mus)
            x += s*ux
            y += s*uy
            z += s*uz

            #Drop photon weight (W) into local bin.
            absorb = W*(1 - albedo)      #photon weight absorbed at this step 
            W -= absorb                  #decrement WEIGHT by amount absorbed 
            
            #spherical 
            r = pow(x*x + y*y + z*z,2)   #current spherical radial position 
            ir = r/dr                    #ir = index to spatial bin 
            if (ir >= NR): ir = NR        #last bin is for overflow 
            Csph[ir] += absorb           #DROP absorbed weight into bin 
            
            #cylindrical 
            r = pow(x*x + y*y,2)         #current cylindrical radial position 
            ir = r/dr                    #ir = index to spatial bin 
            if (ir >= NR): ir = NR        #last bin is for overflow 
            Ccyl[ir] += absorb           #DROP absorbed weight into bin 
            
            #planar 
            r = abs(z)                  #current planar radial position 
            ir = r/dr                   #ir = index to spatial bin 
            if (ir >= NR): ir = NR        #last bin is for overflow 
            Cpla[ir] += absorb           #DROP absorbed weight into bin 

            #SPIN
            #Scatter photon into new trajectory defined by theta and psi.
            #Theta is specified by cos(theta), which is determined 
            #based on the Henyey-Greenstein scattering function.
            #Convert theta and psi into cosines ux, uy, uz. 

            #costheta sample
            rnd = 1-random.random()
            temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
            costheta = (1.0 + g*g - temp*temp)/(2.0*g)
            sintheta = pow(1.0 - costheta*costheta,2)

            #psi sample
            psi = 2.0*pi*random.Random()
            cospsi = math.cos(psi)
            if (psi < pi):
                sinpsi = pow(1.0 - cospsi*cospsi,2)
            else:
                sinpsi = -pow(1.0 - cospsi*cospsi,2)

            #new trajectory
            if (1 - abs(uz) <= one_minus_coszero):    #close to perpendicular. 
                uxx = sintheta * cospsi
                uyy = sintheta * sinpsi
                uzz = costheta * sign(uz)   #SIGN() is faster than division.  
            else:					#usually use this option 
                temp = pow(1.0 - uz * uz,2)
                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
                uzz = -sintheta * cospsi * temp + uz * costheta

            #Update trajectory
            ux = uxx
            uy = uyy
            uz = uzz

            #Check Roulette
            #If photon weight below THRESHOLD, then terminate photon using Roulette technique.
            #Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
            #and 1-CHANCE probability of terminating.
            if (W < threshold):
                if (random.random <= chance):
                    W /= chance
                else: 
                    photon_status = 0


    fprintf(target, "number of photons = %f\n", Nphotons);
    fprintf(target, "bin size = %5.5f [cm] \n", dr);
    fprintf(target, "last row is overflow. Ignore.\n");

    #print column titles
    fprintf(target, "r [cm] \t Fsph [1/cm2] \t Fcyl [1/cm2] \t Fpla [1/cm2]\n")

    #print data:  radial position, fluence rates for 3D, 2D, 1D geometries
    for ir in range(0,NR):
        #r = sqrt(1.0/3 - (ir+1) + (ir+1)*(ir+1))*dr;
        r = (ir + 0.5)*dr;
        shellvolume = 4.0*pi*r*r*dr; #per spherical shell
        Fsph = Csph[ir]/Nphotons/shellvolume/mua;
        shellvolume = 2.0*pi*r*dr;   #per cm length of cylinder
        Fcyl = Ccyl[ir]/Nphotons/shellvolume/mua;
        shellvolume = dr;            #per cm2 area of plane
        Fpla =Cpla[ir]/Nphotons/shellvolume/mua;
        fprintf(target, "%5.5f \t %4.3e \t %4.3e \t %4.3e \n", r, Fsph, Fcyl, Fpla);

    fclose(target);

