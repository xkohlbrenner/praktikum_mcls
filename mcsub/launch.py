import math
import random
import photon as photon

    #/**** LAUNCH 
    #   Initialize photon position and trajectory.
    #   Implements an isotropic point source.
    #*****/

def launch(iphoton, environment):
    mcflag = environment["mcflag"]        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    radius = environment["radius"]        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    tissueRefractive = environment["tissueRefractive"]    #refractive index of tissue, former tissueRefractive
    tissueExtMedium = environment["tissueExtMedium"]      #refractive index of external medium, former tissueExtMedium
    zfocus = environment["zfocus"]        #depth of focus (if Gaussian (if mcflag = 1)
    PI = 3.1415926
    xs = environment["xs"]                #used if mcflag = 2, isotropic pt source
    ys = environment["ys"]                #used if mcflag = 2, isotropic pt source
    zs = environment["zs"]                #used if mcflag = 2, isotropic pt source
    mut = environment["mua"] + environment["mus"]        #mua: absorption coefficient [cm^-1], mus: scattering coefficient [cm^-1]
    boundaryflag = environment["boundaryflag"]    #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    dr = environment["radialBinSize"]                        #radial bin size [cm]
    dz = environment["depthBinSize"]                        #depth bin size [cm]
    NR = environment["radialBins"]                        #number of radial bins
    NZ = environment["depthBins"]                        #number of depth bins
    waist = environment["waist"]                  #1/e radius of Gaussian focus
    g = environment["g"]                          #excitation anisotropy [dimensionless]
    THRESHOLD = environment["THRESHOLD"]          #used in roulette        
    CHANCE = environment["CHANCE"]                #used in roulette
    albedo = environment["mus"]/environment["mut"]

    absorbInfo = [] #matrix of fluence rate [cm^-2] = [W/cm2 per W], former F
    escapeFlux = []  #vector of escaping flux [cm^-2] versus radial position, former J
    tempRsptot = 0   
    Atot = 0

    phot = photon.photon(iphoton)

    if mcflag == 0:
        #/* UNIFORM COLLIMATED BEAM INCIDENT AT z = zs */
        #/* Launch at (r,zz) = (radius*sqrt(rnd), 0).
        # * Due to cylindrical symmetry, radial launch position is
        # * assigned to x while y = 0. 
        # * radius = radius of uniform beam. */
        #/* Initial position */
        rnd = float(random.random()) #random takes uniformly a number of (0,1]
        phot.update_positon(radius*math.sqrt(rnd), 0, zs)
        #x = radius*math.sqrt(rnd)
        #y = 0
        #z = zs
        #/* Initial trajectory as cosines */
        ux = 0
        uy = 0
        uz = 1  
        #/* specular reflectance */
        temp   = tissueRefractive/tissueExtMedium #/* refractive index mismatch, internal/external */
        temp   = (1.0 - temp)/(1.0 + temp)
        rsp    = temp*temp #/* specular reflectance at boundary */
        
    elif (mcflag == 1):
        #/* GAUSSIAN BEAM AT SURFACE */
        #/* Launch at (r,z) = (radius*sqrt(-log(rnd)), 0).
        # * Due to cylindrical symmetry, radial launch position is
        # * assigned to x while y = 0. 
        # * radius = 1/e radius of Gaussian beam at surface. 
        # * waist  = 1/e radius of Gaussian focus.
        # * zfocus = depth of focal point. */
        #/* Initial position */
        #/* avoids rnd = 0 */
        rnd = float(random.random())
        x = radius*math.sqrt(-math.log(rnd))
        phot.update_positon(x, 0, 0)
        #y = 0.0
        #z = 0.0
        #/* Initial trajectory as cosines */
        #/* Due to cylindrical symmetry, radial launch trajectory is
        # * assigned to ux and uz while uy = 0. */
        #/* avoids rnd = 0 */ 
        xfocus = waist*math.sqrt(-math.log(rnd))
        temp = math.sqrt((x - xfocus)*(x - xfocus) + zfocus*zfocus)
        sintheta = -(x - xfocus)/temp
        costheta = zfocus/temp
        ux = sintheta
        uy = 0.0
        uz = costheta
        #/* specular reflectance and refraction */
        (rsp, uz) = RFresnel(tissueExtMedium, tissueRefractive, costheta) #/* new uz */
        ux  = math.sqrt(1.0 - uz*uz) #/* new ux */
    elif  (mcflag == 2):
        #/* ISOTROPIC POINT SOURCE AT POSITION xs,ys,zs */
        #/* Initial position */
        rnd = float(random.random())
        phot.update_positon(xs, ys, zs)
        #x = xs
        #y = ys
        #z = zs
        #/* Initial trajectory as cosines */
        costheta = 1.0 - 2.0*random.random()
        sintheta = math.sqrt(1.0 - costheta*costheta)
        psi = 2.0*PI*random.random()
        cospsi = math.cos(psi)
        if (psi < PI):
            sinpsi = math.sqrt(1.0 - cospsi*cospsi) 
        else:
            sinpsi = -math.sqrt(1.0 - cospsi*cospsi)
        ux = sintheta*cospsi
        uy = sintheta*sinpsi
        uz = costheta
        #/* specular reflectance */
        rsp = 0.0
    else: 
        print("choose mcflag between 0 to 2\n")
        

    phot.set_weight(1.0 - rsp)  #/* set photon initial weight */
    tempRsptot += rsp #/* accumulate specular reflectance per photon */

    #/******************************************
    #****** HOP_ESCAPE_SPINCYCLE **************
    #* Propagate one photon until it dies by ESCAPE or ROULETTE. 
    #*******************************************/

    while phot.get_status():
        #/**** HOP
        # * Take step to new position
        # * s = stepsize
        # * ux, uy, uz are cosines of current photon trajectory
        # *****/
        rnd = random.random()   #/* avoids rnd = 0 */
        s = -math.log(rnd)/mut   #/* Step size.  Note: log() is base e */

        #/* Does photon ESCAPE at surface? ... z + s*uz <= 0? */
        if ((boundaryflag == 1) & (phot.get_position()[2] + s*uz <= 0)):
            rnd = random.random()
            #/* Check Fresnel reflectance at surface boundary */
            rf, uz1 = RFresnel(tissueRefractive, tissueExtMedium, -uz)
            if (rnd > rf): 
                #/* Photon escapes at external angle, uz1 = cos(angle) */
                s  = abs(phot.get_position()[2]/uz) #/* calculate stepsize to reach surface */
                phot.update_positon(s*ux, s*uy, 0)
                #x += s*ux       #/* partial step to reach surface */
                #y += s*uy
                pos = phot.get_position()
                r = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])   #/* find radial position r */
                ir = round((r/dr) + 1) #/* round to 1 <= ir */
                ir = min(ir,NR)  #/* ir = NR is overflow bin */
                escapeFlux.append((ir, phot.get_weight()))      #/* increment escaping flux */
                phot.update_dead()
                
            else:
                pos = phot.get_position()
                phot.update_positon(0, 0, -(pos[2] + s*uz))
                #z = -(z + s*uz)   #/* Total internal reflection. */
                uz = -uz
        else:
            phot.update_positon(s*ux, s*uy, s*uz)
            #x += s*ux           #/* Update positions. */
            #y += s*uy
            #z += s*uz
                

        if phot.get_status():
            #/*********************************************
            # ****** SPINCYCLE = DROP_SPIN_ROULETTE ******
            # *********************************************/

            #/**** DROP
            # * Drop photon weight (W) into local bin.
            # *****/
            absorb = phot.get_weight()*(1 - albedo)       #/* photon weight absorbed at this step */
            phot.update_weight(-absorb)                  #/* decrement WEIGHT by amount absorbed */
            Atot += absorb               #/* accumulate absorbed photon weight */
            #/* deposit power in cylindrical coordinates z,r */
            pos = phot.get_position()
            r  = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])         #/* current cylindrical radial position */
            ir = round((r/dr) + 1)        
            iz = round((abs(pos[2])/dz) + 1)  
            ir = min(ir, NR)    #/* round to 1 <= ir, iz */
            iz = min(iz, NZ)    #/* last bin is for overflow */
            absorbInfo.append((iz*NR + ir, absorb))          #/* DROP absorbed weight into bin */

            #/**** SPIN 
            # * Scatter photon into new trajectory defined by theta and psi.
            # * Theta is specified by cos(theta), which is determined 
            # * based on the Henyey-Greenstein scattering function.
            # * Convert theta and psi into cosines ux, uy, uz. 
            # *****/
            #/* Sample for costheta */
            rnd = random.random()
            if (g == 0.0):
                costheta = 2.0*rnd - 1.0
            elif (g == 1.0):
                costheta = 1.0
            else:
                temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
                costheta = (1.0 + g*g - temp*temp)/(2.0*g)
            sintheta = math.sqrt(1.0 - costheta*costheta)	#/*sqrt faster than sin()*/

            #/* Sample psi. */
            psi = 2.0*PI*random.random()
            cospsi = math.cos(psi)
            if (psi < PI):
                sinpsi = math.sqrt(1.0 - cospsi*cospsi)	#/*sqrt faster */
            else:
                sinpsi = -math.sqrt(1.0 - cospsi*cospsi)

            #/* New trajectory. */
            if (1 - abs(uz) <= 1.0e-12): #/* close to perpendicular. */
                uxx = sintheta*cospsi
                uyy = sintheta*sinpsi
                if uz>=0:
                    uzz = costheta
                else:
                    uzz = -costheta
            else:  #/* usually use this option */
                temp = math.sqrt(1.0 - uz*uz)
                uxx = sintheta*(ux*uz*cospsi - uy*sinpsi)/temp + ux*costheta
                uyy = sintheta*(uy*uz*cospsi + ux*sinpsi)/temp + uy*costheta
                uzz = -sintheta*cospsi*temp + uz*costheta

            #/* Update trajectory */
            ux = uxx
            uy = uyy
            uz = uzz

            #/**** CHECK ROULETTE 
            # * If photon weight below THRESHOLD, then terminate photon using
            # * Roulette technique. Photon has CHANCE probability of having 
            # * its weight increased by factor of 1/CHANCE,
            # * and 1-CHANCE probability of terminating.
            # *****/
            if (phot.get_weight() < THRESHOLD):
                rnd = random.random()
                if (rnd <= CHANCE):
                    phot.set_weight(phot.get_weight() / CHANCE)
                    #W /= CHANCE
                else:
                    phot.update_dead()
    

    #return absorbInfo, escapFlux, tempRspot, Atot
    # [[(1D-number,absorbValue), ...], [(number, weight), ...], value, value]
    return  [absorbInfo, escapeFlux, tempRsptot, Atot]

#/**********************************************
#  **** END of SPINCYCLE = DROP_SPIN_ROULETTE *
#  **********************************************/
