import math
import random
import time
from multiprocessing import Process, Queue, cpu_count
import numpy as np
import photon
import input_data
import manageEnv


#*** LAUNCH 
#   Initialize photon position and trajectory.
#   Implements an isotropic point source.
#****

def launch(queue_photons, environmentGeneral, queue_result, queue_end, number, envManager):
    while True:
        iphoton = queue_photons.get()
        if iphoton == "DONE":
            print("END"+str(number))
            queue_result.put("DONE")
            break
        mcflag = environmentGeneral["mcflag"]        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
        radius = environmentGeneral["radius"]        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
        tissueRefractive = environmentGeneral["tissueRefractive"]    #refractive index of tissue, former tissueRefractive
        tissueExtMedium = environmentGeneral["tissueExtMedium"]      #refractive index of external medium, former tissueExtMedium
        zfocus = environmentGeneral["zfocus"]        #depth of focus (if Gaussian (if mcflag = 1)
        PI = 3.1415926
        xs = environmentGeneral["xs"]                #used if mcflag = 2, isotropic pt source
        ys = environmentGeneral["ys"]                #used if mcflag = 2, isotropic pt source
        zs = environmentGeneral["zs"]                #used if mcflag = 2, isotropic pt source
        boundaryflag = environmentGeneral["boundaryflag"]    #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
        dr = environmentGeneral["dr"]                        #radial bin size [cm]
        dz = environmentGeneral["dz"]                        #depth bin size [cm]
        NR = environmentGeneral["NR"]                        #number of radial bins
        NZ = environmentGeneral["NZ"]                        #number of depth bins
        waist = environmentGeneral["waist"]                  #1/e radius of Gaussian focus
        THRESHOLD = environmentGeneral["THRESHOLD"]          #used in roulette        
        CHANCE = environmentGeneral["CHANCE"]                #used in roulette
    	
        envStart = envManager.envDefault.get_variables()
        mut = envStart["mua"] + envStart["mus"]
        albedo = envStart["mus"]/mut
        g = envStart["excitAnisotropy"]
        name = envStart["name"]

        absorbInfo = [] #matrix of fluence rate [cm^-2] = [W/cm2 per W], former F
        escapeFlux = []  #vector of escaping flux [cm^-2] versus radial position, former J
        tempRsptot = 0   
        Atot = 0

        phot = photon.photon(iphoton)

        if mcflag == 0:
            # UNIFORM COLLIMATED BEAM INCIDENT AT z = zs 
            # Launch at (r,zz) = (radius*sqrt(rnd), 0).
            # * Due to cylindrical symmetry, radial launch position is
            # * assigned to x while y = 0. 
            # * radius = radius of uniform beam. 
            # Initial position 
            rnd = float(random.random()) #random takes uniformly a number of (0,1]
            phot.update_positon(radius*math.sqrt(rnd), 0, zs)
            #x = radius*math.sqrt(rnd)
            #y = 0
            #z = zs
            # Initial trajectory as cosines 
            ux = 0
            uy = 0
            uz = 1  
            # specular reflectance 
            temp   = tissueRefractive/tissueExtMedium # refractive index mismatch, internal/external 
            temp   = (1.0 - temp)/(1.0 + temp)
            rsp    = temp*temp # specular reflectance at boundary 
            
        elif (mcflag == 1):
            # GAUSSIAN BEAM AT SURFACE 
            # Launch at (r,z) = (radius*sqrt(-log(rnd)), 0).
            # * Due to cylindrical symmetry, radial launch position is
            # * assigned to x while y = 0. 
            # * radius = 1/e radius of Gaussian beam at surface. 
            # * waist  = 1/e radius of Gaussian focus.
            # * zfocus = depth of focal point. 
            # Initial position 
            # avoids rnd = 0 
            rnd = float(random.random())
            x = radius*math.sqrt(-math.log(rnd))
            phot.update_positon(x, 0, 0)
            #y = 0.0
            #z = 0.0
            # Initial trajectory as cosines 
            # Due to cylindrical symmetry, radial launch trajectory is
            # * assigned to ux and uz while uy = 0. 
            # avoids rnd = 0  
            xfocus = waist*math.sqrt(-math.log(rnd))
            temp = math.sqrt((x - xfocus)*(x - xfocus) + zfocus*zfocus)
            sintheta = -(x - xfocus)/temp
            costheta = zfocus/temp
            ux = sintheta
            uy = 0.0
            uz = costheta
            # specular reflectance and refraction 
            (rsp, uz) = RFresnel(tissueExtMedium, tissueRefractive, costheta) # new uz 
            ux  = math.sqrt(1.0 - uz*uz) # new ux 
        elif  (mcflag == 2):
            # ISOTROPIC POINT SOURCE AT POSITION xs,ys,zs 
            # Initial position 
            rnd = float(random.random())
            phot.update_positon(xs, ys, zs)
            #x = xs
            #y = ys
            #z = zs
            # Initial trajectory as cosines 
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
            # specular reflectance 
            rsp = 0.0
        else: 
            print("choose mcflag between 0 to 2\n")
            


        phot.set_weight(1.0 - rsp)  # set photon initial weight 
        tempRsptot += rsp # accumulate specular reflectance per photon 

        #*****************************************
        #****** HOP_ESCAPE_SPINCYCLE **************
        #* Propagate one photon until it dies by ESCAPE or ROULETTE. 
        #******************************************

        while phot.get_status():
            #*** HOP
            # * Take step to new position
            # * s = stepsize
            # * ux, uy, uz are cosines of current photon trajectory
            # ****
            rnd = random.random()   # avoids rnd = 0 
            s = -math.log(rnd)/mut   # Step size.  Note: log() is base e 

            # Does photon ESCAPE at surface? ... z + s*uz <= 0? 
            if ((boundaryflag == 1) & (phot.get_position()[2] + s*uz <= 0)):
                rnd = random.random()
                # Check Fresnel reflectance at surface boundary 
                rf, uz1 = RFresnel(tissueRefractive, tissueExtMedium, -uz)
                if (rnd > rf): 
                    # Photon escapes at external angle, uz1 = cos(angle) 
                    s  = abs(phot.get_position()[2]/uz) # calculate stepsize to reach surface 
                    phot.update_positon(s*ux, s*uy, 0)
                    #x += s*ux       # partial step to reach surface 
                    #y += s*uy
                    pos = phot.get_position()
                    r = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])   # find radial position r 
                    ir = round((r/dr) + 1) # round to 1 <= ir 
                    ir = min(ir,NR)  # ir = NR is overflow bin 
                    escapeFlux.append((ir, phot.get_weight()))      # increment escaping flux 
                    phot.update_dead()
                    
                else:
                    pos = phot.get_position()
                    phot.update_positon(0, 0, -(pos[2] + s*uz))
                    #z = -(z + s*uz)   # Total internal reflection. 
                    uz = -uz
            else:
                photPosOld = phot.get_position()
                phot.update_positon(s*ux, s*uy, s*uz)   #update Positions
                photPosNew = phot.get_position()
                #check if boundary over new environment was crossed
                if not checkSameBin(photPosOld[0], photPosOld[1], photPosNew[0], photPosNew[1], dr, dr):
                    newEnv = envManager.findEnv(photPosNew[0], photPosNew[1], photPosNew[2], name)
                    if newEnv != 0:
                        mut = newEnv["mua"] + newEnv["mus"]
                        albedo = newEnv["mus"]/mut
                        g = newEnv["excitAnisotropy"]
                        name = newEnv["name"]
                #x += s*ux           # Update positions. 
                #y += s*uy
                #z += s*uz
                    

            if phot.get_status():
                #********************************************
                # ****** SPINCYCLE = DROP_SPIN_ROULETTE ******
                # ********************************************

                #*** DROP
                # * Drop photon weight (W) into local bin.
                # ****
                absorb = phot.get_weight()*(1 - albedo)       # photon weight absorbed at this step 
                phot.update_weight(-absorb)                  # decrement WEIGHT by amount absorbed 
                Atot += absorb               # accumulate absorbed photon weight 
                # deposit power in cylindrical coordinates z,r 
                pos = phot.get_position()
                r  = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])         # current cylindrical radial position 
                ir = round((r/dr) + 1)        
                iz = round((abs(pos[2])/dz) + 1)  
                ir = min(ir, NR)    # round to 1 <= ir, iz 
                iz = min(iz, NZ)    # last bin is for overflow 
                absorbInfo.append((iz*NR + ir, absorb))          # DROP absorbed weight into bin 

                #*** SPIN 
                # * Scatter photon into new trajectory defined by theta and psi.
                # * Theta is specified by cos(theta), which is determined 
                # * based on the Henyey-Greenstein scattering function.
                # * Convert theta and psi into cosines ux, uy, uz. 
                # ****
                # Sample for costheta 
                rnd = random.random()
                if (g == 0.0):
                    costheta = 2.0*rnd - 1.0
                elif (g == 1.0):
                    costheta = 1.0
                else:
                    temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
                    costheta = (1.0 + g*g - temp*temp)/(2.0*g)
                sintheta = math.sqrt(1.0 - costheta*costheta)	#sqrt faster than sin()

                # Sample psi. 
                psi = 2.0*PI*random.random()
                cospsi = math.cos(psi)
                if (psi < PI):
                    sinpsi = math.sqrt(1.0 - cospsi*cospsi)	#sqrt faster 
                else:
                    sinpsi = -math.sqrt(1.0 - cospsi*cospsi)

                # New trajectory. 
                if (1 - abs(uz) <= 1.0e-12): # close to perpendicular. 
                    uxx = sintheta*cospsi
                    uyy = sintheta*sinpsi
                    if uz>=0:
                        uzz = costheta
                    else:
                        uzz = -costheta
                else:  # usually use this option 
                    temp = math.sqrt(1.0 - uz*uz)
                    uxx = sintheta*(ux*uz*cospsi - uy*sinpsi)/temp + ux*costheta
                    uyy = sintheta*(uy*uz*cospsi + ux*sinpsi)/temp + uy*costheta
                    uzz = -sintheta*cospsi*temp + uz*costheta

                # Update trajectory 
                ux = uxx
                uy = uyy
                uz = uzz

                #*** CHECK ROULETTE 
                # * If photon weight below THRESHOLD, then terminate photon using
                # * Roulette technique. Photon has CHANCE probability of having 
                # * its weight increased by factor of 1/CHANCE,
                # * and 1-CHANCE probability of terminating.
                # ****
                if (phot.get_weight() < THRESHOLD):
                    rnd = random.random()
                    if (rnd <= CHANCE):
                        phot.set_weight(phot.get_weight() / CHANCE)
                        #W /= CHANCE
                    else:
                        phot.update_dead()
        

        #return absorbInfo, escapFlux, tempRspot, Atot
        # [[(1D-number,absorbValue), ...], [(number, weight), ...], value, value]
        queue_result.put([absorbInfo, escapeFlux, tempRsptot, Atot])

def launcher_start(processes, environmentGeneral, queue_result, queue_photons, queue_end, envManager):
    #*** LAUNCH 
    #   Initialize photon position and trajectory.
    #   Implements an isotropic point source.
    #***   
    all_launcher_procs = []
    for i in range(0, processes):
        print("launcher " + str(i)+ " started")
        launcher_p = Process(target=launch, args=((queue_photons), (environmentGeneral), (queue_result),(queue_end),(i), (envManager),))
        launcher_p.daemon = True
        launcher_p.start()

        all_launcher_procs.append(launcher_p)

    return all_launcher_procs

def writer(count, num_of_reader_procs, queue):
    """Write integers into the queue.  A reader_proc() will read them from the queue"""
    for ii in range(0, count):
        queue.put(ii)  # Put 'count' numbers into queue

    ### Tell all readers to stop...
    for ii in range(0, num_of_reader_procs+1):
        queue.put("DONE")

def sort(queue, processes, escapeFlux, absorbInfo, tempRsptot, Atot): 
    #launchReturn contains absorb Information, escaped photons, absorbed photon weight and specular reflectance

    seen_processes = 0
    while seen_processes < processes:
        launchReturn = queue.get()
        if launchReturn == "DONE":
            print("DONE %s" % seen_processes)
            seen_processes += 1
        else:
            for i in launchReturn[0]:
                absorbInfo[i[0]] += i[1]
            for j in launchReturn[1]:
                escapeFlux[j[0]] += j[1]
            tempRsptot += launchReturn[2]
            Atot += launchReturn[3]
    print("sort DONE")
    return [absorbInfo, escapeFlux, tempRsptot, Atot]


def checkSameBin( x1, y1,  x2,  y2,  dx, dy):
    return (x1/dx == x2/dx and y1/dy == y2/dy)
#*********************************************
#  **** END of SPINCYCLE = DROP_SPIN_ROULETTE *
#  *********************************************

def RFresnel(n1,		# incident refractive index.
            n2,		# transmit refractive index.
            ca1): 	# pointer to the cosine 
                    # of the transmission 
                    # angle a2, a2>0. 
    if n1==n2: #* matched boundary. *
        ca2_Ptr = ca1
        r = 0.0
    elif ca1>(1.0 - 1.0e-12): #* normal incidence. *
        ca2_Ptr = ca1
        r = (n2-n1)/(n2+n1)
        r *= r
    elif ca1< 1.0e-6:	#* very slanted. *
        ca2_Ptr = 0.0
        r = 1.0
    else:	  		#* general. *
        sa1 = sa2 = 0 # sine of incident and transmission angles. 
        ca2 = 0    # cosine of transmission angle. 
        sa1 = math.sqrt(1-ca1*ca1)
        sa2 = n1*sa1/n2
        if sa2>=1.0:
            #  check for total internal reflection. 
            ca2_Ptr = 0.0
            r = 1.0
        else:
            cap = cam = 0	# cosines of sum ap or diff am of the two 
                                # angles: ap = a1 + a2, am = a1 - a2. 
            sap = sam = 0	# sines. 
            ca2_Ptr = ca2 = math.sqrt(1-sa2*sa2)
            cap = ca1*ca2 - sa1*sa2 # c+ = cc - ss. 
            cam = ca1*ca2 + sa1*sa2 # c- = cc + ss. 
            sap = sa1*ca2 + ca1*sa2 # s+ = sc + cs. 
            sam = sa1*ca2 - ca1*sa2 # s- = sc - cs. 
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam) 
            # rearranged for speed. 
    return r, ca2_Ptr
#******* END SUBROUTINE *********

#**********************************************************
# * SAVE RESULTS TO FILES 
#**********************************************************
def SaveFile(Nfile,  J,  F,  S,  A,  E, environmentGeneral, Nphotons):

    mua = environmentGeneral["mua"]
    mus = environmentGeneral["mus"]
    mcflag = environmentGeneral["mcflag"]        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    radius = environmentGeneral["radius"]        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    tissueRefractive = environmentGeneral["tissueRefractive"]    #refractive index of tissue, former tissueRefractive
    tissueExtMedium = environmentGeneral["tissueExtMedium"]      #refractive index of external medium, former tissueExtMedium
    zfocus = environmentGeneral["zfocus"]        #depth of focus (if Gaussian (if mcflag = 1)
    PI = 3.1415926
    xs = environmentGeneral["xs"]                #used if mcflag = 2, isotropic pt source
    ys = environmentGeneral["ys"]                #used if mcflag = 2, isotropic pt source
    zs = environmentGeneral["zs"]                #used if mcflag = 2, isotropic pt source
    mut = environmentGeneral["mua"] + environmentGeneral["mus"]        #mua: absorption coefficient [cm^-1], mus: scattering coefficient [cm^-1]
    boundaryflag = environmentGeneral["boundaryflag"]    #boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    dr = environmentGeneral["dr"]                        #radial bin size [cm]
    dz = environmentGeneral["dz"]                        #depth bin size [cm]
    NR = environmentGeneral["NR"]                        #number of radial bins
    NZ = environmentGeneral["NZ"]                        #number of depth bins
    waist = environmentGeneral["waist"]                  #1/e radius of Gaussian focus
    g = environmentGeneral["excitAnisotropy"]                          #excitation anisotropy [dimensionless]
    THRESHOLD = environmentGeneral["THRESHOLD"]          #used in roulette        
    CHANCE = environmentGeneral["CHANCE"]                #used in roulette

    print("mcOUT%d.dat" % Nfile)
    file = open("mcOUT" + str(Nfile) + ".dat", "w")

    # print run parameters 
    file.write("%0.3e\tmua, absorption coefficient [1/cm]\n" % mua)
    file.write("%0.4f\tmus, scattering coefficient [1/cm]\n" % mus)
    file.write("%0.4f\tg, anisotropy [-]\n" % g)
    file.write("%0.4f\tn1, refractive index of tissue\n" % tissueRefractive)
    file.write("%0.4f\tn2, refractive index of outside medium\n" % tissueExtMedium)
    file.write("%d\tmcflag\n" % mcflag)
    file.write("%0.4f\tradius, radius of flat beam or 1/e radius of Gaussian beam [cm]\n" % radius)
    file.write("%0.4f\twaist, 1/e waist of focus [cm]\n" % waist)
    file.write("%0.4f\txs, x position of isotropic source [cm]\n" % xs)
    file.write("%0.4f\tys, y\n" % ys)
    file.write("%0.4f\tzs, z\n" % zs)
    file.write("%d\tNR\n" % NR)
    file.write("%d\tNZ\n" % NZ)
    file.write("%0.5f\tdr\n" % dr)
    file.write("%0.5f\tdz\n" % dz)
    file.write("%0.1e\tNphotons\n" % Nphotons)

    # print SAE values 
    file.write("%1.6e\tSpecular reflectance\n" % S)
    file.write("%1.6e\tAbsorbed fraction\n" % A)
    file.write("%1.6e\tEscaping fraction\n" % E)

    # print r[ir] to row 
    file.write("%0.1f" % 0.0) # ignore upperleft element of matrix 
    for ir in range(1, NR+1):
        r2 = dr*ir
        r1 = dr*(ir-1)
        r = 2.0/3*(r2*r2 + r2*r1 + r1*r1)/(r1 + r2)
        file.write("\t%1.5f" % r)
    file.write("\n")

    # print J[ir] to next row 
    file.write("%0.1f" % 0.0) # ignore this first element of 2nd row 
    for ir in range(1, NR+1):	
        file.write("\t%1.12e" % J[ir])
    file.write("\n")

    # printf z[iz], F[iz][ir] to remaining rows 
    for iz in range(1, NR+1):
        z = (iz - 0.5)*dz # z values for depth position in 1st column 
        file.write("%1.5f" % z)
        for ir in range(1, NR+1):
            file.write("\t %1.6e" % F[iz*NR + ir])
        file.write("\n")
    file.close()




if __name__ == '__main__':
    startTime = time.time()

    environmentGeneral = input_data.environmentGeneral
    envManager = manageEnv.manageEnv()

    for env in input_data.envDetail:
        envManager.addEnvironment(env["name"], env["mua"], env["mus"], env["excitAnisotropy"], env["formula"], env["default"])
    #**********************
    #* MAIN PROGRAM
    #*********************

    #************************
    #*** USER CHOICES *******
    #************************
    #number of file for saving
    Nfile = 0       # saves as mcOUTi.dat, where i = Nfile 
    # Run parameters 
    Nphotons = input_data.Nphotons # number of photons to be launched 
    PRINTOUT = 1
    #***********************
    #***** Setup output parameters, vectors, arrays *********
    #S     specular reflectance at air/tissue boundary 
    #A     total fraction of light absorbed by tissue 
    #E     total fraction of light escaping tissue 
    #J     escaping flux, J[ir], [W/cm2 per W incident] 
    #F     fluence rate, F[iz][ir], [W/cm2 per W incident] 
    escapeFlux       = np.zeros(environmentGeneral["NR"]+1)        # for escaping flux 
    absorbInfo       = np.zeros((environmentGeneral["NZ"]+1)*(environmentGeneral["NR"]+1)) # for absorbed fluence rate 
    #***********************
    #***********************

    #// choose Nphotons
    #THRESH = 1e-4
    #albedo = environmentGeneral["mus"]/(environmentGeneral["mua"] + environmentGeneral["mus"])
    #Nsteps = math.log(THRESH)/math.log(albedo)
    #tperstep = 249e-9
    #t = 30
    #Nphotons = round(t/(tperstep*Nsteps))
    print("Nphotons = %5.4e\n"% Nphotons)

    #**********************
    # DECLARE SUBROUTINES:
    #********************

    #*****
    #* The Monte Carlo subroutine mcsub()
    #*	Tissue properties:
    #*			mua       = absorption coefficient [cm^-1]
    #*			mus       = scattering coefficient [cm^-1]
    #*			g         = anisotropy of scattering [dimensionless]
    #*			n1        = refractive index of tissue
    #*			n2        = refractive index of external medium
    #*	Incident beam characteristics:
    #*			mcflag:   0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
    #*			xs,ys,zs  = position of istropic point source (if mcflag = 2)
    #*			boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
    #*			radius    = radius of beam (1/e width if Gaussian) (if mcflag < 2)
    #*			waist     = 1/e radius of focus (if Gaussian (if mcflag = 1)
    #*			zfocus    = depth of focus (if Gaussian (if mcflag = 1)
    #*	OUTPUT:
    #*			J[ir]     = vector of escaping flux [cm^-2] versus radial position
    #*			F[iz][ir] = matrix of fluence rate [cm^-2] = [W/cm2 per W]
    #*			Sptr      = pointer to S, specular refelctance
    #*			Aptr      = pointer to A, total absorbed fraction
    #*			Eptr      = pointer to E, total escaping fraction 
    #*

    #*** INITIALIZATIONS ****
    tempRsptot = 0.0 # accumulate specular reflectance per photon 
    Atot   = 0.0 # accumulate absorbed photon weight 
    processes = cpu_count()
    print(processes)
    queue_result = Queue()
    queue_photons = Queue()
    queue_end = Queue()
    #============================================================
    #======================= RUN N photons =====================
    # * Launch N photons, initializing each one before progation.
    #============================================================

    launcher_procs = launcher_start(processes, environmentGeneral, queue_result, queue_photons, queue_end, envManager)
    writer(Nphotons, processes, queue_photons)
    return_dict = sort(queue_result, processes, escapeFlux, absorbInfo, tempRsptot, Atot)
    for idx, a_launcher_proc in enumerate(launcher_procs):
        print("    Waiting for reader_p.join() index %s" % idx)

        a_launcher_proc.join()  # Wait for a_launcher_proc() to finish

        print("        reader_p() idx:%s is done" % idx)



    absorbInfo = return_dict[0]
    escapeFlux = return_dict[1]
    tempRsptot = return_dict[2]
    Atot = return_dict[3]


    #***********************
    # * NORMALIZE 
    # *   J[ir]      escaping flux density [W/cm^2 per W incident] 
    # *              where bin = 2.0*PI*r[ir]*dr [cm^2].
    # *	  F[iz][ir]  fluence rate [W/cm^2 per W incident] 
    # *              where bin = 2.0*PI*r[ir]*dr*dz [cm^3].
    # ***********************
    temp = 0.0
    NR = environmentGeneral["NR"]
    NZ = environmentGeneral["NZ"]
    dr = environmentGeneral["dr"]
    PI = environmentGeneral["PI"]
    dz = environmentGeneral["dz"]
    mua = environmentGeneral["mua"]

    for ir in range(1, NR+1):
        r = (ir - 0.5)*dr
        temp += escapeFlux[ir]    # accumulate total escaped photon weight 
        escapeFlux[ir] /= 2.0*PI*r*dr*Nphotons
        for iz in range(1, NZ+1):	                # flux density 
            absorbInfo[iz*NR + ir] /= 2.0*PI*r*dr*dz*Nphotons*mua # fluence rate 

    Sptr = S = tempRsptot/Nphotons
    Aptr = A = Atot/Nphotons
    Eptr = E = temp/Nphotons

    print("Nphotons = %5.1e\n" % Nphotons)
    print("Specular = %5.6f\n" % S)
    print("Absorbed = %5.6f\n" % A)
    print("Escaped  = %5.6f\n" % E)
    SAE= S+A+E
    print("total    = %5.6f\n" % SAE)
    SaveFile(1, escapeFlux, absorbInfo, S, A, E, environmentGeneral, Nphotons)
    endTime = time.time()
    print(endTime - startTime)
