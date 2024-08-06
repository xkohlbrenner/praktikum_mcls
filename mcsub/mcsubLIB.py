import math
import random
import time
from multiprocessing import Process, Queue, Event, Manager, cpu_count
import numpy as np
import photon
import input_data
import manageEnv
from quadTree import QuadTree


#*** LAUNCH 
# Initialize photon position and trajectory.
#***

def launch(queue_photons, environmentGeneral, queue_result, ID, tree):
    ls = 0.1/(10*7)                                 #little step: to move the photon a bit in the voxel face
    while True:                                     #the loop get stoped, when the photon queue is empty
        iphoton = queue_photons.get()
        if iphoton == "DONE":                       #stops the launcher
            print("END "+str(ID))
            queue_result.put("DONE")
            break
        mcflag = environmentGeneral["mcflag"]        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic po
        radius = environmentGeneral["radius"]        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
        tissueExtMedium = environmentGeneral["tissueExtMedium"]      #refractive index of external medium, former tissueExtMedium
        zfocus = environmentGeneral["zfocus"]                       #depth of focus (if Gaussian (if mcflag = 1)
        PI = 3.1415926
        xs = environmentGeneral["xs"]                #used if mcflag = 2, isotropic pt source
        ys = environmentGeneral["ys"]                #used if mcflag = 2, isotropic pt source
        zs = environmentGeneral["zs"]                #used if mcflag = 2, isotropic pt source
        dr = environmentGeneral["radialSize"]/environmentGeneral["bins"] #radial bin size [cm]
        dz = environmentGeneral["depthSize"]/environmentGeneral["bins"] #depth bin size [cm]
        NR = environmentGeneral["bins"]                        #number of radial bins
        NZ = environmentGeneral["bins"]                        #number of depth bins
        waist = environmentGeneral["waist"]                  #1/e radius of Gaussian focus
        THRESHOLD = environmentGeneral["THRESHOLD"]          #used in roulette        
        CHANCE = environmentGeneral["CHANCE"]                #used in roulette
    	


        absorbInfo = []     #matrix of fluence rate [cm^-2] = [W/cm2 per W], former F
        escapeFlux = []     #vector of escaping flux [cm^-2] versus radial position, former J
        tempRsptot = 0   
        Atot = 0

        phot = photon.photon(iphoton)

        if mcflag == 0:
            # UNIFORM COLLIMATED BEAM INCIDENT AT z = zs 
            # Launch at (r,zz) = (radius*sqrt(rnd), 0).
            # Due to cylindrical symmetry, radial launch position is
            # assigned to x while y = 0. 
            # radius = radius of uniform beam. 
            # Initial position 
            rnd = float(1-random.random()) #random takes uniformly a number of (0,1]
            phot.update_positon(radius*math.sqrt(rnd), 0, 0)
            # Initial trajectory as cosines 
            ux = 0
            uy = 0
            uz = 1  
            # specular reflectance
            photpos = phot.get_position()        
            envStart = tree.get_env(math.sqrt(photpos[0]*photpos[0] + photpos[1]*photpos[1]), photpos[2])
            refIndex = envStart["refIndex"]
            temp   = refIndex/tissueExtMedium # refractive index mismatch
            temp   = (1.0 - temp)/(1.0 + temp)
            rsp    = temp*temp # specular reflectance at boundary 
            
        elif (mcflag == 1):
            # GAUSSIAN BEAM AT SURFACE 
            # Launch at (r,z) = (radius*sqrt(-log(rnd)), 0).
            # Due to cylindrical symmetry, radial launch position is
            # assigned to x while y = 0. 
            # radius = 1/e radius of Gaussian beam at surface. 
            # waist  = 1/e radius of Gaussian focus.
            # zfocus = depth of focal point T. 
            # Initial position 
            rnd = float(1-random.random())
            x = radius*math.sqrt(-math.log(rnd))
            phot.update_positon(x, 0, 0)

            # Initial trajectory as cosines 
            # Due to cylindrical symmetry, radial launch trajectory is
            # assigned to ux and uz while uy = 0. 
            xfocus = waist*math.sqrt(-math.log(rnd))
            temp = math.sqrt((x - xfocus)*(x - xfocus) + zfocus*zfocus)
            costheta = zfocus/temp
            # specular reflectance and refraction
            photpos = phot.get_position()        
            envStart = tree.get_env(math.sqrt(photpos[0]*photpos[0] + photpos[1]*photpos[1]), photpos[2])
            refIndex = envStart["refIndex"]
            (rsp, uz) = RFresnel(tissueExtMedium, refIndex , costheta) # new uz 
            ux  = math.sqrt(1.0 - uz*uz) # new ux 
            uy = 0.0

        elif  (mcflag == 2):
            # ISOTROPIC POINT SOURCE AT POSITION xs,ys,zs 
            # Initial position 
            phot.update_positon(xs, ys, zs)
            # Initial trajectory as cosines 
            costheta = 1.0 - 2.0*(1-random.random())
            sintheta = math.sqrt(1.0 - costheta*costheta)
            psi = 2.0*PI*(1-random.random())
            cospsi = math.cos(psi)
            if (psi < PI):
                sinpsi = math.sqrt(1.0 - cospsi*cospsi) 
            else:
                sinpsi = -math.sqrt(1.0 - cospsi*cospsi)
            ux = sintheta*cospsi
            uy = sintheta*sinpsi
            uz = costheta
            # specular reflectance
            photpos = phot.get_position() 
            envStart = tree.get_env(math.sqrt(photpos[0]*photpos[0] + photpos[1]*photpos[1]), photpos[2])
            refIndex = envStart["refIndex"]
            rsp = 0.0
        else: 
            print("choose mcflag between 0 to 2\n")
            break
        
        # set starting environment
        mus = envStart["mus"]
        mut = envStart["mua"] + envStart["mus"]
        albedo = envStart["mus"]/mut
        g = envStart["excitAnisotropy"]
        envName = envStart["name"]


        phot.set_weight(1.0 - rsp)  # set photon initial weight 
        tempRsptot += rsp # accumulate specular reflectance per photon 

        #****** HOP CYCLE **************
        # Propagate one photon until it dies by ESCAPE or ROULETTE. 
        #******************************************

        while phot.get_status():
            #*** HOP
            # Take step to new position
            # sleft = stepsize
            # ux, uy, uz are cosines of current photon trajectory
            #***
            rnd = 1-random.random()   # avoids rnd = 0 
            sleft = -math.log(rnd)/mut   # Step size

            while sleft > 0:
                s = sleft/mus				
                photPosOld = phot.get_position()
                phot.update_positon(s*ux, s*uy, s*uz)   #update position
                photPosNew = phot.get_position()
                # request once the old and new coordinates 
                rOld = math.sqrt(photPosOld[0]*photPosOld[0] + photPosOld[1]*photPosOld[1])
                hOld = photPosOld[2]
                rNew = math.sqrt(photPosNew[0]*photPosNew[0] + photPosNew[1]*photPosNew[1])
                hNew = photPosNew[2]

                sv = checkSameBin(rOld, hOld, rNew, hNew, dr, dz)

                # photon is in the same bin
                if sv:
                    #*** DROP
                    # Drop photon weight in local bin.
                    #***
                    absorb = phot.get_weight()*(1 - albedo)     # photon weight absorbed at this step 
                    phot.update_weight(-absorb)                 # decrement WEIGHT by amount absorbed 
                    Atot += absorb                              # accumulate absorbed photon weight 
                    
                    # deposit power in cylindrical coordinates z,r 
                    pos = phot.get_position()
                    r  = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])       # current cylindrical radial position 
                    ir = round((r/dr) + 1)        
                    iz = round((abs(pos[2])/dz) + 1)  
                    ir = min(ir, NR)    # round to 1 <= ir, iz 
                    iz = min(iz, NZ)    # last bin is for overflow 
                    absorbInfo.append((iz*NR + ir, absorb))             # drop absorbed weight o bin 
                    
                    sleft = 0
                
                #photon has crossed bin boundary
                else:
                    if hNew <= 0:
                        rf, uz1 = RFresnel(refIndex, tissueExtMedium, -uz)
                        rnd = 1-random.random()   # avoids rnd = 0 
                        if (rnd > rf):
                            phot.update_positon(-s*ux, -s*uy, -s*uz)   #return to original positions
                            # Photon escapes at external angle, uz1 = cos(angle) 
                            s  = abs(hNew/uz) # calculate stepsize to reach surface 
                            phot.update_positon(s*ux, s*uy, s*uz)
                            pos = phot.get_position()
                            r = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])   # find radial position r 
                            ir = round((r/dr) + 1) # round to 1 <= ir 
                            ir = min(ir,NR)  # ir = NR is overflow bin 
                            escapeFlux.append((ir, phot.get_weight()))      # increment escaping flux 
                            phot.update_dead()
                            break
                            
                        phot.update_positon(0, 0, -hNew)
                        uz = -uz
                    #check if both bins have the same environment
                    newEnv = tree.get_env(rNew, hNew)
                    if newEnv["name"] == envName:
                        #*** DROP
                        # * Drop photon weight in local bin.
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
                        absorbInfo.append((iz*NR + ir, absorb))          # DROP absorbed weight o bin 
                        
                        sleft = 0

                    else:
                        phot.update_positon(-s*ux, -s*uy, -s*uz)   #update to old positions

                        # take step to next boundary
                        s = ls + find_next_border(rOld, hOld, dr, dz, math.sqrt(ux*ux + uy*uy), uz)
                        phot.update_positon(s*ux, s*uy, s*uz)
                        
                        absorb = phot.get_weight()*(1 - albedo)     # photon weight absorbed at this step 
                        phot.update_weight(-absorb)                 # decrement WEIGHT by amount absorbed 
                        Atot += absorb                              # accumulate absorbed photon weight 
                        # deposit power in cylindrical coordinates z,r 
                        pos = phot.get_position()
                        r  = math.sqrt(pos[0]*pos[0] + pos[1]*pos[1])         # current cylindrical radial position 
                        ir = round((r/dr) + 1)        
                        iz = round((abs(pos[2])/dz) + 1)  
                        ir = min(ir, NR)    # round to 1 <= ir, iz 
                        iz = min(iz, NZ)    # last bin is for overflow 
                        absorbInfo.append((iz*NR + ir, absorb))          # DROP absorbed weight o bin 

                        sleft = (sleft-s)*mut
                        if sleft <= ls: 
                            sleft = 0
                    
                        #update positions until the next 
                        phot.update_positon(s*ux, s*uy, s*uz)
                                                
                        #get and update the environment variables
                        newEnv = tree.get_env(r, pos[2])
                        
                        mut = newEnv["mua"] + newEnv["mus"]
                        albedo = newEnv["mus"]/mut
                        g = newEnv["excitAnisotropy"]
                        envName = newEnv["name"]
                        refIndex = newEnv["refIndex"]

            #*** SPIN 
            # Scatter photon on new trajectory defined by theta and psi.
            # Theta is specified by cos(theta), which is determined 
            # based on the Henyey-Greenstein scattering function.
            # Convert theta and psi o cosines ux, uy, uz. 
            #***
            # Sample for costheta 
            rndCos = 1-random.random()
            if (g == 0.0):
                costheta = 2.0*rndCos - 1.0
            elif (g == 1.0):
                costheta = 1.0
            else:
                temp = (1.0 - g*g)/(1.0 - g + 2*g*rndCos)
                costheta = (1.0 + g*g - temp*temp)/(2.0*g)
            sintheta = math.sqrt(1.0 - costheta*costheta)	

            # Sample psi. 
            psi = 2.0*PI*(1-random.random())
            cospsi = math.cos(psi)
            if (psi < PI):
                sinpsi = math.sqrt(1.0 - cospsi*cospsi)
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
            # If photon weight below THRESHOLD, then terminate photon using
            # Roulette technique. Photon has CHANCE probability of having 
            # its weight increased by factor of 1/CHANCE,
            # and 1-CHANCE probability of terminating.
            #***
            if (phot.get_weight() < THRESHOLD):
                rnd = 1-random.random()
                if (rnd <= CHANCE):
                    phot.set_weight(phot.get_weight() / CHANCE)
                    #W /= CHANCE
                else:
                    phot.update_dead()
        

        #return absorbInfo, escapFlux, tempRspot, Atot
        # [[(1D-number,absorbValue), ...], [(number, weight), ...], value, value]
        queue_result.put([absorbInfo, escapeFlux, tempRsptot, Atot])


def writer(count, num_of_reader_procs, queue):
    """Write egers o the queue.  A reader_proc() will read them from the queue"""
    for ii in range(0, count):
        queue.put(ii)  # Put 'count' numbers to queue

    ### Tell all readers to stop
    for ii in range(0, num_of_reader_procs+1):
        queue.put("DONE")


# sort function sorts the results of the launchers, while these are processing
def sort(queue, processes, escapeFlux, absorbInfo, tempRsptot, Atot): 
    seen_processes = 0
    # sorter stops, when all launcher have send 'DONE'
    while seen_processes < processes:
        #launchReturn contains absorb Information, escaped photons, absorbed photon weight and specular reflectance
        launchReturn = queue.get()
        if launchReturn == "DONE":
            print("DONE %s" % seen_processes)
            seen_processes += 1
        else:
            for i in launchReturn[0]:
                if i[1]<0:
                    print(i[0])
                absorbInfo[i[0]] += i[1]
            for j in launchReturn[1]:
                escapeFlux[j[0]] += j[1]
            tempRsptot += launchReturn[2]
            Atot += launchReturn[3]
    print("sort DONE")
    return [absorbInfo, escapeFlux, tempRsptot, Atot]

# boolean function checks, if both coordinates lay in the same bin
def checkSameBin( x1, y1,  x2,  y2,  dx, dy):
    return (x1/dx == x2/dx and y1/dy == y2/dy)

# calculates the distance to the next bin boundary
def find_next_border(r1, z1, dr, dz, ur, uz):
	
    ir1 = r1/dr
    iz1 = z1/dz
    
    if ur>=0:
        ir2 = ir1+1
    else:
        ir2 = ir1
    
    if uz>=0:
        iz2 = iz1+1
    else:
        iz2 = iz1

    rs = 0
    zs = 0
    if ur != 0:
        rs = abs((ir2*dr - r1)/ur)
    if uz != 0:
        zs = abs((iz2*dz - z1)/uz)
    
    s = min(rs, zs)
    
    return s


def RFresnel(n1,		# incident refractive index.
            n2,		# transmit refractive index.
            ca1): 	# poer to the cosine 
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
            #  check for total ernal reflection. 
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



def SaveFile(Nfile,  J,  F,  S,  A,  E, environmentGeneral, Nphotons):

    mcflag = environmentGeneral["mcflag"]        #0 = collimated uniform, 1 = Gaussian, 2 = isotropic po
    radius = environmentGeneral["radius"]        #radius of beam (1/e width if Gaussian) (if mcflag < 2)
    xs = environmentGeneral["xs"]                #used if mcflag = 2, isotropic pt source
    ys = environmentGeneral["ys"]                #used if mcflag = 2, isotropic pt source
    zs = environmentGeneral["zs"]                #used if mcflag = 2, isotropic pt source
    dr = environmentGeneral["radialSize"]/environmentGeneral["bins"] #radial bin size [cm]
    dz = environmentGeneral["depthSize"]/environmentGeneral["bins"] #depth bin size [cm]
    NR = environmentGeneral["bins"]                        #number of radial bins
    NZ = environmentGeneral["bins"]                        #number of depth bins
    waist = environmentGeneral["waist"]                  #1/e radius of Gaussian focus

    print("mcOUT%d.dat" % Nfile)
    file = open("mcOUT" + str(Nfile) + ".dat", "w")

    # pr run parameters 
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

    # pr SAE values 
    file.write("%1.6e\tSpecular reflectance\n" % S)
    file.write("%1.6e\tAbsorbed fraction\n" % A)
    file.write("%1.6e\tEscaping fraction\n" % E)

    # pr r[ir] to row 
    file.write("%0.1f" % 0.0) # ignore upperleft element of matrix 
    for ir in range(1, NR+1):
        r2 = dr*ir
        r1 = dr*(ir-1)
        r = 2.0/3*(r2*r2 + r2*r1 + r1*r1)/(r1 + r2)
        file.write("\t%1.5f" % r)
    file.write("\n")

    # pr J[ir] to next row 
    file.write("%0.1f" % 0.0) # ignore this first element of 2nd row 
    for ir in range(1, NR+1):	
        file.write("\t%1.12e" % J[ir])
    file.write("\n")

    # prf z[iz], F[iz][ir] to remaining rows 
    for iz in range(1, NR+1):
        z = (iz - 0.5)*dz # z values for depth position in 1st column 
        file.write("%1.5f" % z)
        for ir in range(1, NR+1):
            file.write("\t %1.6e" % F[iz*NR + ir])
        file.write("\n")
    file.close()


if __name__ == '__main__':
    startTime = time.time()

    # set up the environment manager
    environmentGeneral = input_data.environmentGeneral
    envManager = manageEnv.manageEnv()
    for env in input_data.envDetail:
        envManager.add_environment(env["name"], env["mua"], env["mus"], env["excitAnisotropy"], env["refractiveIndex"], env["height"], env["radius"], env["default"], environmentGeneral["bins"]/environmentGeneral["depthSize"])

    Nfile = 0       # saves as mcOUTi.dat, where i = Nfile 
    # Run parameters 
    Nphotons = input_data.Nphotons # number of photons to be launched 
    PROUT = 1

    # data variables
    escapeFlux       = np.zeros(environmentGeneral["bins"]+1)        # for escaping flux 
    absorbInfo       = np.zeros((environmentGeneral["bins"]+1)*(environmentGeneral["bins"]+1)) # for absorbed fluence rate 
    tempRsptot = 0.0 
    Atot   = 0.0 

    
    # create envList
    pixel = environmentGeneral["bins"]
    createEnvListTimeStart = time.time()
    envList = [envManager.get_default_variables()]*(pixel*pixel)
    envList = envManager.assign_env(envList, pixel)     
    createEnvListTimeEnd = time.time()
    envListTime = createEnvListTimeEnd-createEnvListTimeStart

    # set quadtree
    quadTree = QuadTree()
    pixelhalf = pixel/2
    quadTree.create_Tree(-pixelhalf, pixelhalf, -pixelhalf, pixelhalf, pixel, math.log(pixel, 2), envList)
    quadTreeTimeEnd = time.time()
    quadTreeTime = quadTreeTimeEnd-createEnvListTimeEnd

    # create queues for communication with multiprocessing
    queue_result = Queue()
    queue_photons = Queue()

    launcherTimeStart = time.time() 

    # multiprocessing and launcher
    processes = cpu_count()
    launchProcesses = processes

    all_launcher_procs = []
    for i in range(0, launchProcesses):
        launcher_p = Process(target=launch, args=(queue_photons, environmentGeneral, queue_result, i, quadTree,))
        launcher_p.daemon = True
        launcher_p.start()
        print("launcher " + str(i)+ " started")

        all_launcher_procs.append(launcher_p)

    launcherTimeEnd = time.time()
    
    launcherTime = launcherTimeEnd - launcherTimeStart
    writer(Nphotons, launchProcesses, queue_photons)    # creates all photons in a queue
    return_dict = sort(queue_result, launchProcesses, escapeFlux, absorbInfo, tempRsptot, Atot) # sorts all collected data of the launchers
    
    # ensure, that all launcher have finished
    for idx, a_launcher_proc in enumerate(all_launcher_procs):
        print("    Waiting for reader_p.join() index %s" % idx)

        a_launcher_proc.join()  # Wait for a_launcher_proc() to finish

        print("        reader_p() idx:%s is done" % idx)

    absorbInfo = return_dict[0]
    escapeFlux = return_dict[1]
    tempRsptot = return_dict[2]
    Atot = return_dict[3]



    # Normalize
    #   J[ir]      escaping flux density [W/cm^2 per W incident] 
    #              where bin = 2.0*PI*r[ir]*dr [cm^2].
    #   F[iz][ir]  fluence rate [W/cm^2 per W incident] 
    #              where bin = 2.0*PI*r[ir]*dr*dz [cm^3].
    temp = 0.0
    NR = environmentGeneral["bins"]
    NZ = environmentGeneral["bins"]
    dr = environmentGeneral["radialSize"]/NR #radial bin size [cm]
    dz = environmentGeneral["depthSize"]/NR #depth bin size [cm]
    PI = environmentGeneral["PI"]

    for ir in range(1, NR+1):
        r = (ir - 0.5)*dr
        temp += escapeFlux[ir]    # accumulate total escaped photon weight 
        escapeFlux[ir] /= 2.0*PI*r*dr*Nphotons # flux density

    for iz in range(1, NZ+1):
        mua = envList[(iz-1)*NR]["mua"]
        for ir in range(1, NR+1):
            r = (ir - 0.5)*dr        	                 
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
    print("Time to create env list (sec): " + str(envListTime))
    print("Time to create quadTree (sec): " + str(quadTreeTime))
    print("Time to start all launcher (sec): " + str(launcherTime))
    print("Over all time (sec): " + str(endTime - startTime))
