import math
import random
import numpy as np


'''**********************************************************
#mcsubLIB.c
# Monte Carlo simulation of fluence rate F and escaping flux J
# in a semi-infinite medium such as biological tissue,
# with an external_medium/tissue surface boundary.
# *
 * Contains subroutines,
 *		mcsub()
 *		RFresnel()
 *		SaveFile()
 *		RandomGen()
 *	and memory allocation routines,
 *		nerror()
 *		*AllocVector()
 *		**AllocMatrix()
 *		FreeVector()
 *		FreeMatrix()
 **********************************************************'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#/*************
# * callmcsub.c
# * A calling program that
# *	1. defines parameters for Monte Carlo run
# *	2. calls the mcsub() routine
# *	3. saves the results into an output file using SaveFile()
# *************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mcsubLIB.h"

#/*************************/
#/**** USER CHOICES *******/
#/*************************/
#define BINS        101      /* number of bins, NZ and NR, for z and r */
#/*************************/
#/*************************/

if __name__ == '__main__':
	#/**********************
	# * MAIN PROGRAM
	# *********************/

	#/*************************/
	#/**** USER CHOICES *******/
	#/*************************/
	#/* number of file for saving */
	Nfile = 0       #/* saves as mcOUTi.dat, where i = Nfile */
	#/* tissue parameters */
	BINS = 101      #/* number of bins, NZ and NR, for z and r */
	mua  = 1.0      #/* excitation absorption coeff. [cm^-1] */
	mus  = 100      #/* excitation scattering coeff. [cm^-1] */
	g    = 0.90     #/* excitation anisotropy [dimensionless] */
	n1   = 1.40     #/* refractive index of medium */
	n2   = 1.00     #/* refractive index outside medium */
	#/* beam parameters */
	mcflag = 0    #/* 0 = collimated, 1 = focused Gaussian, 
							#2 >= isotropic pt. */
	radius = 0    #/* used if mcflag = 0 or 1 */
	waist  = 0.10 #/* used if mcflag = 1 */
	zfocus = 1.0    #/* used if mcflag = 1 */
	xs     = 0.0    #/* used if mcflag = 2, isotropic pt source */
	ys     = 0.0    #/* used if mcflag = 2, isotropic pt source */
	zs     = 0.0    #/* used if mcflag = 2, isotropic pt source, or mcflag = 0 collimated*/
	boundaryflag = 1 #/* 0 = infinite medium, 1 = air/tissue surface boundary
	#/* Run parameters */
	Nphotons = 10000 #/* number of photons to be launched */
	dr     = 0.0020 #/* radial bin size [cm] */
	dz     = 0.0020 #/* depth bin size [cm] */
	PRINTOUT = 1
	#/*************************/
	#/****** Setup output parameters, vectors, arrays **********/
	#S    /* specular reflectance at air/tissue boundary */
	#A    /* total fraction of light absorbed by tissue */
	#E    /* total fraction of light escaping tissue */
	#J    /* escaping flux, J[ir], [W/cm2 per W incident] */
	#F    /* fluence rate, F[iz][ir], [W/cm2 per W incident] */
	NR = BINS    #/* number of radial bins */
	NZ = BINS    #/* number of depth bins */
	J       = np.zeros(NR)        #/* for escaping flux */
	F       = np.zeros([NZ, NR]) #/* for absorbed fluence rate */
	#/*************************/
	#/*************************/


	#// choose Nphotons
	THRESH = 1e-4
	albedo = mus/(mua + mus)
	Nsteps = math.log(THRESH)/math.log(albedo)
	tperstep = 249e-9
	t = 30
	Nphotons = t/(tperstep*Nsteps)
	print("Nphotons = %5.4e\n", Nphotons)


	#/*************************/
	#/*************************/

	#**********************
	# DECLARE SUBROUTINES:
	#*********************/
	
	'''*****
	* The Monte Carlo subroutine mcsub()
	*	Tissue properties:
	*			mua       = absorption coefficient [cm^-1]
	*			mus       = scattering coefficient [cm^-1]
	*			g         = anisotropy of scattering [dimensionless]
	*			n1        = refractive index of tissue
	*			n2        = refractive index of external medium
	*	Incident beam characteristics:
	*			mcflag:   0 = collimated uniform, 1 = Gaussian, 2 = isotropic point
	*			xs,ys,zs  = position of istropic point source (if mcflag = 2)
	*			boundaryflag = 1 if air/tissue surface, = 0 if infinite medium
	*			radius    = radius of beam (1/e width if Gaussian) (if mcflag < 2)
	*			waist     = 1/e radius of focus (if Gaussian (if mcflag = 1)
	*			zfocus    = depth of focus (if Gaussian (if mcflag = 1)
	*	OUTPUT:
	*			J[ir]     = vector of escaping flux [cm^-2] versus radial position
	*			F[iz][ir] = matrix of fluence rate [cm^-2] = [W/cm2 per W]
	*			Sptr      = pointer to S, specular refelctance
	*			Aptr      = pointer to A, total absorbed fraction
	*			Eptr      = pointer to E, total escaping fraction 
	*'''




	def mcsub( mua,  mus,  g,  n1,  n2, 
				NR, NZ,  dr,  dz,  Nphotons,
				mcflag,  xs,  ys,  zs, boundaryflag,
				radius,  waist,  zfocus,
				J,  F,  Sptr,  Aptr,  Eptr,
				PRINTOUT):
		test=1

	#* Computes internal reflectance at tissue/air interface */

	def RFresnel( n1,  n2,  ca1,  *ca2_Ptr):
		test=1


	'''* Saves OUTPUT file, "mcOUTi.dat" where i = Nfile.
	* The INPUT parameters are saved
	* surface escape R(r) and fluence rate distribution F(z,r) 
	* SAVE RESULTS TO FILES 
	*   to files named �Ji.dat� and �Fi.dat� where i = Nfile.
	* Saves �Ji.dat� in following format:
	*   Saves r[ir]  values in first  column, (1:NR,1) = (rows,cols).
	*   Saves Ji[ir] values in second column, (1:NR,2) = (rows,cols).
	*   Last row is overflow bin.
	* Saves �Fi.dat� in following format:
	*   The upper element (1,1) is filled with zero, and ignored.
	*   Saves z[iz] values in first column, (2:NZ,1) = (rows,cols).
	*   Saves r[ir] values in first row,    (1,2:NZ) = (rows,cols).
	*   Saves Fi[iz][ir] in (2:NZ+1, 2:NR+1).
	*   Last row and column are overflow bins.
	* Saves "mcSAEi.dat" as three tab-delimited values in one row = [S A E].
	*'''
	#void SaveFile(int Nfile,  *J,  **F,  S,  A,  E, 
	#	 mua,  mus,  g,  n1,  n2, 
	#	short mcflag,  radius,  waist,  xs,  ys,  zs, 
	#	short NR, short NZ,  dr,  dz,  Nphotons)

	#* Random number generator 
	#   Initiate by RandomGen(0,1,NULL)
	#   Use as rnd = RandomGen(1,0,NULL) */
	#RandomGen(char Type, long Seed, long *Status) 

	#/* Memory allocation routines 
	# * from MCML ver. 1.0, 1992 L. V. Wang, S. L. Jacques,
	# * which are modified versions from Numerical Recipes in C. */
	#void   nrerror(char error_text[])
	# *AllocVector(short nl, short nh)
	# **AllocMatrix(short nrl,short nrh,short ncl,short nch)
	#void   FreeVector( *v,short nl,short nh)
	#void   FreeMatrix( **m,short nrl,short nrh,short ncl,short nch)



	#/**********************************************************
	# *             list SUBROUTINES:
	# **********************************************************/

	#/**********************************************************
	# * The Monte Carlo SUBROUTINE
	# **********************************************************/
	def mcsub( mua,  mus,  g,  n1,  n2, 
				NR, NZ,  dr,  dz,  Nphotons,
				mcflag,  xs,  ys,  zs, boundaryflag,
				radius,  waist,  zfocus,
				J,  F,  Sptr,  Aptr,  Eptr,
				PRINTOUT):
		test=1
	#/* Constants */
	PI          = 3.1415926
	ALIVE       = 1          #/* if photon not yet terminated */
	DEAD        = 0           #/* if photon is to be terminated */
	THRESHOLD   = 0.0001        #/* used in roulette */
	CHANCE      = 0.1          #/* used in roulette */

	#/* Variable parameters */


	#/**** INITIALIZATIONS *****/
	CNT = 0
	mut    = mua + mus
	albedo = mus/mut
	Rsptot = 0.0 #/* accumulate specular reflectance per photon */
	Atot   = 0.0 #/* accumulate absorbed photon weight */


	#/* initialize arrays to zero */
	for ir in range(1, NR+1):
		J[ir] = 0.0
		for iz in range(1, NZ+1):
			F[iz][ir] = 0.0

	#/*============================================================
	#======================= RUN N photons =====================
	# * Launch N photons, initializing each one before progation.
	#============================================================*/
		
	#/**** LAUNCH 
	#   Initialize photon position and trajectory.
	#   Implements an isotropic point source.
	#*****/

	if mcflag == 0:
		#/* UNIFORM COLLIMATED BEAM INCIDENT AT z = zs */
		#/* Launch at (r,zz) = (radius*sqrt(rnd), 0).
		# * Due to cylindrical symmetry, radial launch position is
		# * assigned to x while y = 0. 
		# * radius = radius of uniform beam. */
		#/* Initial position */
		rnd = float(random.random()) #random takes uniformly a number of (0,1]
		x = radius*math.sqrt(rnd)
		y = 0
		z = zs
		#/* Initial trajectory as cosines */
		ux = 0
		uy = 0
		uz = 1  
		#/* specular reflectance */
		temp   = n1/n2 #/* refractive index mismatch, internal/external */
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
		y = 0.0
		z = 0.0
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
		(rsp, uz) = RFresnel(n2, n1, costheta, uz) #/* new uz */
		ux  = math.sqrt(1.0 - uz*uz) #/* new ux */
	elif  (mcflag == 2):
		#/* ISOTROPIC POINT SOURCE AT POSITION xs,ys,zs */
		#/* Initial position */
		rnd = float(random.random())
		x = xs
		y = ys
		z = zs
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
		print("choose mcflag between 0 to 3\n")
		

	W = 1.0 - rsp  #/* set photon initial weight */
	Rsptot += rsp #/* accumulate specular reflectance per photon */
	photon_status = ALIVE

	#/******************************************
	#****** HOP_ESCAPE_SPINCYCLE **************
	#* Propagate one photon until it dies by ESCAPE or ROULETTE. 
	#*******************************************/

	while (photon_status == ALIVE):
		#/**** HOP
		# * Take step to new position
		# * s = stepsize
		# * ux, uy, uz are cosines of current photon trajectory
		# *****/
		rnd = random.Random()   #/* avoids rnd = 0 */
		s = -math.log(rnd)/mut   #/* Step size.  Note: log() is base e */
		x += s*ux           #/* Update positions. */
		y += s*uy
		z += s*uz

		#/* Does photon ESCAPE at surface? ... z <= 0? */
		if ( (boundaryflag == 1) & (z <= 0)):
			rnd = random.random()
			#/* Check Fresnel reflectance at surface boundary */
			rf, uz1 = RFresnel(n1, n2, -uz, uz1)
			if (rnd > rf): 
				#/* Photon escapes at external angle, uz1 = cos(angle) */
				x -= s*ux       #/* return to original position */
				y -= s*uy
				z -= s*uz
				s  = abs(z/uz) #/* calculate stepsize to reach surface */
				x += s*ux       #/* partial step to reach surface */
				y += s*uy
				r = math.sqrt(x*x + y*y)   #/* find radial position r */
				ir = (r/dr) + 1 #/* round to 1 <= ir */
				if (ir > NR): ir = NR  #/* ir = NR is overflow bin */
				J[ir] += W      #/* increment escaping flux */
				photon_status = DEAD
				
			else:
				z = -z   #/* Total internal reflection. */
				uz = -uz
				

		if (photon_status  == ALIVE):
			#/*********************************************
			# ****** SPINCYCLE = DROP_SPIN_ROULETTE ******
			# *********************************************/

			#/**** DROP
			# * Drop photon weight (W) into local bin.
			# *****/
			absorb = W*(1 - albedo)       #/* photon weight absorbed at this step */
			W -= absorb                  #/* decrement WEIGHT by amount absorbed */
			Atot += absorb               #/* accumulate absorbed photon weight */
			#/* deposit power in cylindrical coordinates z,r */
			r  = math.sqrt(x*x + y*y)         #/* current cylindrical radial position */
			ir = (r/dr) + 1        #/* round to 1 <= ir */
			iz = (abs(z)/dz) + 1  #/* round to 1 <= iz */
			if (ir >= NR): ir = NR        #/* last bin is for overflow */
			if (iz >= NZ): iz = NZ        #/* last bin is for overflow */
			F[iz][ir] += absorb          #/* DROP absorbed weight into bin */

			#/**** SPIN 
			# * Scatter photon into new trajectory defined by theta and psi.
			# * Theta is specified by cos(theta), which is determined 
			# * based on the Henyey-Greenstein scattering function.
			# * Convert theta and psi into cosines ux, uy, uz. 
			# *****/
			#/* Sample for costheta */
			rnd = random.Random()
			if (g == 0.0):
				costheta = 2.0*rnd - 1.0
			elif (g == 1.0):
				costheta = 1.0
			else:
				temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
				costheta = (1.0 + g*g - temp*temp)/(2.0*g)
			sintheta = math.sqrt(1.0 - costheta*costheta)	#/*sqrt faster than sin()*/

			#/* Sample psi. */
			psi = 2.0*PI*random.Random()
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
			if (W < THRESHOLD):
				rnd = random.Random()
				if (rnd <= CHANCE):
					W /= CHANCE
				else:
					photon_status = DEAD

			#/**********************************************
			#  **** END of SPINCYCLE = DROP_SPIN_ROULETTE *
			#  **********************************************/

		
	#/******************************************
	# ****** END of HOP_ESCAPE_SPINCYCLE ******
	# ****** when photon_status == DEAD) ******
	# ******************************************/

	#/* If photon dead, then launch new photon. */
	#/*======================= End RUN N photons =====================
	#=====================================================================*/

	#/************************
	# * NORMALIZE 
	# *   J[ir]      escaping flux density [W/cm^2 per W incident] 
	# *              where bin = 2.0*PI*r[ir]*dr [cm^2].
	# *	  F[iz][ir]  fluence rate [W/cm^2 per W incident] 
	# *              where bin = 2.0*PI*r[ir]*dr*dz [cm^3].
	# ************************/
	temp = 0.0
	for ir in range(1, NR+1):
		r = (ir - 0.5)*dr
		temp += J[ir]    #/* accumulate total escaped photon weight */
		J[ir] /= 2.0*PI*r*dr*Nphotons
		for iz in range(1, NZ+1):	                #/* flux density */
			F[iz][ir] /= 2.0*PI*r*dr*dz*Nphotons*mua #/* fluence rate */

	Sptr = S = Rsptot/Nphotons
	Aptr = A = Atot/Nphotons
	Eptr = E = temp/Nphotons

	print("Nphotons = %5.1e\n", Nphotons)
	print("Specular = %5.6f\n", S)
	print("Absorbed = %5.6f\n", A)
	print("Escaped  = %5.6f\n", E)
	print("total    = %5.6f\n", S+A+E)
	#/******** END SUBROUTINE **********/
 

#/***********************************************************
# *	FRESNEL REFLECTANCE
# * Computes reflectance as photon passes from medium 1 to 
# * medium 2 with refractive indices n1,n2. Incident
# * angle a1 is specified by cosine value ca1 = cos(a1).
# * Program returns value of transmitted angle a1 as
# * value in *ca2_Ptr = cos(a2).
# ****/
def RFresnel(n1,		#/* incident refractive index.*/
			n2,		#/* transmit refractive index.*/
			ca1,		#/* cosine of the incident */									#/* angle a1, 0<a1<90 degrees. */
			ca2_Ptr): 	#/* pointer to the cosine */
									#/* of the transmission */
									#/* angle a2, a2>0. */
	r

	if(n1==n2): #/** matched boundary. **/
		ca2_Ptr = ca1
		r = 0.0
	elif(ca1>(1.0 - 1.0e-12)): #/** normal incidence. **/
		ca2_Ptr = ca1
		r = (n2-n1)/(n2+n1)
		r *= r
	elif(ca1< 1.0e-6):	#/** very slanted. **/
		ca2_Ptr = 0.0
		r = 1.0
	else:	  		#/** general. **/
		sa1, sa2 #/* sine of incident and transmission angles. */
		ca2      #/* cosine of transmission angle. */
		sa1 = math.sqrt(1-ca1*ca1)
		sa2 = n1*sa1/n2
		if(sa2>=1.0):	
			#/*  check for total internal reflection. */
			ca2_Ptr = 0.0
			r = 1.0
		else:
			cap, cam	#/* cosines of sum ap or diff am of the two */
								#/* angles: ap = a1 + a2, am = a1 - a2. */
			sap, sam	#/* sines. */
			ca2_Ptr = ca2 = math.sqrt(1-sa2*sa2)
			cap = ca1*ca2 - sa1*sa2 #/* c+ = cc - ss. */
			cam = ca1*ca2 + sa1*sa2 #/* c- = cc + ss. */
			sap = sa1*ca2 + ca1*sa2 #/* s+ = sc + cs. */
			sam = sa1*ca2 - ca1*sa2 #/* s- = sc - cs. */
			r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam) 
			#/* rearranged for speed. */
	return r, ca2_Ptr
#/******** END SUBROUTINE **********/


#/***********************************************************
# * SAVE RESULTS TO FILES 
#***********************************************************/
def SaveFile(Nfile,  J,  F,  S,  A,  E, 
	mua,  mus,  g,  n1,  n2, 
	mcflag,  radius,  waist,  xs,  ys,  zs, 
	NR, NZ,  dr,  dz,  Nphotons):

	print("mcOUT%d.dat", Nfile)
	file = open("mcOUT" + Nfile + ".dat", "w")

	#/* print run parameters */
	file.write("%0.3e\tmua, absorption coefficient [1/cm]\n", mua)
	file.write("%0.4f\tmus, scattering coefficient [1/cm]\n", mus)
	file.write("%0.4f\tg, anisotropy [-]\n", g)
	file.write("%0.4f\tn1, refractive index of tissue\n", n1)
	file.write("%0.4f\tn2, refractive index of outside medium\n", n2)
	file.write("%d\tmcflag\n", mcflag)
	file.write("%0.4f\tradius, radius of flat beam or 1/e radius of Gaussian beam [cm]\n", radius)
	file.write("%0.4f\twaist, 1/e waist of focus [cm]\n", waist)
	file.write("%0.4f\txs, x position of isotropic source [cm]\n", xs)
	file.write("%0.4f\tys, y\n", ys)
	file.write("%0.4f\tzs, z\n", zs)
	file.write("%d\tNR\n", NR)
	file.write("%d\tNZ\n", NZ)
	file.write("%0.5f\tdr\n", dr)
	file.write("%0.5f\tdz\n", dz)
	file.write("%0.1e\tNphotons\n", Nphotons)

	#/* print SAE values */
	file.write("%1.6e\tSpecular reflectance\n", S)
	file.write("%1.6e\tAbsorbed fraction\n", A)
	file.write("%1.6e\tEscaping fraction\n", E)

	#/* print r[ir] to row */
	file.write("%0.1f", 0.0) #/* ignore upperleft element of matrix */
	for ir in range(1, NR+1):
		r2 = dr*ir
		r1 = dr*(ir-1)
		r = 2.0/3*(r2*r2 + r2*r1 + r1*r1)/(r1 + r2)
		file.write("\t%1.5f", r)
	file.write("\n")

	#/* print J[ir] to next row */
	file.write("%0.1f", 0.0) #/* ignore this first element of 2nd row */
	for ir in range(1, NR+1):	
		file.write("\t%1.12e", J[ir])
	file.write("\n")

	#/* printf z[iz], F[iz][ir] to remaining rows */
	for iz in range(1, NR+1):
		z = (iz - 0.5)*dz #/* z values for depth position in 1st column */
		file.write("%1.5f", z)
		for ir in range(1, NR+1):
			file.write("\t %1.6e", F[iz][ir])
		file.write("\n")
	file.close()

