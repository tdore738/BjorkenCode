from scipy import *
from scipy import integrate
from scipy import optimize
from scipy.integrate import ode
from scipy import interpolate as terp
from pylab import *
from matplotlib import pyplot
import sys
import gc
import dill
import os
import errno
from operator import itemgetter
from matplotlib.ticker import ScalarFormatter


# These EoS files were stored in binary (I don't remember how neessary this was..)
# so they need to be opened using 'dill.load', where dill is an import that can manage
# binary stroage

#might need to change open path!!!!!!!!
f_p = dill.load(open('../../EoS/P_of_T_Mu','rb'))
f_e = dill.load(open('../../EoS/En_of_T_Mu','rb'))
f_nB = dill.load(open('../../EoS/nB_of_T_Mu','rb'))
f_c = dill.load(open('../../EoS/Sound_of_T_Mu','rb'))
# f_corr = dill.load(open('CorLen_of_T_Mu','rb')) #only used for critcal scaled
f_eta = dill.load(open('../../EoS/Eta_of_T_Mu','rb'))

ACCURACY = 1e-6
MAXITER  = 1000

# rootfinder for EoS, originally written by Chun Shen
# the search utilizes monotonicity of the EoS in (T,mu),(e,rho)
# first doing a 1-D search to find max and min values of
# temperature, then finding the right pair of T and mu

def binary_search_1d(ed_local, muB_local):
    iteration = 0
    T_min = .030; T_max = .8
    e_low = f_e(muB_local, T_min)
    e_up  = f_e(muB_local, T_max)
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(muB_local, T_mid)
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(muB_local, T_mid)
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            iteration += 1
        return(T_mid)

def binary_search_2d(ed_local, nB_local):
    iteration = 0
    muB_min = 0.0; muB_max = .600
    T_max = binary_search_1d(ed_local, muB_min)
    nB_min = f_nB(muB_min, T_max)
    T_min = binary_search_1d(ed_local, muB_max)
    nB_max = f_nB(muB_max, T_min)
    if (nB_local < nB_min):
        return(T_max, muB_min)
    elif (nB_local > nB_max):
        return(T_min, muB_max)
    else:
        muB_mid = (muB_min + muB_max)/2.
        T_mid = binary_search_1d(ed_local, muB_mid)
        nB_mid = f_nB(muB_mid, T_mid)
        abs_err = abs(nB_mid - nB_local)
        rel_err = abs_err/abs(nB_mid + nB_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (nB_local < nB_mid):
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            muB_mid = (muB_max + muB_min)/2.
            T_mid = binary_search_1d(ed_local, muB_mid)
            nB_mid = f_nB(muB_mid, T_mid)
            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/abs(nB_mid + nB_local)
            iteration += 1

        return(T_mid, muB_mid) 


# Every iteration of this code will only run one initial condition for e and rho
# these need to be passed in as arguments

rho_0 = float(sys.argv[1])
E_0 = float(sys.argv[2])
#FOtype = str(sys.argv[3])

#if (FOtype=='away'):
#	T_fo = .145
#	Mu_fo = .124
#elif(FOtype=='critical'):
#	T_fo = .14
#	Mu_fo = .348
#else:
#	sys.exit('Uncertain Freeze Out Point')


hbarC = .19733
# corr0 = 2 #arbitrarily chosen scaling, only used for critical scaling

#thermodynamic relation
def entropy(en,rho,p,Mu_B,T):
    return (en - Mu_B*rho + p)/T
VecEntropy = vectorize(entropy)

# determine the shear relaxation time from the eta/tau_pi = e+p/5 = w/5
# remember that f_eta is actually eta*T/w
def sRelax(Mu_B,T):
    return (5.*f_eta(Mu_B,T)/T)*hbarC
VecSrelax = vectorize(sRelax)


# no critical scaling for the bulk viscosity in this file
def bulkV(muB,T):
    return 36*((1/3 - f_c(muB,T))/(8.*pi))
VecBulk = vectorize(bulkV)

# a relation for zeta/tau_PI, calculated at leading order in mass corrections
# to ultra-relativistic kinetic theory; Denicol, Gale, Jeon, PRC (2014)
def bRelax(muB,T):
    return (bulkV(muB,T)/(T*15*(1/3-f_c(muB,T))**2))*hbarC
VecBrelax = vectorize(bRelax)




# The hydrodynamic equations of motion for bjorken flow
def functionVectorUNITS(t,F,rho):
#     print(F,t,a,b,c,d,e,f,g,h,i,j,k)
    
    en, S, B = F
    
    
    n = rho(t)
    (T,Mu_B) = binary_search_2d(en,n) 
#     (T,Mu_B) = (f_T(en,rho(t)),f_Mu(en,rho(t)))    

    p = f_p(Mu_B,T)
    cs2 = f_c(Mu_B,T)
    lPp = (8/5)*((1/3)-cs2)
    

    den_dt = -(en + p - S + B)/t
    
    # again we use relations from kientic theory for the ratios of transport coefficients
    # to the relaxation times, this is where the constant numbers come from
    dS_dt = -S/sRelax(Mu_B,T) + 4/(15*t)*(en+p) - (58*S)/(21*t)\
                + (4*B)/(5*t)


    # again using the ratio relations from kinetic theory
    dB_dt = -B/bRelax(Mu_B,T) - 15*(((1/3)-cs2)**2)*(en + p)/t \
            - (2*B)/(3*t) + (2*lPp*S)/(3*t)


    
    return [den_dt, dS_dt, dB_dt]




t0 = .6
tf = 500
visc0 = [(0.5,-0.5)]#,(.5,-.3),(.5,-.1),(.5,0.),(.3,0.),(.1,0.),
#(0.,0.),(-.1,0),(-.3,0.),(-.5,0.),(-.5,0.1),(-.5,0.3),(-.5,0.5)]
rho_0 = rho_0*t0
#FOregion = .0025

def Rho(t):
	return rho_0/t



# easiest to give one array with all combinations of initial shear and bulk you would like to run
for Pis in visc0:
    dt = .001
    shear0 = Pis[0]
    bulk0 = Pis[1]
    (temp0,chem0) = binary_search_2d(E_0,Rho(t0))

# anything to do with FO types is no longer done in this file, can be done after
#    FOdistance = sqrt((temp0-T_fo)**2 + (chem0-Mu_fo)**2)

#    if(FOdistance > FOregion):
#        inRegion = False
#    else:
#        inRegion = True

    shear0 = shear0*(E_0 + f_p(chem0,temp0))[0]
    bulk0 = bulk0*(E_0 + f_p(chem0,temp0))[0]

    initialCon = [E_0,shear0,bulk0]
    res = [[t0,initialCon[0],initialCon[1],initialCon[2],temp0,chem0]]

    integ = ode(functionVectorUNITS).set_integrator('vode',nsteps=50000,method='bdf')

    integ.set_initial_value(initialCon, t0)
    integ.set_f_params(Rho)


    integr = initialCon
    curTemp = temp0

    decreaseRes = False
    while((curTemp >= .110 and integ.t < tf ) and (integ.successful() )):
        if(curTemp < .135 and not decreaseRes):
            dt = 10*dt
            decreaseRes = True
        integr = integ.integrate(integ.t+dt)
        (curTemp,chem) =  binary_search_2d(integr[0],Rho(integ.t))

#        newFOd = sqrt((curTemp-T_fo)**2 + (chem-Mu_fo)**2)

#        if((newFOd < FOregion) and not inRegion):
#            inRegion = True
#        if((not inRegion) and (curTemp<=T_fo)):
#            break
#        else:
#            FOdistance = newFOd

        res.append([integ.t+dt,integr[0],integr[1],integr[2],curTemp,chem])

# This 'if' statement is leftover from when the trajectory selection was done while running the code
# this is now done in a separate piece of code
    if(True):
        results = array(res)
 #       print(rho_0/t0,'[',shear0,',',bulk0,']',': made it')

        time = results[::,0]
        enAr = results[::,1]
        sh = results[::,2]
        bu = results[::,3]
        rhoAr = array([Rho(time[i]) for i in range(time.size)])


        tempOft = results[::,4]
        chemOft = results[::,5]
        temptoFourth = np.power(tempOft/.19733,4)
        presrOft = terp.dfitpack.bispeu(f_p.tck[0], f_p.tck[1], f_p.tck[2], f_p.tck[3], f_p.tck[4], chemOft, tempOft)[0]
        entrOft = VecEntropy(enAr,rhoAr,presrOft,chemOft,tempOft)
        etOft = terp.dfitpack.bispeu(f_eta.tck[0], f_eta.tck[1], f_eta.tck[2], f_eta.tck[3], f_eta.tck[4], chemOft, tempOft)[0]#VecEta(enAr,chemOft,tempOft)
        bOft = VecBulk(chemOft,tempOft)
        tau_piOft = VecSrelax(chemOft,tempOft)
        tau_PIOft = VecBrelax(chemOft,tempOft)

        chi = sh/(enAr + presrOft)
        omega = bu/(enAr + presrOft)

        sndOft = terp.dfitpack.bispeu(f_c.tck[0], f_c.tck[1], f_c.tck[2], f_c.tck[3], f_c.tck[4], chemOft, tempOft)[0]
        effPres = array([omega[i]+(presrOft[i]/(enAr[i]+presrOft[i]))-chi[i] for i in range(time.size)])
        trPres = array([omega[i]+(presrOft[i]/(enAr[i]+presrOft[i]))+(chi[i]/2) for i in range(time.size)])
        plOVpt = effPres/trPres
        enArFM = enAr/.19733

        shearLabel = round(chi[0],2)
        bulkLabel = round(omega[0],2)
        rhoLabel = round(rho_0/t0,3)


        typeDir = "../DNMR"    
#        FOpointDir = "{0}/FO:[{1},{2}]".format(typeDir,Mu_fo,T_fo)                    
        bulkDir = "{0}/NoCS".format(typeDir)
        EnDir = "{0}/En0:{1}".format(bulkDir,E_0)
        viscDir = "{0}/{1}_{2}".format(EnDir,shearLabel,bulkLabel)
        rhoDir = "{0}/{1}".format(viscDir,rhoLabel)

        saveDirec = '../DNMR/NoCS/En0:{0}/{1}_{2}/{3}'.format(E_0,shearLabel,bulkLabel,rhoLabel)

        os.makeDirs(saveDirec,exist_ok=true)



        plot(time,etOft)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\eta/s$',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/etaVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,bOft)
        pyplot.xlabel('Time ',size=20)
        pyplot.ylabel(r'$\zeta/s$',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/zetaVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,entrOft)
        pyplot.xlabel('Time ',size=20)
        pyplot.ylabel('Entropy (fm^-3)',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/entropyVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,sndOft)
        pyplot.xlabel('Time ',size=20)
        pyplot.ylabel('Sound',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/soundVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,presrOft/temptoFourth)
        # plot(time,presrOft)
        pyplot.xlabel('Time ',size=20)
        pyplot.ylabel(r'Pressure $(GeV/fm^3$)',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/presVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,effPres)
        pyplot.xlabel('Time ',size=20)
        pyplot.ylabel(r'$\Pi + p - \pi$',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/effectivePressureVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,trPres)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\Pi + p + \pi/2$',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/trPressureVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,plOVpt)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$p_L/p_T$',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/plOVpt.png".format(rhoDir))
        pyplot.clf()
        plot(time,tempOft)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel('T (GeV)',size = 26)
        pyplot.title(r'$\chi_0=${0} $\Omega_0=${1} $\rho_0=${2}'.format(shearLabel,bulkLabel,rhoLabel),size = 18)
        pyplot.savefig("{0}/tempVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,chemOft)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\mu_B$ (GeV)',size = 26)
        pyplot.title(r'En0={0} Rho0={1}'.format(E_0,Rho(t0)),size = 22)
        pyplot.savefig("{0}/chemVtime.png".format(rhoDir))
        pyplot.clf()
        plot(chemOft,tempOft)
        pyplot.xlabel(r'$\mu_B$ (GeV)',size=20)
        pyplot.ylabel('T (GeV)',size = 26)
        pyplot.title(r'En0={0} Rho0={1}'.format(E_0,Rho(t0)),size = 22)
        pyplot.savefig("{0}/TMuTraj.png".format(rhoDir))
        pyplot.clf()
        plot(time,rhoAr)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\rho_B (fm^{-3})$',size = 26)
        pyplot.title(r'En0={0} Rho0={1}'.format(E_0,Rho(t0)),size = 22)
        pyplot.savefig("{0}/rhoVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,tau_piOft)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\tau_\pi (fm)$',size = 26)
        pyplot.title(r'En0={0} Rho0={1}'.format(E_0,Rho(t0)),size = 22)
        pyplot.savefig("{0}/SrelaxVtime.png".format(rhoDir))
        pyplot.clf()
        plot(time,tau_PIOft)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel(r'$\tau_\Pi (fm)$',size = 26)
        pyplot.title(r'En0={0} Rho0={1}'.format(E_0,Rho(t0)),size = 22)
        pyplot.savefig("{0}/BrelaxVtime.png".format(rhoDir))
        pyplot.clf()


        savetxt("{0}/time.dat".format(rhoDir),time)
        savetxt("{0}/energy.dat".format(rhoDir),enAr)
        savetxt("{0}/shearT.dat".format(rhoDir),sh)
        savetxt("{0}/bulkT.dat".format(rhoDir),bu)
        savetxt("{0}/temperature.dat".format(rhoDir),tempOft)
        savetxt("{0}/entropy.dat".format(rhoDir),entrOft)
        savetxt("{0}/eta.dat".format(rhoDir),etOft)
        savetxt("{0}/zeta.dat".format(rhoDir),bOft)
        savetxt("{0}/pressure.dat".format(rhoDir),presrOft)
        savetxt("{0}/sound.dat".format(rhoDir),sndOft)
        savetxt("{0}/effPress.dat".format(rhoDir),effPres)
        savetxt("{0}/trPress.dat".format(rhoDir),trPres)
        savetxt("{0}/plOVpt.dat".format(rhoDir),plOVpt)
        savetxt("{0}/omega.dat".format(rhoDir),omega)
        savetxt("{0}/chi.dat".format(rhoDir),chi)
        savetxt("{0}/chem.dat".format(rhoDir),chemOft)
        savetxt("{0}/sRelax.dat".format(rhoDir),tau_piOft)
        savetxt("{0}/bRelax.dat".format(rhoDir),tau_PIOft)
        savetxt("{0}/rho.dat".format(rhoDir),rhoAr)


        toPlot = enArFM/temptoFourth
        toPlot2 = sh
        toPlot3 = bu
        NSshearUL = (4/(15*time))*tau_piOft
        NSbulkUL = (-bOft/time)*tau_PIOft
        NSshearUF = ((4*(enAr+presrOft))/(15*time))*tau_piOft
        NSbulkUF = (-(bOft*(enAr+presrOft))/time)*tau_piOft
        savetxt("{0}/NSshearUL.dat".format(rhoDir),NSshearUL)
        savetxt("{0}/NSbulkUL.dat".format(rhoDir),NSbulkUL)
        savetxt("{0}/NSshearUF.dat".format(rhoDir),NSshearUF)
        savetxt("{0}/NSbulkUF.dat".format(rhoDir),NSbulkUF)


        plot(time,chi)
        pyplot.xlabel('Time (fm)',size=20)
        plot(time,NSshearUL,color='k')
        pyplot.ylabel(r'$\pi/(\epsilon + p)$',size = 26)
        pyplot.title(r'$\chi$ 0={0} $\Omega$ 0={1}'.format(shearLabel,bulkLabel),size = 18)
        pyplot.savefig("{0}/chiVtime.png".format(rhoDir))
        pyplot.clf()


        plot(time,omega)
        pyplot.xlabel('Time (fm)',size=20)
        plot(time,NSbulkUL,color='k')
        pyplot.ylabel(r'$\Pi/(\epsilon + p)$',size = 26)
        pyplot.title(r'$\chi$ 0={0} $\Omega$ 0={1}'.format(shearLabel,bulkLabel),size = 18)
        pyplot.savefig("{0}/omegaVtime.png".format(rhoDir))
        pyplot.clf()


        plot(time,toPlot)
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel('Energy',size = 26)
        pyplot.title(r'$\chi$ 0={0} $\Omega$ 0={1}'.format(shearLabel,bulkLabel),size = 18)
        pyplot.savefig("{0}/energyVtime.png".format(rhoDir))
        pyplot.clf()
        # # show()

        plot(time,toPlot2)
        plot(time,NSshearUF,color='k')
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel('Shear Stress',size = 26)
        pyplot.title(r'$\chi$ 0={0} $\Omega$ 0={1}'.format(shearLabel,bulkLabel),size = 18)
        pyplot.savefig("{0}/shearVtime.png".format(rhoDir))
        pyplot.clf()
        # show()

        plot(time,toPlot3)
        plot(time,NSbulkUF,color='k')
        pyplot.xlabel('Time (fm)',size=20)
        pyplot.ylabel('Bulk',size = 26)
        pyplot.title(r'$\chi$ 0={0} $\Omega$ 0={1}'.format(shearLabel,bulkLabel),size = 18)
        pyplot.savefig("{0}/bulkVtime.png".format(rhoDir))
        pyplot.clf()
    else:
        print(rho_0/t0,'[',shear0,',',bulk0,']',': DID NOT make it, final chem:', chem)
