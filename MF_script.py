# second order MF for thalamic TC and RE cell populations


import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from scipy.special import erfc
from mytools import progressBar,double_gaussian,ornstein_uhlenbeck

#TAG functions

from cell_library import loadparams


def TF(typ,fexc,finh,adapt,diff=False,record=False,I=0):

    p = params[ typ ]
    P,Nexc,Ninh,Qe,Qi,Cm,El = p.P,p.Nexc,p.Ninh,p.Qe,p.Qi,p.Cm,p.El
    Gl,Te,Ti,Ee,Ei = p.Gl,p.Te,p.Ti,p.Ee,p.Ei

    if fexc<1e-9: fe=1e-9
    else: fe = fexc*Nexc
    if finh<1e-9: fi=1e-9
    else: fi = finh*Ninh


    muGi = Qi*Ti*fi;
    muGe = Qe*Te*fe;
    muG = Gl+muGe+muGi;
    muV = (muGe*Ee+muGi*Ei+Gl*El - adapt + I)/muG;
    
    
    muGn = muG/Gl;
    Tm = Cm/muG;
    
    Ue =  Qe/muG*(Ee-muV);
    Ui = Qi/muG*(Ei-muV);

    sV = np.sqrt(fe*(Ue*Te)*(Ue*Te)/2./(Te+Tm)+fi*(Ui*Ti)*(Ui*Ti)/2./(Ti+Tm));

    Tv = ( fe*(Ue*Te)*(Ue*Te) + fi*(Qi*Ui)*(Qi*Ui)) /( fe*(Ue*Te)*(Ue*Te)/(Te+Tm) + fi*(Qi*Ui)*(Qi*Ui)/(Ti+Tm) );
    TvN = Tv*Gl/Cm;
    
    muV0=-60e-3;
    DmuV0=10e-3;
    sV0=4e-3;
    DsV0=6e-3;
    TvN0=0.5;
    DTvN0=1.;

    vthre = P[0] + P[1]*(muV-muV0)/DmuV0 + P[2]*(sV-sV0)/DsV0 + P[3]*(TvN-TvN0)/DTvN0 \
        + P[4]*((muV-muV0)/DmuV0)*((muV-muV0)/DmuV0) + P[5]*(muV-muV0)/DmuV0*(sV-sV0)/DsV0 + P[6]*(muV-muV0)/DmuV0*(TvN-TvN0)/DTvN0 \
        + P[7]*((sV-sV0)/DsV0)*((sV-sV0)/DsV0) + P[8]*(sV-sV0)/DsV0*(TvN-TvN0)/DTvN0 + P[9]*((TvN-TvN0)/DTvN0)*((TvN-TvN0)/DTvN0);


    frout = 1/(2.*Tv)*erfc( (vthre-muV)/(np.sqrt(2)*sV) );
    # if frout<1e-8: frout=1e-8
    
    if record: return frout,muV,sV
    else: return frout;


def MFw(typ, adapt, ff, fexc, finh):

    p = params[ typ ]
    a,b,Tw = p.a,p.b,p.Tw
    Nexc,Ninh,Qe,Qi,El = p.Nexc,p.Ninh,p.Qe,p.Qi,p.El
    Gl,Te,Ti,Ee,Ei = p.Gl,p.Te,p.Ti,p.Ee,p.Ei 

    if fexc<1e-9: fe=1e-9
    else: fe = fexc*Nexc
    if finh<1e-9: fi=1e-9
    else: fi = finh*Ninh

    muGi = Qi*Ti*fi;
    muGe = Qe*Te*fe;
    muG = Gl+muGe+muGi;
    muV = (muGe*Ee+muGi*Ei+Gl*El - adapt)/muG;

    adapt = -adapt/Tw + b*ff + a*(muV-El)/Tw;

    return adapt;


# load transfer function fit parameters

params = loadparams('thalamus_ACh')
# params = loadparams('thalamus_control')


# simulation setup

T=5e-3

tfinal=2
dt=5e-4
df=1e-32

tsteps=int(tfinal/dt)
t=np.linspace(0, tfinal, tsteps)


#TAG-inputs
# external_input=np.zeros(tsteps)
external_input=np.full(tsteps, 10.)
# external_input[:int(50e-3/dt)]+=2
# external_input[int(1/dt):int(1.5/dt)]+=4
# external_input+=double_gaussian(t, 1.0, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.1, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.2, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.3, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.4, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.5, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.6, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.7, 0.005, 0.05, 10)
# external_input+=double_gaussian(t, 1.8, 0.005, 0.05, 10)
# ampl,freq=2,10
# external_input=ampl/2*(1-np.cos(freq*2*np.pi*t))
# noise
# external_input=ornstein_uhlenbeck(tsteps,tfinal, .1, 1/5e-2, 10.5, start=0,seed=1)

stim=np.zeros(tsteps)
# stim=np.full(tsteps, 6.)
stim[int(1/dt):int(1.5/dt)]+=40
# stim[int(1/dt):int(2/dt)]=10.
# stim=double_gaussian(t, 1, 0.002, 0.2, 20)
# amplitude,frequency,offset=10,10,0
# stim=offset+amplitude/2*(1-np.cos(2*frequency*np.pi*t))
# stim=ornstein_uhlenbeck(tsteps,tfinal, 10, 1/5e-2, 10.5, start=0,seed=1)



#TAG-initial conds
fecont=0;
ficont=30;
we=fecont*params['TC'].b*params['TC'].Tw
wi=ficont*params['RE'].b*params['RE'].Tw
# wi=-10e-12
cee,cei,cii=1,1,1


LSwe,LSwi=[],[]
LSfe,LSfi=[],[]
LScee,LScii=[],[]
test,test2=[],[]
LSmuVe,LSsVe,LSmuVi,LSsVi=[],[],[],[]

for i in progressBar(range(len(t))):

    # if i>10*int(50e-3/dt): params['RE'].Ee=5e-3
    
    fecontold=fecont
    ficontold=ficont
    weold,wiold=we,wi
    ceeold,ceiold,ciiold=cee,cei,cii
    TCfe = external_input[i]+stim[i]/8
    REfe = external_input[i]+fecont/16

    #-TFs
    Fe,muVe,sVe = TF('TC',TCfe,ficont,we,record=True)
    Fi,muVi,sVi = TF('RE',REfe,ficont,wi,record=True)


    #TAG TF derivatives

    # first order
    dveFe = 0
    dviFe = (TF('TC',TCfe,ficont+df/2,we,True)-TF('TC',TCfe,ficont-df/2,we,True))/df
    dveFi = (TF('RE',REfe+df/32,ficont,wi,True)-TF('RE',REfe-df/32,ficont,wi,True))/df
    dviFi = (TF('RE',REfe,ficont+df/2,wi,True)-TF('RE',REfe,ficont-df/2,wi,True))/df

    # second order
    dvedveFe = 0
    dvidveFe = 0
    dvidviFe = ( TF('TC',TCfe,ficont+df,we,True) - 2*TF('TC',TCfe,ficont,we,True) + TF('TC',TCfe,ficont-df,we,True) )/df**2
    dvedviFe = 0
    dvedveFi = ( TF('RE',REfe+df/16,ficont,wi,True) - 2*TF('RE',REfe,ficont,wi,True) + TF('RE',REfe-df/16,ficont,wi,True) )/df**2
    dvidveFi = ( (TF('RE',REfe+df/32,ficont+df/2,wi,True)-TF('RE',REfe-df/32,ficont+df/2,wi,True))\
                - (TF('RE',REfe+df/32,ficont-df/2,wi,True)-TF('RE',REfe-df/32,ficont-df/2,wi,True)) )/df**2
    dvidviFi = ( TF('RE',REfe,ficont+df,wi,True) - 2*TF('RE',REfe,ficont,wi,True) + TF('RE',REfe,ficont-df,wi,True) )/df**2
    dvedviFi = ( (TF('RE',REfe+df/32,ficont+df/2,wi,True)-TF('RE',REfe+df/32,ficont-df/2,wi,True))\
                - (TF('RE',REfe-df/32,ficont+df/2,wi,True)-TF('RE',REfe-df/32,ficont-df/2,wi,True)) )/df**2


    #TAG INTEGRATION

    #-first order EULER
    # fecont += dt/T*( (Fe-fecont) )
    # ficont += dt/T*( (Fi-ficont) )

    #-first order HEUN
    # fecont += dt/T*(Fe-fecont)
    # fecont = fecontold + dt/T/2*(Fe-fecontold + TF('TC',TCfe,ficont,we)-fecont)
    # ficont += dt/T*(Fi-ficont)
    # ficont = ficontold + dt/T/2*(Fi-ficontold + TF('RE',REfe,ficont,wi)-ficont)

    #-second order EULER
    fecont += dt/T*( (Fe-fecont) + (cee*dvedveFe+cei*dvedviFe+cii*dvidviFe+cei*dvidveFe)/2 )
    ficont += dt/T*( (Fi-ficont) + (cee*dvedveFi+cei*dvedviFi+cii*dvidviFi+cei*dvidveFi)/2 )

    #-second order HEUN
    # fecont += dt/T*( (Fe-fecont) + (cee*dvedveFe+cei*dvedviFe+cii*dvidviFe+cei*dvidveFe)/2 )
    # fecont = fecontold + dt/T/2*( (Fe-fecontold) + (TF('TC',TCfe,ficont,we)-fecont) + (cee*dvedveFe+cei*dvedviFe+cii*dvidviFe+cei*dvidveFe) )
    # ficont += dt/T*( (Fi-ficont) + (cee*dvedveFi+cei*dvedviFi+cii*dvidviFi+cei*dvidveFi)/2 )
    # ficton = ficontold + dt/T/2*( (Fi-ficontold) + (TF('RE',REfe,ficont,wi)-ficont) + (cee*dvedveFi+cei*dvedviFi+cii*dvidviFi+cei*dvidveFi) )


    #-adaptation EULER
    we += dt*MFw('TC',we,fecontold,TCfe,ficontold)
    wi += dt*MFw('RE',wi,ficontold,REfe,ficontold)
    # wi += dt*MFw('RE',wi,ficontold,ficontold/4,fecontold*4)

    #-adaptation HEUN
    # we += dt*MFw('TC',we,fecontold,0)
    # we = weold + dt/2*( MFw('TC',weold,fecontold,0) + MFw('TC',we,fecontold,0) )
    # wi += dt*MFw('RE',wi,REfe,ficontold)
    # wi = wiold + dt/2*( MFw('RE',wiold,REfe,ficontold) + MFw('RE',wi,REfe,ficontold) )

    if fecont<1e-9: fecont=1e-9
    if ficont<1e-9: ficont=1e-9
    if fecont>200: fecont=200
    if ficont>200: ficont=200
    # if we<0 or wi<0:
    #     print('w<0 !!')

    LSfe.append(float(fecont))
    LSfi.append(float(ficont))
    LSwe.append(float(we))
    LSwi.append(float(wi))


    #-covariances EULER
    cee += dt/T*( Fe*(1/T-Fe)/500 + (Fe-fecontold)**2 + 2*cee*dveFe + 2*ceiold*dviFe - 2*cee)
    cei += dt/T*( (Fe-fecontold)*(Fi-ficontold) + cei*dveFe + ceeold*dveFi + ciiold*dviFe + cei*dviFi - 2*cei)
    cii += dt/T*( Fi*(1/T-Fi)/500 + (Fi-ficontold)**2 + 2*cii*dviFi + 2*ceiold*dveFi - 2*cii)

    #-covariances HEUN
    # cee += dt/T*( Fe*(1/T-Fe)/500 + (Fe-fecontold)**2 + 2*cee*dveFe + 2*cei*dviFe - 2*cee)
    # cee = ceeold + dt/T*( Fe*(1/T-Fe)/500 + (Fe-fecontold)**2 + ceeold*dveFe + 2*cei*dviFe - ceeold + cee*dveFe - cee)
    # cei += dt/T*( (Fe-fecontold)*(Fi-ficontold) + cee*dveFi + cei*dveFe + cei*dviFi + cii*dviFe - 2*cei)
    # cei = ceiold + dt/T*( (Fe-fecontold)*(Fi-ficontold) + cee*dveFi + ceiold*dveFe/2 + ceiold*dviFi/2 + cii*dviFe - ceiold + cei*dveFe/2 + cei*dviFi/2 - cei)
    # cii += dt/T*( Fi*(1/T-Fi)/500 + (Fi-ficontold)**2 + 2*cei*dveFi + 2*cii*dviFi - 2*cii)
    # cii = ciiold + dt/T*( Fi*(1/T-Fi)/500 + (Fi-ficontold)**2 + 2*cei*dveFi + ciiold*dviFi - ciiold + cii*dviFi - cii)

    if cee<1e-9: cee=1e-9
    if cii<1e-9: cii=1e-9
    if cei<1e-9: cei=1e-9

    # cee=np.sqrt(cee)
    # cei=np.sqrt(cei)
    # cii=np.sqrt(cii)
    # LScee.append(cee)
    # LScii.append(cii)
    LScee.append(np.sqrt(cee))
    LScii.append(np.sqrt(cii))

    #-test
    # test.append(muV)
    # test2.append(sV)
    LSmuVe.append(muVe)
    LSsVe.append(sVe)
    LSmuVi.append(muVi)
    LSsVi.append(sVi)

print(fecont,we)
print(ficont,wi)

#-end of loop

LSfe=np.array(LSfe)
LSfi=np.array(LSfi)
LSwe=np.array(LSwe)
LSwi=np.array(LSwi)
LScee=np.array(LScee)
LScii=np.array(LScii)


#TAG SAVE

np.save('data\\MF_out', np.vstack((LSfe,LSfi)))
np.save('data\\MF_out_cov', np.vstack((LScee,LScii)))
# np.savetxt('test.txt',test)
# np.save('data\\MF_out_adaptnew', np.vstack((LSfe,LSfi)))

# np.save('data\\MF_out_w')

np.save('data\\MF_Vexc',np.vstack((LSmuVe,LSsVe)))
np.save('data\\MF_Vinh',np.vstack((LSmuVi,LSsVi)))


#TAG PLOTS

#-testplots
# plt.plot(test)
# plt.plot(test2)
# plt.show()
# plt.plot(LScee, 'b')
# plt.plot(LScii, 'r')
# plt.show()

# plt.plot(t,LSmuVe,'b')
# plt.fill_between(t, np.array(LSmuVe)-LSsVe,np.array(LSmuVe)+LSsVe,color='b',alpha=.2)
# plt.show()

#-main plot
fig = plt.figure()
fig.subplots_adjust(hspace=0.001)
gs = gridspec.GridSpec(2,1, height_ratios=[3,1])
ax3=fig.add_subplot(gs[0])
ax2=fig.add_subplot(gs[1],sharex=ax3)
ax1=ax3.twinx()

ax3.set_axis_off()
ax1.tick_params(labelright=False,labelbottom=False,labelleft=True,labeltop=False,which='both',
                left=True,right=True,bottom=False, top=False)
ax2.tick_params(which='both',right=True,top=True,grid_alpha=0.3)
ax2.tick_params(axis='y', labelsize=8, size=2,grid_alpha=0)

ax1.set_xlim(0,tfinal)
maxpoint=max(np.concatenate([LSfe,LSfi]))
ax1.set_ylim(-maxpoint/10,maxpoint+maxpoint/5)

ax1.plot(t, LSfi, c='r', label=r'$\nu_{\mathrm{RE}}$')
ax1.fill_between(t, LSfi-LScii, LSfi+LScii, color='r', label=r'$\sigma_{\mathrm{RE}}$', alpha=0.2)
ax1.plot(t, LSfe, c='b', label=r'$\nu_{\mathrm{TC}}$')
ax1.fill_between(t, LSfe-LScee, LSfe+LScee, color='b', label=r'$\sigma_{\mathrm{TC}}$', alpha=0.2)
ax1.plot(t,external_input, c='black', label=r'$P_C$')
ax1.plot(t,stim, c='black', ls='--', label=r'$P_S$')

ax2.grid()

ax2.plot(t, LSwe*1e12, c='b', label=r'$\omega_{\mathrm{TC}}$')
ax2.plot(t, LSwi*1e12, c='r', label=r'$\omega_{\mathrm{RE}}$')


ax1.yaxis.set_label_position('left')
ax1.set_ylabel(r'frequency $\nu$ [Hz]',fontsize=12)

ax2.set_xlabel(r'time $t$ [s]',fontsize=12)
ax2.set_ylabel(r'adaptation $\omega$ [pA]',fontsize=10,position=(0,0.5))

leg1 = ax1.legend(bbox_to_anchor=(1.205, 1.0), loc=1, borderaxespad=0.)
leg2 = ax2.legend(bbox_to_anchor=(1.215, 1.0), loc=1, borderaxespad=0.)
ax1.add_artist(leg1)

plt.savefig('gfx\\MF_PLOT.png', dpi=200, bbox_inches='tight')
