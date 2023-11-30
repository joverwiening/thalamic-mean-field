from brian2 import *
prefs.codegen.target = "numpy"
from mytools import ornstein_uhlenbeck

start_scope()
control=1

DT=0.01 # time step
defaultclock.dt = DT*ms
N_inh = 500 # number of inhibitory neurons
N_exc = 500 # number of excitatory neurons
N_ped = 8000 # external pop

TotTime=3000 #Simulation duration (ms)
duration = TotTime*ms
tsteps=int(TotTime/DT)
tt = np.linspace(0,TotTime, tsteps)

isyn = TimedArray((tt>=1500)*(tt<=2000), defaultclock.dt)
# sigma=0 #0.2 # 0.025
# isyn=TimedArray(ornstein_uhlenbeck(tsteps,TotTime/1e3, 0, 1/5e-1, sigma, start=0,seed=None,nonzero=False), defaultclock.dt)
# *nA

eqs='''
dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
dw/dt = (a*(v-El)-w)/tau_w:ampere
dGsynI/dt = -GsynI/Tsyn : siemens
dGsynE/dt = -GsynE/Tsyn : siemens
Is=A*isyn(t) - In : ampere
A:ampere
In:ampere
Cm:farad
gl:siemens
El:volt
a:siemens
b:ampere
tau_w:second
Dt:volt
Vt:volt
Ee:volt
Ei:volt
Tsyn:second
'''

# Population 1 [inhibitory] - RE - Reticular

G_inh = NeuronGroup(1, eqs, threshold='v > -20*mV', reset='v = -55*mV; w += b', refractory='5*ms', method='heun')
# init:
G_inh.v = -55.*mV
G_inh.w = 0.*pA
# synaptic parameters
G_inh.GsynI = 0.0*nS
G_inh.GsynE = 0.0*nS
G_inh.Ee = 0.*mV
G_inh.Ei = -80.*mV
G_inh.Tsyn = 5.*ms
# cell parameters
G_inh.Cm = 200.*pF
G_inh.gl = 10.*nS
G_inh.Vt = -45.*mV
G_inh.Dt = 2.5*mV
G_inh.tau_w = 200.*ms
# G_inh.Is = 0.0*nA # external input
G_inh.A = 450*pA
G_inh.In = 0*pA
G_inh.El = -75.*mV
G_inh.a = 8.0*nS
G_inh.b = 10.*pA

if control:
    G_inh.gl = 13*nS
    G_inh.El = -85*mV #-85
    G_inh.a = 28*nS # 28
    G_inh.b = 20.*pA
    G_inh.tau_w = 230.*ms
    G_inh.A = 900*pA


# Population 2 [excitatory] - TC - Thalamocortical

G_exc = NeuronGroup(1, eqs, threshold='v > -20.0*mV', reset='v = -50*mV; w += b', refractory='5*ms',  method='heun')
# init
G_exc.v = -50.*mV
G_exc.w = 0.*pA
# synaptic parameters
G_exc.GsynI = 0.0*nS
G_exc.GsynE = 0.0*nS
G_exc.Ee = 0.*mV
G_exc.Ei = -80.*mV
G_exc.Tsyn = 5.*ms
# cell parameters
G_exc.Cm = 160.*pF
G_exc.gl = 10.*nS
G_exc.Vt = -50.*mV
G_exc.Dt = 4.5*mV
G_exc.tau_w = 200.*ms
# G_exc.Is = 0.0*nA # ext inp
G_exc.A = 150*pA
G_exc.In = 0*pA
G_exc.El = -65.*mV # -55
G_exc.a = 0.*nS
G_exc.b = 10.*pA

if control:
    G_exc.gl = 9.5*nS
    G_exc.El = -70*mV # -63(dom) -> -73(me)
    G_exc.a = 24*nS
    # G_exc.a = 14*nS
    G_exc.b = 20*pA
    G_exc.tau_w = 270.*ms
    G_exc.A = 350*pA


# external drive--------------------------------------------------------------------------

PED = PoissonGroup(1, rates=0*Hz)

SED = PoissonGroup(8000, rates=0*Hz)

# Network-----------------------------------------------------------------------------

# S_ei = Synapses(G_exc, G_inh, on_pre='GsynE_post+=20*nS')
# S_ei.connect()

# S_ii = Synapses(G_inh, G_inh, on_pre='GsynI_post+=50*nS')
# S_ii.connect()

# S_ie = Synapses(G_inh, G_exc, on_pre='GsynI_post+=30*nS')
# S_ie.connect()


S_pe = Synapses(PED, G_exc, on_pre='GsynE_post+=24*nS') # hyp
# S_pe = Synapses(PED, G_exc, on_pre='GsynE_post+=14*nS') # rest
# S_pe.connect(p=0.10)
S_pe.connect()

# S_se = Synapses(SED, G_exc, on_pre='GsynE_post+=1*nS')
# S_se.connect(p=0.20)

S_pi = Synapses(PED, G_inh, on_pre='GsynE_post+=28*nS')
# S_pi.connect(p=0.05)
S_pi.connect()


# Recording tools -------------------------------------------------------------------------------

# FRG_inh = PopulationRateMonitor(G_inh)
# FRG_exc = PopulationRateMonitor(G_exc)

SM_P = SpikeMonitor(PED)
SM_exc = SpikeMonitor(G_exc)

# Useful trick to record global variables ------------------------------------------------------

Gw_inh = NeuronGroup(1, 'Wtot : ampere', method='rk4')
Gw_exc = NeuronGroup(1, 'Wtot : ampere', method='rk4')

SwInh1=Synapses(G_inh, Gw_inh, 'Wtot_post = w_pre : ampere (summed)')
SwInh1.connect(p=1)
SwExc1=Synapses(G_exc, Gw_exc, 'Wtot_post = w_pre : ampere (summed)')
SwExc1.connect(p=1)

MWinh = StateMonitor(Gw_inh, 'Wtot', record=0)
MWexc = StateMonitor(Gw_exc, 'Wtot', record=0)



GV_inh = NeuronGroup(1, 'Vtot : volt', method='rk4')
GV_exc = NeuronGroup(1, 'Vtot : volt', method='rk4')

SvInh1=Synapses(G_inh, GV_inh, 'Vtot_post = v_pre : volt (summed)')
SvInh1.connect(p=1)
SvExc1=Synapses(G_exc, GV_exc, 'Vtot_post = v_pre : volt (summed)')
SvExc1.connect(p=1)

MVinh = StateMonitor(GV_inh, 'Vtot', record=0)
MVexc = StateMonitor(GV_exc, 'Vtot', record=0)


# Run simulation -------------------------------------------------------------------------------

run(duration)

# Plots -------------------------------------------------------------------------------



# # prepare firing rate
# def bin_array(array, BIN, time_array):
#     N0 = int(BIN/(time_array[1]-time_array[0]))
#     N1 = int((time_array[-1]-time_array[0])/BIN)
#     return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

# BIN=5
# time_array = arange(int(TotTime/DT))*DT



# LfrG_exc=array(FRG_exc.rate/Hz)
# TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)

# LfrG_inh=array(FRG_inh.rate/Hz)
# TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)


rateP,rateE = np.zeros((2,tsteps))

rateP[(SM_P.t/ms/DT).astype('int')]=1
rateE[(SM_exc.t/ms/DT).astype('int')]=1


# np.save('data\\singlecell_rates_N.npy',np.vstack((rateP,rateE)))


# np.save('Wtot.npy',[np.array(MWinh.Wtot[0]/mamp),np.array(MWexc.Wtot[0]/mamp)])

# np.save('Vtot.npy',[np.array(MVinh.Vtot[0]/mV),np.array(MVexc.Vtot[0]/mV)])

# plt.subplot(211)
# plt.plot(tt,np.array(MVexc.Vtot[0]/mV),'-b')
# plt.plot(tt,np.array(MWexc.Wtot[0]/pA/10),'-b',alpha=.3)
# plt.subplot(212)
# plt.plot(tt,np.array(MVinh.Vtot[0]/mV),'-r')
# plt.plot(tt,np.array(MWinh.Wtot[0]/pA/10),'-r',alpha=.3)
# plt.show()

# ------------------------------------

plt.rcParams.update({'font.size': 15})

fig,ax = plt.subplots(2,1,sharex=True,figsize=(5,7))
ax0 = ax[0].twinx()
ax1 = ax[1].twinx()

ax[0].set_title('Awake (ACh present)\n')
# ax[0].set_title('Sleep (ACh absent)\n')

ax[0].plot(tt,np.array(MVexc.Vtot[0]/mV),'-b',lw=2,label='$v_{\mathrm{TC}}$')
ax0.plot(tt,np.array(MWexc.Wtot[0]/pA),'-b',lw=2,alpha=.4,label='$\omega_{\mathrm{TC}}$')
ax0.plot(tt,(tt>=1500)*(tt<=2000)*G_exc.A/pA,'-k',lw=3,alpha=.3)
ax[1].plot(tt,np.array(MVinh.Vtot[0]/mV),'-r',lw=2,label='$v_{\mathrm{RE}}$')
ax1.plot(tt,np.array(MWinh.Wtot[0]/pA),'-r',lw=2,alpha=.4,label='$\omega_{\mathrm{RE}}$')
ax1.plot(tt,(tt>=1500)*(tt<=2000)*G_inh.A/pA,'-k',lw=3,alpha=.3)

ax[0].set_xlim(1300,2300)
ax0.set_xlim(1300,2300)
ax[1].set_xlim(1300,2300)

ax[0].legend(loc='upper left')
ax0.legend(loc='upper right')
ax[1].legend(loc='upper left')
ax1.legend(loc='upper right')

ax[1].set_xlabel('time (s)')
# fig.supylabel('membrane potential (mV)')
fig.supylabel('ext. current $I_\mathrm{ext}$ / adaptation $\omega$ (pA)')

def align_yaxis(ax1, v1, ax2, v2):
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2,(y1-y2)/2,v2)
    adjust_yaxis(ax1,(y2-y1)/2,v1)

def adjust_yaxis(ax,ydif,v):
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny>maxy or (-miny==maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)

if not control:
    align_yaxis(ax[0], -65, ax0, 0)
    align_yaxis(ax[1], -75, ax1, 0)
elif control:
    align_yaxis(ax[0], -70, ax0, 0)
    align_yaxis(ax[1], -85, ax1, 0)

fig.tight_layout()
# plt.show()
plt.savefig('gfx\\paper_singlecells_sleep.png',dpi=300,bbox_inches='tight')