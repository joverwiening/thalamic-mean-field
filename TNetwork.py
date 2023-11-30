from brian2 import *
from mytools import double_gaussian,ornstein_uhlenbeck
prefs.codegen.target = "numpy"

# import sys
# if not sys.warnoptions:
#     import warnings
#     warnings.simplefilter("ignore")

start_scope()

control=0
oscpop=0

DT=0.1 # time step
defaultclock.dt = DT*ms
N_inh = 500 # number of inhibitory neurons
N_exc = 500 # number of excitatory neurons

TotTime=2000 # Simulation duration (ms)
duration = TotTime*ms
tsteps=int(TotTime/DT)
tt = np.linspace(0,TotTime, tsteps)


# Equations ----------------------------------------------------------------------------------
eqs='''
dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
dw/dt = (a*(v-El)-w)/tau_w:ampere
dGsynI/dt = -GsynI/Tsyn : siemens
dGsynE/dt = -GsynE/Tsyn : siemens
Is:ampere
Cm:farad
gl:siemens
El:volt
a:siemens
tau_w:second
Dt:volt
Vt:volt
Ee:volt
Ei:volt
Tsyn:second
Vr:volt
b:ampere
'''

# Populations----------------------------------------------------------------------------------

if not oscpop:
    # Population 1 [inhibitory] - RE - Reticular

    G_inh = NeuronGroup(N_inh, eqs, threshold='v > -10*mV', reset='v = Vr; w += b', refractory='2.5*ms', method='heun')
    # init:
    G_inh.v = -75.*mV
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
    G_inh.Vt = -45.*mV # -45
    G_inh.Dt = 2.5*mV
    G_inh.tau_w = 200.*ms
    G_inh.Is = 0.*nA # external input
    G_inh.El = -75.*mV
    G_inh.a = 8.0*nS
    G_inh.b = 10.*pA
    G_inh.Vr = -55*mV

    if control:
        G_inh.gl = 13*nS
        G_inh.El = -85*mV #-85
        G_inh.a = 28*nS # 28
        G_inh.b = 20*pA
        G_inh.Vr = -48*mV
        G_inh.tau_w = 230.*ms
    
    # domenico ACh+
    G_inh.Cm = 0.08*nF
    G_inh.gl = 0.01*uS
    G_inh.El = -85*mV
    G_inh.Dt = .5*mV
    G_inh.Vt = -60*mV
    G_inh.Vr = -65*mV
    G_inh.a = 0.01*uS
    G_inh.b = 0.01*nA

    # Population 2 [excitatory] - TC - Thalamocortical

    G_exc = NeuronGroup(N_exc, eqs, threshold='v > -10.0*mV', reset='v = Vr; w += b', refractory='2.5*ms',  method='heun')
    # init
    G_exc.v = -65.*mV
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
    G_exc.Vt = -50.*mV # -50
    G_exc.Dt = 4.5*mV
    G_exc.tau_w = 200.*ms
    G_exc.Is = 0*pA # ext inp
    # G_exc.Is = -50*pA # ext inp
    G_exc.El = -65.*mV # -55
    G_exc.a = 0.*nS
    G_exc.b = 10.*pA
    G_exc.Vr = -50*mV

    # G_exc.gl[0] = 9.5*nS
    # G_exc.El[0] = -70*mV
    # # G_exc.a[0] = 24.*nS
    # G_exc.b[0] = 300.*pA
    # G_exc.tau_w[0] = 270.*ms

    # G_exc.gl = 9.5*nS
    # G_exc.El = -70*mV
    # G_exc.a = 10.*nS
    # G_exc.b = 300.*pA
    # G_exc.tau_w = 270.*ms

    if control:
        # G_exc.gl = 9.5*nS
        # G_exc.El = -70*mV # -63(dom) -> -73(me)
        # # G_exc.El = -68*mV # -63(dom) -> -73(me)
        # # G_exc.a = 24*nS # 14
        # G_exc.a = 24*nS # 14
        # # G_exc.Vr = -48*mV
        # G_exc.b = 20*pA
        # G_exc.tau_w = 270.*ms

        G_exc.gl = 9.5*nS
        G_exc.El = -70*mV
        G_exc.a = 10.*nS
        G_exc.b = 300.*pA
        G_exc.tau_w = 270.*ms

    # domenico ACh+
    G_exc.El = -55*mV
    G_exc.b = 0.01*nA


# Populations oscillations ----------------------------------------------------------------------------------
if oscpop:

    # Population 1 [inhibitory] - RE - Reticular

    b_inh = 20.*pA
    G_inh = NeuronGroup(N_inh, eqs, threshold='v > -20*mV', reset='v = -41*mV; w += b_inh', refractory='5*ms', method='heun')
    # init:
    G_inh.v = -70.*mV
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
    G_inh.Vt = -50.*mV # -45
    G_inh.Dt = 2.5*mV
    G_inh.tau_w = 270.*ms
    G_inh.Is = 0*pA # external input
    G_inh.El = -70.*mV
    G_inh.a = 28.0*nS # Spindles
    # G_inh.a = 22.0*nS # Spindles
    # G_inh.a = 10.0*nS # Delta


    # Population 2 [excitatory] - TC - Thalamocortical

    b_exc = 20.*pA
    G_exc = NeuronGroup(N_exc, eqs, threshold='v > -20.0*mV', reset='v = -48*mV; w += b_exc', refractory='5*ms',  method='heun')
    # init
    G_exc.v = -55.*mV
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
    G_exc.Vt = -50.*mV # -50
    G_exc.Dt = 2.5*mV
    G_exc.tau_w = 270.*ms
    G_exc.Is = 0*pA # ext inp
    G_exc.El = -55.*mV # -63 !!!!!!!!!!!!!!!!!!!!! -55
    G_exc.a = 26.*nS #26


#TAG external drive-----------------------------------------------------------------------
ext_inp = 2*Hz
P_ed = PoissonGroup(8000, rates=ext_inp)

# var_P = TimedArray([2*Hz,2*Hz,30*Hz,2*Hz], duration/4)
# amplitude,frequency,offset=50,2,0
# var_P = TimedArray(ornstein_uhlenbeck(tsteps, TotTime*1e-3, 10, 2, 10,start=0,seed=2)*Hz, defaultclock.dt)
# vararray=np.zeros(tsteps)
# vararray[:int(50/DT)]=100
# var_P = TimedArray(vararray*Hz, defaultclock.dt)
# P_ed=PoissonGroup(8000, rates='var_P(t)')

# amplitude,frequency,offset=50,2,0
# var_S = TimedArray((offset+amplitude/2*(1-np.cos(2*frequency*np.pi*tt/1e3)))*Hz, defaultclock.dt)
# var_S = TimedArray(double_gaussian(tt, 1e3, 0.002e3, 0.02e3, 6)*Hz
#                    +double_gaussian(tt, 1.1e3, 0.002e3, 0.02e3, 6)*Hz
#                    +double_gaussian(tt, 1.2e3, 0.002e3, 0.02e3, 6)*Hz
#                    +double_gaussian(tt, 1.3e3, 0.002e3, 0.02e3, 6)*Hz
#                    +double_gaussian(tt, 1.4e3, 0.002e3, 0.02e3, 6)*Hz
#                    +double_gaussian(tt, 1.5e3, 0.002e3, 0.02e3, 6)*Hz
#                 +0.5*Hz, defaultclock.dt)
# P_ed=PoissonGroup(8000, rates='var_S(t)')
# STIM_ed=PoissonGroup(500, rates='var_S(t)')

# var_STIM = TimedArray([20*Hz,20*Hz,20*Hz,20*Hz], duration/4)
# stim=double_gaussian(tt, 1e3, 0.002e3, 0.2e3, 20)
# var_STIM = TimedArray(stim*Hz, defaultclock.dt)
STIM_ed=PoissonGroup(500, rates=0*Hz)
# STIM_ed=PoissonGroup(50, rates=10*Hz)
# STIM_ed=PoissonGroup(500, rates='var_STIM(t)')


# Network-----------------------------------------------------------------------------

# quantal increment in synaptic conductances:
# from P_ed to G_exc (p -> e)
scale=1

Qpe = 1*nS*scale
Ppe = .10/scale

Qpi = 2*nS*scale
Ppi = .05/scale

Qei = 2*nS*scale
Pei = .10/scale

Qii = 2*nS*scale
Pii = .30/scale

Qie = 4*nS*scale
Pie = .05/scale

# probability of connection
# prbC= 0.01
# TC-RE delay
delTH= 0*ms

# create synapses
S_12 = Synapses(G_inh, G_exc, 's=1:1',on_pre='GsynI_post+=Qie', delay=delTH)
# S_12.connect('i!=j',p=Pie)
S_12.connect(p=Pie)

S_11 = Synapses(G_inh, G_inh, 's=1:1',on_pre='GsynI_post+=Qii')
# S_11.connect('i!=j',p=Pii)
S_11.connect(p=Pii)

S_21 = Synapses(G_exc, G_inh, 's=1:1',on_pre='GsynE_post+=Qei', delay=delTH)
# S_21.connect('i!=j',p=Pei)
S_21.connect(p=Pei)

S_ed_in = Synapses(P_ed, G_inh, on_pre='GsynE_post+=Qpi')
# S_ed_in = Synapses(P_ed, G_inh[:int(N_inh/20)], on_pre='GsynE_post+=Qpi')
S_ed_in.connect(p=Ppi)

S_ed_ex = Synapses(P_ed, G_exc, on_pre='GsynE_post+=Qpe')
# S_ed_ex = Synapses(P_ed, G_exc[:int(N_exc/20)], on_pre='GsynE_post+=Qpe')
S_ed_ex.connect(p=Ppe)

# Synapses(P_ed, G_exc[0],on_pre='GsynE_post+=Qpe').connect(p=2*Ppe)

# S_st_ex0 = Synapses(STIM_ed, G_exc[0], on_pre='GsynE_post+=Qpe')
# S_st_ex0.connect(p=.5)
S_st_ex = Synapses(STIM_ed, G_exc, on_pre='GsynE_post+=Qpe')
S_st_ex.connect(p=0.2) # 0.2
# S_st_ex.connect(p=0.4) # 0.2


#########OSCILLATIONS########################
# from igraph import Graph
# import igraph as ig

W12 = np.full((len(G_inh), len(G_exc)), 0)
W12[S_12.i[:], S_12.j[:]] = S_12.s[:]
W11 = np.full((len(G_inh), len(G_inh)), 0)
W11[S_11.i[:], S_11.j[:]] = S_11.s[:]
W1=np.concatenate((W11,W12), axis=1)
# W1=np.concatenate((np.zeros((500,500)),W12), axis=1)
# W1=np.concatenate((np.zeros((25,25)),W12[:25,:25]), axis=1)

W21 = np.full((len(G_exc), len(G_inh)), 0)
W21[S_21.i[:], S_21.j[:]] = S_21.s[:]
W2=np.concatenate((W21,np.zeros((500,500))), axis=1)
# W2=np.concatenate((W21[:25,:25],np.zeros((25,25))), axis=1)

W=np.concatenate((W1,W2), axis=0)/1
# print(W.shape)
# print(W[:5,:5])

eigenval, _ = np.linalg.eig(W)
print(np.amax(eigenval))
# print(eigenval.size)
# plt.plot(eigenval.real)
# plt.show()

# # print(W.shape)
# G = Graph.Adjacency(W)
# # m=G.motifs_randesu(size=4)
# # print(f'4-loop: {m[128]}')



# # C=G.community_leading_eigenvector()
# # print('mod=', C.modularity)
# # print('size=', len(C))

# C=G.community_walktrap(steps=10)
# print(len(C.as_clustering()))


# fig, ax = plt.subplots()
# # C=G.connected_components()
# ig.plot(G, target=ax)
# plt.show()

###########################################

# Recording tools -------------------------------------------------------------------------------

M1G_inh = SpikeMonitor(G_inh)
FRG_inh = PopulationRateMonitor(G_inh)
M1G_exc = SpikeMonitor(G_exc)
FRG_exc = PopulationRateMonitor(G_exc)
FRG_0 = PopulationRateMonitor(G_exc[0])

# FRG_ed = PopulationRateMonitor(P_ed)
FRG_ed = PopulationRateMonitor(STIM_ed)
FRG_stim = PopulationRateMonitor(STIM_ed)


MVexc = StateMonitor(G_exc, 'v', record=True)
MVinh = StateMonitor(G_inh, 'v', record=True)

MV0 = StateMonitor(G_exc[0], 'v', record=0)
MV1 = StateMonitor(G_exc[1], 'v', record=0)
MV2 = StateMonitor(G_exc[2], 'v', record=0)
MV3 = StateMonitor(G_exc[3], 'v', record=0)
MV4 = StateMonitor(G_exc[4], 'v', record=0)


# Useful trick to record global variables ------------------------------------------------------

# Gw_inh = NeuronGroup(1, 'Wtot : ampere', method='rk4')
# Gw_exc = NeuronGroup(1, 'Wtot : ampere', method='rk4')

# SwInh1=Synapses(G_inh, Gw_inh, 'Wtot_post = w_pre : ampere (summed)')
# SwInh1.connect(p=1)
# SwExc1=Synapses(G_exc, Gw_exc, 'Wtot_post = w_pre : ampere (summed)')
# SwExc1.connect(p=1)

# MWinh = StateMonitor(Gw_inh, 'Wtot', record=0)
# MWexc = StateMonitor(Gw_exc, 'Wtot', record=0)



# GV_inh = NeuronGroup(1, 'Vtot : volt', method='rk4')
# GV_exc = NeuronGroup(1, 'Vtot : volt', method='rk4')

# SvInh1=Synapses(G_inh, GV_inh, 'Vtot_post = v_pre : volt (summed)')
# SvInh1.connect(p=1)
# SvExc1=Synapses(G_exc, GV_exc, 'Vtot_post = v_pre : volt (summed)')
# SvExc1.connect(p=1)

# MVinh = StateMonitor(GV_inh, 'Vtot', record=0)
# MVexc = StateMonitor(GV_exc, 'Vtot', record=0)


# Run simulation -------------------------------------------------------------------------------

print('--##Start simulation##--')
run(duration)
print('--##End simulation##--')

# Plots -------------------------------------------------------------------------------


# prepare raster plot
RasG_inh = array([M1G_inh.t/ms, [i+N_exc for i in M1G_inh.i]])
RasG_exc = array([M1G_exc.t/ms, M1G_exc.i])
# print(RasG_exc[0]/DT)


# binning and seting a frequency that is time dependant
def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)

BIN = 5 # Size of the time windows in ms
time_array = arange(int(TotTime/DT))*DT



LfrG_exc=array(FRG_exc.rate/Hz)
TimBinned,popRateG_exc=bin_array(time_array, BIN, time_array),bin_array(LfrG_exc, BIN, time_array)
LfrG_0=array(FRG_0.rate/Hz)
TimBinned_0,popRateG_0=bin_array(time_array, 50, time_array),bin_array(LfrG_0, 50, time_array)

LfrG_inh=array(FRG_inh.rate/Hz)
TimBinned,popRateG_inh=bin_array(time_array, BIN, time_array),bin_array(LfrG_inh, BIN, time_array)

LfrG_ed=array(FRG_ed.rate/Hz)
TimBinned,popRateG_ed=bin_array(time_array, BIN, time_array),bin_array(LfrG_ed, BIN, time_array)
LfrG_stim=array(FRG_stim.rate/Hz)
TimBinned,popRateG_stim=bin_array(time_array, BIN, time_array),bin_array(LfrG_stim, BIN, time_array)

# meanRate_inh, meanRate_exc = np.mean(popRateG_inh[int(len(popRateG_inh)/2):]), np.mean(popRateG_exc[int(len(popRateG_exc)/2):])
meanRate_inh, meanRate_exc = np.mean(popRateG_inh), np.mean(popRateG_exc)


# plt.hist(popRateG_exc[150::])
# plt.show()
np.save('data\\popRateG_exc',popRateG_exc)
np.save('data\\popRateG_inh',popRateG_inh)

# prepare membrane potential stuff
MPexc=MVexc.v/mV
# MPexc[RasG_exc[1].astype(int),(RasG_exc[0]/DT).astype(int)]=np.nan
# print(MPexc[RasG_exc[1].astype(int),(RasG_exc[0]/DT).astype(int)])
# MPexc=np.zeros(Nexc,tsteps)
# for i in range(N_exc):
#     for t in range(tsteps):
#         if (i in RasG_exc[1]) and (t in RasG_exc[0]/DT):
#             MPexc[i,t:t+5]=None

MPinh=MVinh.v/mV
# MPinh[RasG_inh[1].astype(int)-N_exc,(RasG_inh[0]/DT).astype(int)]=np.nan
# MPexc=np.zeros(Ninh,tsteps)
# for i in range(N_inh):
#     for t in range(tsteps):
#         if (i in RasG_inh[1]) and (t in RasG_inh[0]/DT):
#             MPinh[i,t:t+5]=None


# print(MVexc.v.shape)
# plt.hist(np.mean(MVexc.v/mV,axis=1))
# plt.hist(np.mean(MVinh.v/mV,axis=1))
# # plt.plot(MVexc.v[0]/mV)
# plt.show()
np.save('data\\MVexc',MPexc)
np.save('data\\MVinh',MPinh)
# np.save('data\\MVexc',MVexc.Vtot/mV)
# np.save('data\\MVinh',MVinh.Vtot/mV)

# oscillations-------------------------
# data=popRateG_exc
# norm=data-data.mean()

# size=2**np.ceil(np.log2(2*len(data) - 1)).astype('int')
# fft=np.fft.fft(norm,size)
# power=np.abs(fft)**2
# acorr=np.fft.ifft(power).real/np.var(data)/len(data)
# plt.plot(np.fft.fftshift(acorr))


# fourier------------------------
# fft=np.fft.fft(norm)
# freq=np.fft.fftfreq(len(data), 1/(DT*2e3))
# plt.plot(freq,np.abs(fft.real))
# plt.xlim(0,10)

# plt.show()
# np.save('data\\fft_SN',np.vstack((freq,fft.real)))

#TAG SAVE-------------------------------------------------------------

np.save('data\\TNetwork_out', np.vstack((popRateG_exc,popRateG_inh)))

# np.save('data\\TN_single_cell_tonic',np.vstack((TimBinned_0,popRateG_0)))

# create the figure

# plt.plot(np.array(MVinh.v[0]/mV),c='r')
plt.figure(figsize=(6,5))
plt.subplot(511)
plt.plot(tt,np.array(MV0.v[0]/mV),c='b')
plt.gca().get_xaxis().set_visible(False)
plt.subplot(512)
plt.plot(tt,np.array(MV1.v[0]/mV),c='b')
plt.gca().get_xaxis().set_visible(False)
plt.subplot(513)
plt.plot(tt,np.array(MV2.v[0]/mV),c='b')
plt.gca().get_xaxis().set_visible(False)
plt.ylabel('membrane potential (mV)')
plt.subplot(514)
plt.plot(tt,np.array(MV3.v[0]/mV),c='b')
plt.gca().get_xaxis().set_visible(False)
plt.subplot(515)
plt.plot(tt,np.array(MV4.v[0]/mV),c='b')
plt.xlabel('time (ms)')
plt.savefig('gfx\\TN_singleMPs.png',dpi=200,bbox_inches='tight')

# plt.plot(np.array(MVinh.Vtot[0]/mV),c='r')
# plt.plot(np.array(MVexc.Vtot[0]/mV),c='b')
# plt.plot(np.array(MWinh.Wtot[0]/mamp),c='r')
# plt.plot(np.array(MWexc.Wtot[0]/mamp),c='b')
# plt.show()

# ----

fig=figure(figsize=(8,12))
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)

ax1.plot(RasG_inh[0], RasG_inh[1], ',r')
ax1.plot(RasG_exc[0], RasG_exc[1], ',b')

ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Neuron index')

ax2.plot(TimBinned,popRateG_inh, 'r',label='mean firing rate RE')
# ax2.axhline(meanRate_inh, c='r',ls='--', label=f'mean inh: {meanRate_inh:.2f}')
ax2.plot(TimBinned,popRateG_exc, 'b',label='mean firing rate TC')
# ax2.axhline(meanRate_exc, c='b',ls='--', label=f'mean exc: {meanRate_exc:.2f}')

ax2.plot(TimBinned,popRateG_ed, c='black', label='ext. drive')


# ax2.set_title(f'mean inh: {meanRate_inh:.2f} | mean exc: {meanRate_exc:.2f}')
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Firing Rate (Hz)')
plt.legend()

name_fig='gfx\\TNetwork_PLOT.png'
plt.savefig(name_fig,dpi=200,bbox_inches='tight')


# name_rates='FR_2pop_Reg_ext_'+str(ext_inp)+'.npy'
# np.save(name_rates,np.array([BIN, TimBinned, popRateG_inh,popRateG_exc,LfrG_inh,LfrG_exc], dtype=object))

# np.save('Wtot.npy',[np.array(MWinh.Wtot[0]/mamp),np.array(MWexc.Wtot[0]/mamp)])

# np.save('Vtot.npy',[np.array(MVinh.Vtot[0]/mV),np.array(MVexc.Vtot[0]/mV)])

# plt.show()


#==============================

# plt.rcParams.update({'font.size': 15})
# fig=figure(figsize=(6,4.5))

# plt.plot(RasG_inh[0], RasG_inh[1], '.r',ms=1.5,alpha=1)
# plt.plot(RasG_exc[0], RasG_exc[1], '.b',ms=1.5,alpha=1)
# plt.plot(0,0, 'or',ms=4,label='RE')
# plt.plot(0,0, 'ob',ms=4,label='TC')

# plt.xlabel('time (ms)')
# plt.ylabel('neuron index')

# plt.legend(loc='upper right')
# plt.xlim(500,1800)
# plt.ylim(0,1000)

# name_fig='gfx\\paper_test.png'
# plt.savefig(name_fig,dpi=300,bbox_inches='tight')