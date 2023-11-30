import numpy as np
from mytools import AttrDict


P = {}

def loadparams(scenario):
    return P[scenario]

# load thalamus fitting params
PTC=np.load('NEW6params_TC.npy')
PRE=np.load('NEW6params_RE.npy')

# THALAMUS =============================================


# domenico with ACh ------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800, #400
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    'El' : -65e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 0,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -75e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 8e-9,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_ACh'] = params


# domenico control ------------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800, # 400
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    # 'El' : -63e-3,
    'El' : -70e-3, #-63
    'Gl' : 9.5e-9,
    'Tw' : 270e-3,
    # 'Tw' : 500e-3,
    # 'a' : 14e-9,
    'a' : 24e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -85e-3,
    'Gl' : 13e-9,
    'Tw' : 230e-3,
    'a' : 28e-9,
    # 'a' : 18e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_control'] = params


# TC only sleep ------------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800, # 400
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    # 'El' : -63e-3,
    'El' : -70e-3, #-63
    'Gl' : 9.5e-9,
    'Tw' : 270e-3,
    # 'Tw' : 500e-3,
    # 'a' : 14e-9,
    'a' : 24e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -75e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 8e-9,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_TCsleep'] = params


#  oscillations (sleep like/control) ------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800,
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    'El' : -50e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 14e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -75e-3,
    'Gl' : 10e-9,
    'Tw' : 100e-3,
    'a' : 30e-9,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_oscillations'] = params


# domenico + destexhe? osc ------------------------------------------

params = {}

params['TC'] = AttrDict({
    # 'P' : np.load('data\\FITparams_TC_osc3.npy'),
    # 'P':PTC,
    # 'P' : np.load('data\\FITparams_TC_osc_new.npy'),
    'Nexc' : 800, # 400
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    'El' : -55e-3, #-63 !!!!!!!!!!
    'Gl' : 10e-9,
    'Tw' : 270e-3,
    'a' : 28e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    # 'P' : np.load('data\\FITparams_RE_osc3.npy'),
    # 'P':PRE,
    # 'P' : np.load('data\\FITparams_RE_osc_new.npy'),
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -70e-3,
    'Gl' : 10e-9,
    'Tw' : 270e-3,
    'a' : 28e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_spindles'] = params


# destexhe 2009 osc ------------------------------------------

params = {}

params['TC'] = AttrDict({
    # 'P' : np.load('data\\FITparams_TC_osc3.npy'),
    # 'P':PTC,
    # 'Nexc' : 320, # 400
    # 'Ninh' : 20,
    'Nexc' : 640, # '*2' included from network!!
    'Ninh' : 40,
    'Qe' : 1e-9,
    'Qi' : 10e-9,
    'Cm' : 160e-12,
    'El' : -55e-3, #-63 !
    'Gl' : 10e-9,
    'Tw' : 270e-3,
    'a' : 24e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    # 'P' : np.load('data\\FITparams_RE_osc3.npy'),
    # 'P':PRE,
    # 'Nexc' : 320,
    # 'Ninh' : 80,
    'Nexc' : 640, # times 2 !!
    'Ninh' : 160,
    'Qe' : 2e-9,
    'Qi' : 0.5e-9,
    'Cm' : 200e-12,
    'El' : -70e-3,
    'Gl' : 10e-9,
    'Tw' : 230e-3,
    'a' : 28e-9,
    'b' : 20e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_destexhe'] = params


# Wolfart inputs ------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800,
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    'El' : -65e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 0,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -75e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 8e-9,
    'b' : 10e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_Wolfart'] = params


# no adaptation ------------------------------------

params = {}

params['TC'] = AttrDict({
    'P' : PTC,
    'Nexc' : 800,
    'Ninh' : 25,
    'Qe' : 1e-9,
    'Qi' : 6e-9,
    'Cm' : 160e-12,
    'El' : -65e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 0,
    'b' : 0e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
params['RE'] = AttrDict({
    'P' : PRE,
    'Nexc' : 400,
    'Ninh' : 150,
    'Qe' : 4e-9,
    'Qi' : 1e-9,
    'Cm' : 200e-12,
    'El' : -75e-3,
    'Gl' : 10e-9,
    'Tw' : 200e-3,
    'a' : 0e-9,
    'b' : 0e-12,
    'Ti' : 5e-3,
    'Te' : 5e-3,
    'Ee' : 0,
    'Ei' : -80e-3
})
P['thalamus_noadapt'] = params

