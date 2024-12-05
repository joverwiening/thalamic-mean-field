# Thalamic-mean-field

The numerical codes for simulating a thalamic spiking network and mean-field as used in the paper 'A Multi-scale study of Thalamic Responsiveness' are represented here. Two populations of thalamocortical relay neurons (TC) and thalamic reticular neurons (RE) are simulated. Two paramater sets are initiated corresponding to the 'awake' and 'sleep' states defined via the concentration of the neuromodulator Acetycholine (ACh) in the thalamus (see paper for more detail).

The files for simulating spiking network and mean-field respectively are: __TNetwork.py__, __MF_sctipt.py__.

The file for fitting the transfer function F on (simulated) single cell data is: __fitTF.ipynb__, with __ExpTF.ipynb__ the single cell scan simulation for the fit data.
