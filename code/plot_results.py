## 
## plot_results.py 
## 
## Renders plots of results.  
## 
## Tim Farrell, tmf@bu.edu
## QNL, BU 
## 20160415
##
import numpy as np
from sys import argv
from os import listdir
import scipy.io as scio
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

## Functions 
def plot_sv_spectrum(S, subject_run, pdf): 
   fig, ax1 = plt.subplots() 
   xs = np.arange(len(S)) 
   ys = np.array([(S[i]**2/sum(S**2)) for i in xs])
   ax1.plot(xs, ys, '-o')
   ax1.set_ylim((0, 0.04)) 
   ax2 = ax1.twinx() 
   ax2.plot(xs, np.cumsum(ys))
   ax2.set_ylim((0, 1.0))
   plt.xlabel('Singular Value')
   ax1.set_ylabel('Relative Variance')
   ax2.set_ylabel('Cumulative Relative Variance') 
   plt.title('Singular Value Spectrum ' + subject_run) 
   pdf.savefig(fig) 
   plt.close(fig) 
   return

def plot_temporal_comps(V, subject_run, pdf, first_n=5): 
   fig = plt.figure() 
   plt.plot(V[1:first_n, :].T) 
   plt.xlabel('Time')
   plt.title('Temporal Components ' + subject_run + ' first_n=' + str(first_n)) 
   pdf.savefig(fig) 
   plt.close(fig) 
   return

## Main
# args 
spectrum_only = False
if argv[1] == '--spectrum-only':
   spectrum_only = True
# dirs 
datadir = '/projectnb/bohland/HCP/data/' 
svds = listdir(datadir + '/svds') 
svs = PdfPages(datadir + 'plots/singular_value_spectrums.pdf') 
if not spectrum_only: tcs = PdfPages(datadir + 'plots/temporal_components.pdf') 
# plot for each svd
for svd in svds:
    print("Loading " + svd + "...") 
    mat = scio.loadmat(datadir + 'svds/' + svd) 
    print("Plotting...") 
    print("\tSingular value spectrum...") 
    plot_sv_spectrum(mat['S'][0], '-'.join(svd.split('-')[:2]), svs)
    if not spectrum_only: 
       print("\tTemporal components...") 
       plot_temporal_comps(mat['V'], '-'.join(svd.split('-')[:2]), tcs)
    print("Done.")
svs.close()
if not spectrum_only: tcs.close() 
print("Done.") 
