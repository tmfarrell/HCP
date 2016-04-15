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
   ax2 = ax1.twinx() 
   ax2.plot(xs, np.cumsum(ys)) 
   plt.xlabel('Singular Value')
   ax1.set_ylabel('Relative Variance')
   ax2.set_ylabel('Cumulative Relative Variance') 
   plt.title('Singular Value Spectrum ' + subject_run) 
   pdf.savefig(fig) 
   plt.close(fig) 
   return

def plot_temporal_comps(V, subject_run, pdf, first_n=5): 
   fig = plt.figure() 
   xs = np.arange(V.shape[0])
   for ix in range(first_n):  
       ys = V[:, ix].T
       plt.plot(xs, ys)
   plt.xlabel('Time')
   #plt.ylabel('') 
   plt.title('Temporal Components ' + subject_run + ' first_n=' + str(first_n)) 
   pdf.savefig(fig) 
   plt.close(fig) 
   return

## Main
datadir = '/projectnb/bohland/HCP/data/' 
svds = listdir(datadir + 'svds/') 
svs = PdfPages(datadir + 'plots/singular_value_spectrums.pdf') 
tcs = PdfPages(datadir + 'plots/temporal_components.pdf') 
for svd in svds:
    print("Loading " + svd + "...") 
    mat = scio.loadmat(datadir + 'svds/' + svd) 
    print("Plotting...") 
    print("\tSingular value spectrum...") 
    plot_sv_spectrum(mat['S'][0], '-'.join(svd.split('-')[:2]), svs)
    print("\tTemporal components...") 
    plot_temporal_comps(mat['V'], '-'.join(svd.split('-')[:2]), tcs)
    print("Done.")
pdf.close() 
print("Done.") 
