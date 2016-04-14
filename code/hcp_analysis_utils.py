## 
## hcp_analysis_utils.py 
## 
## Various utility functions that support 
## precprocessing of HCP dtseries data. 
## 
## Tim Farrell, tmf@bu.edu
## QNL, BU 
## 20160210
## 
import numba
import numpy as np 
from math import pi
import nibabel as nb
from os import remove
from scipy.linalg import svd
import matplotlib.pyplot as plt
from subprocess import check_output
from scipy.interpolate import interp1d
from scipy.ndimage import binary_erosion
from sklearn.linear_model import LinearRegression 
from scipy.signal import detrend, butter, filtfilt
from matplotlib.backends.backend_pdf import PdfPages

# Useful dir/ filename fcns 
def subject_datadir(subject): 
   return '/projectnb/connectomedb/Q6/'+ subject +'/MNINonLinear/'

def fnf_ts(subject, run, dense): 
   # gives name for subject-run timeseries file, either dense or not 
   if dense: 
      return subject_datadir(subject) +'Results/rfMRI_'+ run +\
        '/rfMRI_'+ run +'_Atlas_hp2000_clean.dtseries.nii'
   else: 
      return subject_datadir(subject) +'Results/rfMRI_'+ run +\
        '/rfMRI_'+ run +'_hp2000_clean.nii.gz'

# Simple numba-compiled array fcns 
@numba.jit("f8(f8[:])") 
def sum_c(arr): 
   total = 0 
   for i in xrange(len(arr)): 
      total += arr[i]
   return(total)

@numba.jit("f8(f8[:])") 
def mean_c(arr): 
   return(float(sum_c(arr))/float(len(arr)))

@numba.jit("f8[:,:](f8[:,:])") 
def doubly_center_c(mat): 
   #[(xij - mean_c(mat[i,:]) - mean_c(mat[:,j]) + mean_c(mat))
   # for xij in mat] 
   nrow, ncol = mat.shape 
   # get means 
   total = 0
   rmeans = np.zeros(nrow)
   cmeans = np.zeros(ncol) 
   for r in xrange(nrow): 
      rmeans[r] = mean_c(mat[r, :]) 
      total += sum_c(mat[r, :])  # sum up rows to get total
   for c in xrange(ncol): 
      cmeans[c] = mean_c(mat[:, c]) 
   mat_mean = (float(total)/float(nrow*ncol))
   # do computation
   for r in xrange(nrow):
      for c in xrange(ncol): 
         mat[r, c] = mat[r, c] - rmeans[r] - cmeans[c] + mat_mean 
   return(mat)

def plot(plots, subject_run): 
   pdf = PdfPages('/projectnb/bohland/HCP/data/plots/preprocessing-' + subject_run + '.pdf')
   for p in plots: 
      if len(p) > 2:
         before, after, step = p 
         plot_before_after(before, after,\
                           np.random.randint(after.shape[0], size=5), step, pdf)
      else: 
         ts_mat, step = p
         plot_ts_mat(ts_mat, step, pdf) 
   print("Saving plot...") 
   pdf.close() 
   print("Done.") 
         
def plot_ts_mat(ts_mat, step, pdf): 
   # assumes timeseries arranged (ordinate, timeseries) 
   fig = plt.figure()
   if len(ts_mat.shape) == 2:  plt.plot(ts_mat[:, :].T)
   else:                       plt.plot(ts_mat) 
   plt.title(step) 
   pdf.savefig(fig) 
   plt.close(fig) 

def plot_before_after(before, after, rand_idx, step, pdf): 
   # plots two dense timeseries matrices before and after some processing step
   # assumes matrices arranged (ordinate, timeseries) 
   assert before.shape == after.shape, "Timeseries must be of the same shape." 
   fig = plt.figure()
   print("Plotting " + step + "...")
   # lines 
   plt.subplot(2, 1, 1) 
   plt.plot(before[rand_idx, :].T) 
   plt.title(step)
   plt.ylabel('Before') 
   plt.subplot(2, 1, 2) 
   plt.plot(after[rand_idx, :].T) 
   plt.ylabel('After') 
   pdf.savefig(fig)
   plt.close(fig) 
   # images 
   fig = plt.figure() 
   plt.subplot(2, 1, 1)
   plt.imshow(before, cmap='gray', extent=[-4, 4, -1.5, 1.5])
   plt.title(step)
   plt.ylabel('Before')
   #plt.colorbar()
   plt.subplot(2, 1, 2) 
   plt.imshow(after, cmap='gray', extent=[-4, 4, -1.5, 1.5]) 
   #plt.colorbar()
   plt.ylabel('After') 
   pdf.savefig(fig)
   plt.close(fig) 

# 
# For a given HCP subject (e.g. '111312') and resting state run (e.g. 'REST1_LR), gives
# back the corresponding cleaned dense timeseries. If full_cifti, gives as <Cifti2Image object> 
# (i.e. with header/etc included); else just the matrix.     
# 
def get_timeseries(subject, run, full_cifti): 
    if not full_cifti:
        return np.asarray(nb.load(fnf_ts(subject, run, dense=True)).get_data()[0][0][0][0].T,\
                             dtype=np.dtype('f8'))
    else: 
        return nb.load(fnf_ts(subject, run, dense=True))

#
# For a given HCP subject (e.g. '100307') and resting state run (e.g. 'REST1_LR'), 
# get movement regressors. 
#
def get_motion_regressors(subject, run): 
   return np.array([map(float, filter(lambda x: x != '', l.strip().split(' ')))\
                       for l in open(subject_datadir(subject) + '/Results/rfMRI_' + run\
                                        +'/Movement_Regressors.txt').readlines()],\
                      dtype=np.dtype('f8'))

def get_framewise_disp(motion, radius=50): 
   # calculates framewise displacment from motion regressors 
   motion = np.copy(motion[:, 0:6])
   motion[:, 3:6] = (float(1)/float(360)) * motion[:, 3:6] * (2*pi*radius) 
   motion_dt = [m1 - m2 for m1, m2 in zip(motion[:-1, :], motion[1:, :])] 
   fd = np.asarray([0] +\
        [np.apply_along_axis(sum, 0, np.apply_along_axis(abs, 0, d)) for d in motion_dt]) 
   return fd

# gets noise mask from the data of a segmentation image
def get_noise_mask(seg_img, noise_regions, size_limit=None, save_as=None): 
   seg_img_data = seg_img.get_data() 
   # get values that represent noise regions (e.g. white matter, etc.) 
   # using FreeSurfer color coding 
   noise_vals = check_output(["egrep '" + '|'.join(noise_regions) + "'"  +\
                              " $FREESURFER_HOME/FreeSurferColorLUT.txt" +\
                              " | awk '{print $1, $2}'"], shell=True)
   noise_vals = [float(l.split(' ')[0]) for l in noise_vals.split('\n') if not l == ''] 
   # get indices of segmentation file that contain those noise values
   # Note: noise_idx is a 3-tuple like (0_idx, 1_idx, 2_idx) 
   #       where i_idx is an index along axis i of voxels containing a noise_val
   #       and len(i_idx) == #(noise voxels) for all i
   noise_idx = np.where(np.apply_along_axis(lambda xs:\
                        map(lambda x: x in noise_vals, xs), 1, seg_img_data))
   noise_idx = np.array(zip(noise_idx[0], noise_idx[1], noise_idx[2]))
   # construct 3D matrix indicator matrix as noise mask (i.e. noise_mask[i, j] in {0,1})   
   noise_mask = np.zeros(seg_img_data.shape)
   for x, y, z in noise_idx: 
      noise_mask[x, y, z] = 1
   # erode 
   noise_mask_eroded = binary_erosion(noise_mask).astype(noise_mask.dtype)
   noise_idx = np.where(noise_mask_eroded > 0)                         # update noise_idx
   if not size_limit: 
      noise_idx = np.array(zip(noise_idx[0], noise_idx[1], noise_idx[2])) 
   else: # if size limit, sample randomly from the eroded mask
      noise_idx = np.array(zip(noise_idx[0], noise_idx[1], noise_idx[2]))[\
                                   np.random.randint(len(noise_idx[0]), size=size_limit)]
   if save_as: 
      for mask, label in [(noise_mask if size_limit else noise_mask_eroded,\
                              (str(size_limit) if size_limit else 'eroded'))]:  
         noise_mask_nii = nb.Nifti1Image(mask, seg_img.affine, header=seg_img.header) 
         noise_mask_nii.update_header() 
         nb.save(noise_mask_nii, '/projectnb/bohland/HCP/data/imgs/' + save_as + '-' +\
                 '-'.join(noise_regions) + '-' + label + '.nii.gz') 
   return noise_idx

#
# For a given HCP subject (e.g. '100307') and resting state run (e.g. 'REST1_LR), 
# gets the principal components (PCs) of the corresponding noise matrix, 
# a voxel by timeseries of CSF and white-matter (i.e. those areas that contain just noise); 
# and returns the first_n of those PCs as noise regressors.   
# 
def get_noise_regressors(subject, run, first_n=5, plots=None, save_noise_mask=False,\
                         save_resampled_seg=False, save_noise_dts=False, size_limit=None,\
                         noise_regions=['White-Matter','Ventricle']): 
   # resample segmentation file to timeseries space 
   seg_f = subject_datadir(subject) + 'aparc+aseg.nii.gz'
   ts_f = fnf_ts(subject, run, dense=False)
   resampled_seg_f = '/projectnb/bohland/HCP/data/imgs/resampled_aparc+aseg-'\
                       + subject + '-' + run + '.nii.gz'
   r = check_output(['$FREESURFER_HOME/bin/mri_vol2vol --mov ' + seg_f +\
                     ' --targ ' + ts_f + ' --o ' + resampled_seg_f +\
                     ' --regheader --interp nearest'], shell=True) 
   # load the resampled segmentation 
   seg = nb.load(resampled_seg_f) 
   ts = nb.load(ts_f) 
   # get noise mask 
   noise_mask = get_noise_mask(seg, noise_regions=noise_regions, size_limit=size_limit,\
                               save_as=(subject+'-'+run if save_noise_mask else None))
   if not save_resampled_seg:  
      remove(resampled_seg_f)
   # get noise img and dts, size limiting if desired 
   ts_img = ts.get_data()
   ix = 0 
   noise_dts = np.zeros((noise_mask.shape[0], 1200))
   for x, y, z in noise_mask: 
      noise_dts[ix, :] = np.copy(ts_img[x, y, z, :])
      ix = ix + 1
   # center dts
   noise_dts_ = doubly_center_c(np.copy(noise_dts))
   if plots:
      plots = plots + [(noise_dts, noise_dts_, 'noise centering ' + subject + '-' + run)] 
   # compute SVD 
   _, _, noise_V = svd(noise_dts_)   
   # return first_n PCs 
   if plots: 
      return (noise_V[:, 0:first_n], plots)
   else: 
      return noise_V[:, 0:first_n]

# 
# For a given HCP subject and resting-state run: 
# 
#     a) Detrends each ordinate (i.e. along axis 0). 
#
#     b) "Regresses out" noise and motion from each timepoint (i.e. 
# 	 along axis 1). Specifically, it subtracts the prediction of
# 	 a multivariate linear model of ordinate signal (at time t) as a 
# 	 function of noise and motion from each ordinate signal (at t).  
#
#     c) Censors timepoints that have framewise displacements(FDs) above
# 	 some threshold, via linear interpolation. 
#
#     d) Filters each doubly-detrended, censor-interpolated timeseries
# 	 with butterworth bandpass filter. 
# 
def get_preprocessed_ts(subject, run, fd_threshold=0.5, filter_order=2,\
                        filter_band=(0.01, 0.4), TR=0.720, filtfilt_padtype='odd', 
                        save_plots=False, save_noise_mask=False, noise_size_limit=None, 
                        include_args_in_plot=None): 
   dense_ts = get_timeseries(subject, run, full_cifti=False) 
   ## a) detrend along axis 0 
   if save_plots: 
      if include_args_in_plot: 
         import inspect
         args, _, _, vals = inspect.getargvalues(inspect.currentframe())
         assert all([(a in args) for a in include_args_in_plot]), "Args to be included in plots, " +\
                         "must match with 'get_preprocessed_ts' function args: " + str(args) + "."  
         argvals = ';'.join(map(lambda p: '='.join(p),\
                                   [(a, str(vals[a])) for a in include_args_in_plot])).replace(' ','') 
         subject_run = subject + '-' + run + '-' + argvals
      else: 
         subject_run = subject + '-' + run 
      dense_ts_dt = detrend(np.copy(dense_ts), axis=0) 
      dts_origin = np.copy(dense_ts)
      dts_detrend = np.copy(dense_ts_dt)

      plots = [(dense_ts, dense_ts_dt, 'detrending ' + subject_run)] 
      
      dense_ts = np.copy(dense_ts_dt) 
   else: 
      plots = None
      dense_ts = detrend(dense_ts, axis=0) 
   # to line up with regressors; consider transposing regressors instead
   dense_ts = dense_ts.T
   
   ## b) detrend noise and motion along axis 1
   # get regressors
   motion = get_motion_regressors(subject, run)
   if save_plots:  noise, plots = get_noise_regressors(subject, run, save_noise_mask=save_noise_mask,\
                                                          plots=plots, size_limit=noise_size_limit) 
   else:           noise        = get_noise_regressors(subject, run, save_noise_mask=save_noise_mask,\
                                                          size_limit=noise_size_limit)
   regressors = np.concatenate((noise, motion), axis=1)
   if save_plots:  
      plots = plots + [(noise.T, 'noise regressors ' + subject_run),\
                       (motion.T, 'motion regressors ' + subject_run)] 
   # build linear model 
   regression = LinearRegression().fit(regressors, dense_ts)
   # subtract out model, leaving residuals 
   for time_pt in xrange(regressors.shape[0]): 
      dense_ts[time_pt, :] = dense_ts[time_pt, :] - regression.predict(regressors[time_pt, :]) 
   
   if save_plots: 
      plots = plots + [(dts_detrend, dense_ts.T, 'regress ' + subject_run)]
      dts_regressed = np.copy(dense_ts.T)
   
   ## c) censor timepoints with FD > fd_threshold
   # get FDs from motion_regressors
   framewise_disp = get_framewise_disp(motion)
   if save_plots: 
      plots = plots + [(framewise_disp, 'framewise displacement ' + subject_run)] 
   # get good_idx (i.e. where FD < fd_threshold)   
   good_idx = [i for i in xrange(len(framewise_disp)) if framewise_disp[i] < fd_threshold] 
   # linearly interpolate over not good_idx   
   ## and d) filter
   nyquist_freq = (float(1)/float(TR))/float(2)                 # get nyquist
   filter_band = [float(f)/nyquist_freq for f in filter_band]   # convert filter band
   fB, fA = butter(filter_order, Wn=filter_band, btype='band')  
   for vox in xrange(dense_ts.shape[1]): 
      dense_ts[:, vox] = interp1d(good_idx, dense_ts[good_idx, vox], bounds_error=False,\
                                     fill_value=np.mean(dense_ts[:, vox]))(range(len(framewise_disp)))
      dense_ts[:, vox] = filtfilt(fB, fA, dense_ts[:, vox], padtype=filtfilt_padtype) + np.mean(dense_ts[:, vox]) 

   if save_plots:  
      plots = plots + [(dts_regressed, dense_ts.T, 'censoring/interpolation + filtering ' + subject_run),\
                       (dts_origin, dense_ts.T, 'whole preprocessing ' + subject_run)] 
      plot(plots, subject_run) 
   return dense_ts.T 

# Unit Tests
def main(): 
   '''for subject, run in [(s, r) for s in ['100307', '111312'] for r in ['REST1_LR']]: 
      print("\nGetting dtseries for subject " + subject + " and run " + run + "...")
      d = get_timeseries(subject, run, full_cifti=False) 
      print("Done. Shape of matrix is " + str(d.shape) + ".") 
      
      print("\nTesting centering...")
      d_centered = doubly_center_c(d)
      print("Shape of centered matrix " + str(d_centered.shape) + "...")
      row_means, col_means = d_centered.mean(0), d_centered.mean(1) 
      epsilon = 10**(-5) 
      if row_means[np.where(row_means > epsilon)].any()\
      or col_means[np.where(col_means > epsilon)].any(): 
         print("Centered matrix means were not within " + str(epsilon) + " of 0.")
      else: 
         print("Matrix was centered properly.")

      print("\nTesting dask svd...") 
      U, S, V = [op.compute() for op in da.linalg.svd(\
                 da.from_array(d_centered, chunks=(4000, 1200)))]
      print("Shapes of U, S, V " + ','.join([str(U.shape), str(S.shape), str(V.shape)]) + "...")
      print("Calculating U * S * V...") 
      d_svd = np.dot(np.dot(U, np.diag(S)), V)
      diff = d_centered - d_svd
      if diff[np.where(diff > epsilon)].any(): 
         print("Matrix gotten back from SVD was not within " + str(epsilon) + " of original.")
      else: 
         print("SVD worked properly.")
   '''
   return
      
if __name__ == '__main__': 
   main() 
