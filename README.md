## HCP

####Description

A preprocessing pipeline for fMRI data from the Human Connectome Project (WU-Minn). 

####Directory Structure  

    /projectnb/bohland/HCP/
    
	code/
	# see each script for more detailed documentation 
	
	      hcp_analysis_utils.py 
	      # contains functions that perform preprocessing
	      
	      svd_dtseries.sh 
	      # calls svd_dtseries.py with supplied args

	      svd_dtseries.py 
	      # computes svd of a partition in $HCP/data/id_partitions.txt

	      svd_dtseries_parallel.sh 
	      # parallelizes svd_dtseries.sh across several partitions
	      # in accordance with supplied args 

	data/ 
	      err/
	      # stderr of each array job (of qsub) 
	      
	      id_partitions.txt
	      # list of subjects partitioned based on args supplied to svd_dtseries_parallel.sh
	      # where each row is a partition

	      imgs/
	      # noise mask imgs (saved as .nii.gz), if saved during svd computation
	      # resampled parcellation imgs (aparc+aseg; saved as .nii.gz) 

	      out/
	      # stdout of each array job (of qsub) 

	      plots/ 
	      # preprocessing plots, if saved during svd computation
	      # one per subject 

	      svds/
	      # computed svds, saved as .mat

