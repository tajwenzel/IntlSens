#!/bin/bash                                                                                
#SBATCH --job-name=countryrun    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)          
#SBATCH --mail-user=nwenzel@uw.edu     # Where to send mail                    
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)       
#SBATCH --ntasks=1         # We only run one R instance = 1 task               
#SBATCH --cpus-per-task=4 # number of threads we want to run on
#SBATCH --time=48:00:00               # Time limit hrs:min:sec                 
#SBATCH --output=serial_test_%j.log   # Standard output and error log                  
# Load R (version 3.3.2)                                                                   
                                            
ml R/3.6.0-foss-2016b-fh1

R CMD BATCH --no-save --no-restore FUNC_FITS4.R
