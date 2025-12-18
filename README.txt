HOW TO RUN:

ensure the backages listed in yml file are installed, preferably by creating a new environment from the yml.
ensure coordinates dir is in files dir.



initialize worskpace for spce water by :

python init.py

and run by :

python project.py submit

cutoff jobs tend to be faster, so init.py can be run first without PME in the cutoff settings, 
and when complete and analysis done, re-initialized with PME and run again : this works as long as new
types of statepoints are not added, but only existing types are appended. Old work will not be overwritten 
by initialize.

----------------------------------------------------------------------------------------------------------------

troubleshooting :

1 ) a common error is that signac sometimes poorly interacts with node configuration: 

python project.py submit --force

tends to submit jobs when error claiming that "non-gpu nodes cannot submit gpu-based jobs" shows up.

2 ) try running 1000 timesteps for both equilibrations and production, changing all saving to be every 100 steps
