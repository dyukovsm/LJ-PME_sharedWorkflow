HOW TO RUN:

```bash
git clone https://github.com/dyukovsm/LJ-PME_sharedWorkflow.git
```

ensure the packages listed in LJPMECUTenv.yml file are installed, preferably by creating a new environment from the yml.

```bash
gh release download v1.0.0-trajectories
```
downloads 5Gb of trajectory data that is used as a starting point. Concatinate parts of the downalod and extract the result:

```bash
cat trajectories_full.tar.xz.part* | tar -xJvf -
```

initialize worskpace for spce water by :

```bash
python init.py
```

and run by :

```bash
python project.py submit
```

cutoff jobs tend to be faster, so init.py can be run first without PME in the cutoff settings, 
and when complete and analysis done, re-initialized with PME and run again : this works as long as new
types of statepoints are not added, but only existing types are appended. Old work will not be overwritten 
by initialize.

----------------------------------------------------------------------------------------------------------------

troubleshooting :

1 ) a common error is that signac sometimes poorly interacts with node configuration: 
```bash
python project.py submit --force
```

tends to submit jobs when error claiming that "non-gpu nodes cannot submit gpu-based jobs" shows up.

2 ) try running 1000 timesteps for both equilibrations and production, changing all saving to be every 100 steps

3 ) some grid nodes like rom does not run gmx, exiting with error : 
	subprocess.CalledProcessError: Command '['gmx', 'trjconv', '-s', 'INPUT_TEMPLATE_SLAB.tpr', '-f', 'INPUT_TEMPLATE_SLAB.trr', '-o', 'ELONGATED_BOX_PLACEHOLDER.gro', '-dump', '45.0']' died with <Signals.SIGILL:$
	
	my solution is to ask for a different node: voh works : 

	srun -q gpu --gres=gpu:1 --constraint=avx512 -N 1 -n 4 --mem 8G -t 11:00:00 --pty bash

	then try running project.py manually :
		
	python project.py run -o <job-that-breaks>
		
	doing conda pre-compiled gromacs also works but I was not gettign a better performance than 1 ns/day for LJ-PME.
