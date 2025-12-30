# LJ-PME Shared Workflow

## Setup

```bash
git clone https://github.com/dyukovsm/LJ-PME_sharedWorkflow.git
cd LJ-PME_sharedWorkflow
git checkout wat_met_noTess
```

ensure compatible mamba environments, preferably by installing fresh from LJPMECUTenv.yml

```bash
mamba env create -f LJPMECUTenv.yml
mamba activate LJPMECUTenv
```

## Running the Workflow

The `flows_of_work` directory contains separate subdirectories for methane and water workflows.
Both run the same way, with the major difference being that water requires pre-equilibrated frames.

### Water Setup (Additional Step: water only)

Download and extract initial frames:
```bash
cd flows_of_work/water/
gh release download v1.0.0-trajectories
cat trajectories_full.tar.xz.part* | tar -xJvf -
```

### Initialize and Run (both methane and water)
```bash
python init.py
python project.py submit
```

### Output

Analysis results are written to:
- `aggregate_dens_Data_fromMIN-MAX.txt` — gas and liquid coexistence densities
- `aggregate_general_Data.txt` — all other properties of interest

Cutoff jobs run faster than PME. You can run cutoff first by excluding PME from init.py settings, 
then re-initialize with PME after cutoff analysis completes. Old work is not overwritten when only 
appending existing statepoint types.

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
```bash
srun -q gpu --gres=gpu:1 --constraint=avx512 -N 1 -n 4 --mem 8G -t 11:00:00 --pty bash
```
then try running project.py manually :
```bash		
python project.py run -o `<job-that-breaks`>
```
		
doing conda pre-compiled gromacs also works but I was not gettign a better performance than 1 ns/day for LJ-PME.
