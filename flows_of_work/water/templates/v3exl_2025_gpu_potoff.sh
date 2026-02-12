{% extends "slurm.sh" %}
{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
######SBATCH --gpus-per-node=1
#SBATCH --gres=gpu#:2 --oversubscribe
#SBATCH --exclusive

{%- endif %}

{% set walltime = operations |calc_walltime(parallel) %}
{% if walltime %}
#SBATCH --time {{walltime|format_timedelta}}
{% endif %}


{% block tasks %}
#SBATCH --ntasks={{operations|calc_tasks('np',parallel, force) }}
{% endblock tasks %}

#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=reslab35ai8111
######SBATCH --exclude=reslab32ai8111
######SBATCH --exclude=ressrv4ai8111,ressrv6ai8111,reslab32ai8111
######SBATCH --nodelist=res-lab42-ai8111
######SBATCH --ntasks=1
######SBATCH --ntasks-per-core=1
######SBATCH --ntasks-per-node={{operations|calc_tasks('np',parallel, force) }}

echo  "Running on host" $HOSTNAME
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo  "Time is" $current_time
source ~/.bashrc
#mamba init
mamba activate ham_clone_with_mda

{% if gpus %}
#module load cuda/11.0
{%- endif %}
{% endblock header %}
{% block body %}
    {{- super () -}}
{% endblock body %} (edited) 

