#!/bin/zsh

### job name
#SBATCH --job-name=AD

### File / path where STDOUT will be written
### %J is the job id
#SBATCH --output=output_%J.txt

### ask for 10 processes (num_tasks = num_processes)
#SBATCH --ntasks=10

### ask for 1 node
#SBATCH --nodes=1

### memory usage asked for per MPI process (default in MB)
#SBATCH --mem-per-cpu=512M

### Use exclusive for timings (Each job only has ONE node)
#SBATCH --exclusive

### execute commands
module load gcc/11

export NAG_KUSARI_FILE=/home/eu657030/Software/dco_cpp/nag_key.txt

make run_p_multi_var_draft_new

