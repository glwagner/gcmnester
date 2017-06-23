#!/bin/bash

email="glwagner@mit.edu"

tasks=30
modelName="bb_z100_test"
runName="dt240"
jobName="bbTest"

partition="sched_mit_hill"
#partition="sched_mit_raffaele"
#partition="sched_any_quicktest"

# Do stuff  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Set tasksPerNode and simulation time to maximum
if [ "$partition" = "sched_mit_raffaele" ]; then
    simtime="48:00:00"
    tasksPerNode=20
elif [ "$partition" = "sched_mit_hill" ]; then
    simtime="12:00:00"
    tasksPerNode=16
elif [ "$partition" = "sched_any_quicktest" ]; then
    simtime="00:15:00"
    tasksPerNode=16
fi

# Number of nodes, rounding up to nearest integer
nodes=$(((tasks+tasksPerNode-1)/tasksPerNode))

# Crucial directories
baseDir="/home/glwagner/patches/models/$modelName"
scratchDir="/nobackup1/glwagner/patches/models/$modelName/runs/$runName"

codeDir=$baseDir/code
namelistDir=$baseDir/namelists
inputDir=$baseDir/input

# Clean and link
if [ ! -d $scratchDir ]; then
  mkdir -p $scratchDir;
  mkdir -p $scratchDir/diags;
fi

if [ ! -h jra55 ]; then
  ln -s /nobackup1/glwagner/patches/forcing/jra55 $scratchDir/jra55
fi

# Copy code, input, and executable
printf '%s' "Copying input and grid files, code, namelists, and executable... "
cp $inputDir/*.bin $scratchDir
cp $inputDir/tile*.mitgrid $scratchDir
cp -rf $codeDir/ $scratchDir
cp -f $namelistDir/* $scratchDir
cp -f $baseDir/build/mitgcmuv $scratchDir/mitgcmuv
echo "done."

# Remove run.slurm if it exists, and then reprint it.
if [ -f run.slurm ]; then 
    rm run.slurm
fi

echo "#!/bin/bash

#SBATCH --partition $partition
#SBATCH --exclude node258
#SBATCH --nodes $nodes
#SBATCH --ntasks $tasks
#SBATCH --exclusive
#SBATCH --mem 64000
#SBATCH --time $simtime
#SBATCH --error stderr
#SBATCH --output stdout
#SBATCH --job-name $jobName
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user $email

module add engaging/intel/2013.1.046
module add netcdf netcdff

ulimit -s unlimited

cd $scratchDir

date > run.MITGCM.timing

mpirun -n $tasks ./mitgcmuv

date >> run.MITGCM.timing" >> run.slurm 

# Submit the job
sbatch run.slurm
