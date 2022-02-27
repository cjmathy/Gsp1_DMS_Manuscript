#!/bin/bash                        #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
# Run 2021-04-20, chris.mathy@ucsf.edu
#$ -S /bin/bash                    #-- the shell for the job
#$ -o ddg_out/                     #-- output directory (fill in)
#$ -e ddg_out/                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=5G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=1G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=4:00:00                 #-- runtime limit 
#$ -t 1-4140                       #-- specify the number of tasks

./run_ddg.sh ${SGE_TASK_ID}
