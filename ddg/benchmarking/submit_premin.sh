#!/bin/bash                        #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    #-- the shell for the job
#$ -o premin_console_out/          #-- output directory (fill in)
#$ -e premin_console_out/          #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=5G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=1G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=12:00:00                #-- runtime limit 
#$ -t 1-169                        #-- specify the number of tasks

./run_preminimize.sh ${SGE_TASK_ID}
