When you want to restart a job, follow the following steps:

0) Assuming you are in the run directory

1) go the ncfiles/restart_files/
2) figure out which restart files are "safe":
	a) type: ncdump -v dt_t r_xxxx_1_xx_xx_xx.nc
           (of course, choose a file that exists)
        
           the result of this command is a long text at the end of which
	   you will see something like:  
              
              dt_t = 0.005, 70 ;

           this tells you that the r_xxxx_1_xx_xx_xx.nc files have
           the data for t=70

           now type: ncdump -v dt_t r_xxxx_0_xx_xx_xx.nc

              dt_t = 0.005, 50 ;

           this tells you that the r_xxxx_0_xx_xx_xx.nc files have
           the data for t=50

        b) in the outputfiles, check the last iteration. Is it above 70?
           if it is then the time_dump 70 is safe, so you can use the '1'            restart files.
    	   if the last iteration is exactly 70, you should rather use the 	     '0' files, as the data for time_dump 70 are probably not safe

3) (optional, if you want to change the time_dumps, etc)
   a) open restart_defs_tregisters.F90 and edit the next lines as you see fit:  
  integer::div1=20
  integer::until1=190
  integer::div2=5
  integer::until2=215
  integer::div3=1
  integer::until3=245
  integer::div4=5
  integer::until4=270
  integer::div5=10
  integer::until5

   b) open sub_Lathys.slurm, and make sure that in the following line, -nhm corresponds to the last time_dump you want:
  ### on fait le calcul proprement dit
  srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./quiet_plasma -ncxyz 1500 32 32 -dt 0.005 --gstep 1 -nhm 14000 -pn earth -r 0 -ug 1

  For example, say we want to go up to tmax=110, we should write -nhm 22000

4) open restart.sh
    write the time of the last safe dump, and the corresponding restart
    files number. In the example above, this would be:
	LASTDUMP=50
	RF=0  #(0 or 1)

5) execute restart.sh (/!\ in the Code/Soumission/ directory, e.g. "source ./restart.sh", if already there. Without the right path, the version in 
   ~/bin will be executed, which will not work)


