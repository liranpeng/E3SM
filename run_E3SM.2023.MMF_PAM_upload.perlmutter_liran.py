#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, numpy as np
from shutil import copy2
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = 'm3312'
src_dir = os.getenv('HOME')+'/repositories/E3SM' # branch => ?

# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True

debug_mode = False
noaccel = False
gcm_dt = 200

#queue = 'debug'  # regular / debug 
queue = 'regular'

arch = 'GNUCPU' # GNUCPU / GNUGPU

if queue=='debug'  : stop_opt,stop_n,resub,walltime = 'ndays',1, 0,'00:30:00'
if queue=='regular': stop_opt,stop_n,resub,walltime = 'ndays',6,0,'24:00:00'

grid = 'ne30pg2_EC30to60E2r2'; num_nodes = 32
#grid = 'ne4pg2_oQU480'; num_nodes = 1

# compset = 'F2010-MMF1' # MMF+SAM
#compset = 'F2010-MMF2' # MMF+PAM
compset = 'F2010-MMF2-ECPP' # MMF+PAM+ECPP

case = '.'.join(['E3SM','2023-MMF2_TEST_ECPP_PAM-00757',arch,grid,compset])
crm_nx = 64
crm_ny = 1
crm_dt = 5
crm_dx = 2000
if debug_mode: case += '.regular'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

if 'CPU' in arch: max_mpi_per_node,atm_nthrds  = 128,1 ; max_task_per_node = 128
if 'GPU' in arch: max_mpi_per_node,atm_nthrds  =   4,8 ; max_task_per_node = 32
atm_ntasks = max_mpi_per_node*num_nodes

if 'CPU' in arch: case_root = f'{os.getenv("PSCRATCH")}/e3sm_scratch/pm-cpu/{case}'
if 'GPU' in arch: case_root = f'{os.getenv("PSCRATCH")}/e3sm_scratch/pm-gpu/{case}'

#---------------------------------------------------------------------------------------------------
if newcase :
   if os.path.isdir(case_root): exit(f'\n{clr.RED}This case already exists!{clr.END}\n')
   cmd = f'{src_dir}/cime/scripts/create_newcase'
   cmd += f' --case {case}'
   cmd += f' --output-root {case_root} '
   cmd += f' --script-root {case_root}/case_scripts '
   cmd += f' --handle-preexisting-dirs u '
   cmd += f' --compset {compset}'
   cmd += f' --res {grid} '
   cmd += f' --project {acct} '
   cmd += f' --walltime {walltime} '
   if arch=='GNUCPU' : cmd += f' -mach pm-cpu -compiler gnu    -pecount {atm_ntasks}x{atm_nthrds} '
   if arch=='GNUGPU' : cmd += f' -mach pm-gpu -compiler gnugpu -pecount {atm_ntasks}x{atm_nthrds} '
   run_cmd(cmd)
   # # Copy this run script into the case directory
   # timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d.%H%M%S')
   # run_cmd(f'cp {os.path.realpath(__file__)} {case_dir}/{case}/run_script.{timestamp}.py')
#---------------------------------------------------------------------------------------------------
os.chdir(f'{case_root}/case_scripts')
#---------------------------------------------------------------------------------------------------
if config :
   run_cmd(f'./xmlchange EXEROOT={case_root}/bld ')
   run_cmd(f'./xmlchange RUNDIR={case_root}/run ')
   #-------------------------------------------------------
   # if specifying ncdata, do it here to avoid an error message
   if 'init_file_atm' in locals():
      file = open('user_nl_eam','w');file.write(f' ncdata = \'{init_file_atm}\' \n');file.close()
   #-------------------------------------------------------
   if 'nlev'   in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -nlev {nlev} \" ')
   if 'crm_dx' in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_nz {crm_dx} \" ')
   if 'crm_nz' in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_nz {crm_nz} \" ')
   if 'crm_nx' in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_nx {crm_nx} \" ')
   if 'crm_dt' in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_dt {crm_dt} \" ')
   if 'crm_ny' in locals(): run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_ny {crm_ny} \" ')
#-------------------------------------------------------
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')
#---------------------------------------------------------------------------------------------------
if build : 
   if debug_mode: run_cmd('./xmlchange --file env_build.xml --id DEBUG --val TRUE ')
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
if submit : 
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   #------------------------------
   # Specify history output frequency and variables
   #------------------------------
   file.write(' nhtfrq    = 0,1,-6 \n')
   file.write(' mfilt     = 1,8,4 \n')
   file.write(" fincl1 = 'Z3'") # this is for easier use of height axis on profile plots
   file.write('\n')
   file.write(" fincl2 = 'PS','TS','PSL'")
   file.write(          ",'PRECT','TMQ','SO2','so4_a3','CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'")
   file.write(          ",'PRECC','PRECL','totna500','TOT_CLD_VISTAU','TOT_ICLD_VISTAU'")
   file.write(          ",'LHFLX','SHFLX','ABSORB','EXTINCT'")             # surface fluxes
   file.write(          ",'FSNT','FLNT','FLUT','Z3'")        # Net TOM heating rates
   file.write(          ",'FLNS','FSNS','CLOUD','CLDTOT','CLDLOW','CLDHGH'")               # Surface rad for total column heating
   file.write(          ",'FSNTC','FLNTC','OMEGA500'")             # clear sky heating rates for CRE
   file.write(          ",'TGCLDLWP','TGCLDIWP'")       # liq & ice water path
   file.write(          ",'TUQ','TVQ','TMSO2','TMDMS','TMdst_a1','TMso4_a3'")                 # vapor transport for AR tracking
   # variables for tracking stuff like hurricanes
   file.write(          ",'TBOT:I','QBOT:I','UBOT:I','VBOT:I'") # lowest model leve
   file.write(          ",'T900:I','Q900:I','U900:I','V900:I'") # 900mb data
   file.write(          ",'T850:I','Q850:I','U850:I','V850:I'") # 850mb data
   file.write(          ",'Z300:I','Z500:I'")
   file.write(          ",'OMEGA850:I','OMEGA500:I'")
   file.write(          ",'U200:I','V200:I'")
   file.write('\n')
   # 3D variables
   file.write(" fincl3 = 'PS','TS','PSL'")
   file.write(          ",'T','Q','Z3'")                      # 3D thermodynamic budget components
   file.write(          ",'U','V','OMEGA'")                    # 3D velocity components
   file.write(          ",'QRL','QRS'")                        # 3D radiative heating profiles
   file.write("linoz_data_path='/global/homes/h/heroplr/inputfile/' \n")
   file.write("linoz_data_file='linoz2004_2006jpl_c081216.nc' \n") 
   # Other namelist stuff
   #------------------------------  
   # at the beginning of the case creation
   if 'gcm_dt' in locals():
     dtime = gcm_dt; ncpl = int( 86400 / dtime )

   if 'noaccel' in locals(): file.write(f' use_crm_accel = .false. \n')


# when writing the namelist
   if gcm_dt < 20*60 :
         file.write(f'dt_tracer_factor    = 1 \n')
         file.write(f'dt_remap_factor     = 1 \n')
         file.write(f'hypervis_subcycle_q = 1 \n')
         file.write(f'se_tstep = {gcm_dt} \n')

   # file.write(f' cosp_lite = .true. \n')
   if 'init_file_atm' in locals(): file.write(f' ncdata = \'{init_file_atm}\' \n')
   # file.write(" inithist = \'ENDOFRUN\' \n")
   file.close()
   #-------------------------------------------------------
   # LND namelist
   #-------------------------------------------------------
   if 'init_file_lnd' in locals() or 'data_file_lnd' in locals():
      nfile = 'user_nl_elm'
      file = open(nfile,'w')
      if 'init_file_lnd' in locals(): file.write(f' finidat = \'{init_file_lnd}\' \n')
      if 'data_file_lnd' in locals(): file.write(f' fsurdat = \'{data_file_lnd}\' \n')
      # file.write(f' check_finidat_fsurdat_consistency = .false. \n')
      file.close()
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
# when setting the run time options
   if 'ncpl' in locals(): run_cmd(f'./xmlchange ATM_NCPL={str(ncpl)}')

   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_QUEUE={queue},JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')

   if continue_run :
      run_cmd('./xmlchange CONTINUE_RUN=TRUE')   
   else:
      run_cmd('./xmlchange CONTINUE_RUN=FALSE')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   run_cmd('./case.submit')

#---------------------------------------------------------------------------------------------------
# Print the case name again
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n') 
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
