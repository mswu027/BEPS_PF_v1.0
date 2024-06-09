#!/bin/sh
#===========================================================
#
#    file:       beps_build_nml.sh
#
#    purpose:    shell script to generate 'beps.stdin' namelist file
#                required to run BEPS at site-level.
#                This script is initially based on file 'beps_build_control'
#                which is part of the original BEPS package by Mousong Wu.
#
#    created:    2020/02
#
#    author:     Michael Vossbeck (The Inversion Lab)
#
#=============================


# setting namelists for beps run
input_dir=""
output_dir=""
ic_yyyymmdd="20040101"
ndays=5844
nmlfile=""
#-- default settings for original configuration with two sites
#   (might be changed with option --npts)
npts=1
beps_lai_site_prefix="US_MOZ_LAI_V20240513_P1"
site_bound_prefix="US-MOz_boundary"
meteo_site_flnm_prefix="US-MOz_meteo_2004_2019"

#-- get fully qualified name of *this* script
script_name=$(readlink -f $0)
script_name=$0

#=============================
#
#--        c o m m a n d   l i n e
#
while [ $# -gt 0 ]
do
  arg="$1"
  case $arg in
    -iptdir)
      shift
      input_dir=$1
      ;;
    -outdir)
      shift
      output_dir=$1
      ;;
    -i_yyyymmdd)
      shift
      ic_yyyymmdd=$1
      ;;
    -ndays)
      shift
      ndays=$1
      ;;
    -npts)
      shift
      npts=$1
      ;;
    -*)
      echo "  ERROR::ignore unexpected option ->>>${arg}<<<-"
      ;;
    *)
      echo "  ERROR::unexpected argument ->>>${arg}<<<-"
      exit 1
      ;;
  esac
  shift
done


#-- namelist file to be generated
#   ('beps.stdin' is hard-coded in BEPS source code ('controlInput_mod.F90')
if [ "${nmlfile}" = "" ]
then
  nmlfile=beps.stdin_${ic_yyyymmdd}_ndays-${ndays}
fi

if [ ${npts} != 2 ]
then
  site_bound_prefix="US-MOz_boundary"
  beps_lai_site_prefix="US_MOZ_LAI_V20240513_P1"
  meteo_site_flnm_prefix="US-MOz_meteo_2004_2019"
fi

if [ "${output_dir}" = "" ]
then
  echo "${script_name}::INFO: -output_dir was *not* given on command line"
  exit 1
elif [ "${input_dir}" = "" ]
then
  echo "${script_name}::INFO: -iptdir was *not* given on command line"
  exit 1
else
  #-- dump settings
  echo "${script_name}::INFO: nmlfile       *****${nmlfile}*****"
  echo "${script_name}::INFO: input_dir     *****${input_dir}*****"
  echo "${script_name}::INFO: output_dir    *****${output_dir}*****"
  echo "${script_name}::INFO: init_date     *****${ic_yyyymmdd}*****"
  echo "${script_name}::INFO: ndays         *****${ndays}*****"
  echo "${script_name}::INFO: npts          *****${npts}*****"
  echo
  echo "${script_name}::INFO: beps_lai_site_prefix       *****${beps_lai_site_prefix}*****"
  echo "${script_name}::INFO: site_bound_prefix          *****${site_bound_prefix}*****"
  echo "${script_name}::INFO: meteo_site_flnm_prefix     *****${meteo_site_flnm_prefix}*****"
fi

(
cat <<EOF
  &NLS
   nlat           = 180
   nlon           = 360
   nscale         = 1            !! 0 for global simulation, 1 for site simulation
   calendar       = "GREGORIAN"  !! Valid in "GREGORIAN" and "NO_LEAP"
   icdate         = ${ic_yyyymmdd}
   icsec          = 0            !! 0 s
   sim_type       = 0            !! 0=>"initial",1=>"restart",2=>"branch"
   sim_duration   = -${ndays}    !! >=0 : simulation steps; <0 : days for simulation
   nhtfrq         = -1           !! <0  : hours for average,=0 : monthly
   restart_frq    = -${ndays}    !! >0  : simulation steps; =0 : monthly ;<0 : days for simulation
   beps_domain        = "$input_dir/beps_domain.nc"
   meteo_input        = 1        !! <0 : daily input, >=0 : hourly input
   meteo_path         = "$input_dir/meteo/beps2000-2015"
   meteo_flnm_prefix  = "era_interim_beps.1x1_"
   meteo_site_flnm_prefix  = "${meteo_site_flnm_prefix}"
   surface_data_path  = "$input_dir/beps_surface_data.nc"
   beps_yrdata_path   = "$input_dir/beps_yrdata.nc"
   n_site             = ${npts}
   beps_site_path     = "$input_dir/beps_site/"
   site_bound_prefix  = "${site_bound_prefix}"
   lai_input          = 1        !! <0 : lai is simulated, >=0 : lai as forcing
   beps_lai_path      = "$input_dir/lai/"
   beps_lai_prefix    = "beps_lai_"
   beps_lai_site_prefix = "${beps_lai_site_prefix}"
   beps_Vcmax_path    = "$input_dir/Vcmax/"
   beps_Vcmax_site_path = "$input_dir/beps_site/"
   prior_para_prefix = "Prior_Para_20240607"  !! parameter file for forward modelling
   PF_prior_para_prefix = "US-MOz_ParaBoundary_20240129"   !! parameter file for particle filter assimilation
   MC_prior_para_prefix = "US-MOz_PF_Para_20231215"   !! parameter file for Monte Carlo 
   beps_cpools        = "$input_dir/cpools/cpools_2010.nc"
   beps_PF_obs_path = "$input_dir/beps_site/"
   prior_PF_obs_prefix = "US-MOz_obs_20240129" !! observed for particle filter assimilation
   beps_MC_obs_path = "$input_dir/beps_site/"
   prior_MC_obs_prefix = "US-MOz_PF_obs_20231210"
   beps_out_dir       = "$output_dir/"
   beps_rst_dir       = "$output_dir/"
EOF
) > ${nmlfile} 

#-- link to BEPS standard name
ln -fs ${nmlfile} beps.stdin
echo "${script_name}::INFO: generated ***$(readlink -f ${nmlfile})***"
