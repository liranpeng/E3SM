#pragma once

#include "pam_coupler.h"
#include "Dycore.h"
#include "sat.h"
// #include "ecppvars.h"
// =========================================================================|ecpp_crm_init|=========
// allocate and initialize variables
inline void ecpp_crm_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  printf("Liran check start ECPP init\n");
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto &dm_host     = coupler.get_data_manager_host_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");  // Note that nz   = crm_nz
  auto nx           = coupler.get_option<int>("crm_nx");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");
  auto itavg1       = coupler.get_option<int>("ecpp_itavg1");
  auto itavg2       = coupler.get_option<int>("ecpp_itavg2");
  auto ntavg1       = coupler.get_option<int>("ecpp_ntavg1");
  auto ntavg2       = coupler.get_option<int>("ecpp_ntavg2");
  printf("Liran check start ECPP init 2\n");
  auto gcm_dt       = coupler.get_option<real>("gcm_dt");
  auto crm_dt       = coupler.get_option<real>("crm_dt");
  
  int kbase, ktop;
  int m;
  int nup, ndn, icrm;
  //int ntavg1, ntavg2;
  std::string msg;

  int nxstag = nx + 1;
  int nystag = ny + 1;
  int nzstag = nz + 1;

  int mode_updnthresh = 16;
/*
    1 = method originally implemented by Bill G
        wup_thresh   =  wup_stddev * abs(upthresh);
        wdown_thresh = -wdown_stddev * abs(downthresh);
    
    2 = similar to 1, but include the mean wup and wdown
        wup_thresh   = wup_bar   + wup_stddev * abs(upthresh);
        wdown_thresh = wdown_bar - wdown_stddev * abs(downthresh);
    
    3 = user specifies an absolute threshold
        wup_thresh   =  abs(upthresh);
        wdown_thresh = -abs(downthresh);
    
    4 = similar to 1, but do
        wup_thresh   =  wup_rms * abs(upthresh);
        wdown_thresh = -wdown_rms * abs(downthresh);

    5     = see description  in module_ecpp_stats.cpp
    6,  7 = see descriptions in module_ecpp_stats.cpp
    8,  9 = see descriptions in module_ecpp_stats.cpp
    10, 11 = see descriptions in module_ecpp_stats.cpp
    12, 13 = see descriptions in module_ecpp_stats.cpp
*/

  real upthresh = 1.0;
  real downthresh = 1.0;
  real upthresh2 = 0.5;
  real downthresh2 = 0.5;
  real cloudthresh = 1e-6;
  real prcpthresh = 1e-6;

  real cloudthresh_trans = 1e-5;
  real precthresh_trans = 1e-4;

  int areaavgtype = 1;
  int plumetype = 1;
  bool allcomb = false;

  int nupdraft = 0;
  int ndndraft = 0;
  int ndraft_max = 0;

  int nupdraft_max = 0;  // Note: This variable is initialized but not used in the given code
  int ndndraft_max = 0;  // Note: This variable is initialized but not used in the given code

  // variables should inside ecppvars.h
  int DN1 = 0; // !First index of downward classes
  int NCLASS_TR = 0; // !Num. of transport classes
  int ncc_in       = 2; // Nnumber of clear/cloudy sub-calsses
  int nprcp_in     = 2; // Number of non-precipitating/precipitating sub-classes.
  int NCLASS_CL = ncc_in; // Number of cloud classes
  int NCLASS_PR = nprcp_in; // Number of precipitaion classes
   
 
// Sanity check... <not included NEED TO ADD LATER!!!!>
// Line 182 to line 210

/* Determine number of updrafts and downdrafts

 Updraft kbase & ktop definition:
   ww(i,j,k) > wup_thresh for k=kbase+1 to ktop
   ww(i,j,k) <= wup_thresh at k=kbase and k=ktop+1
 These identify the "T-points" which enclose the updraft "W-points"
 and are affected by the subgrid transport of this updraft

 Downdraft kbase & ktop definition:
   ww(i,j,k) < wdown_thresh for k=kbase+1 to ktop
   ww(i,j,k) >= wdown_thresh at k=kbase and k=ktop+1
 These identify the "T-points" which enclose the downdraft "W-points"
 and are affected by the subgrid transport of this downdraft

 For both updrafts and downdrafts:
   1 <= kbase < ktop < nzstag
*/

nupdraft = 1;
ndndraft = 1;

nupdraft_max = std::max(nupdraft_max, nupdraft);
ndndraft_max = std::max(ndndraft_max, ndndraft);

DN1 = nupdraft + 2;  // Setup index of first downdraft class
NCLASS_TR = nupdraft + ndndraft + 1;

ndraft_max = 1 + nupdraft_max + ndndraft_max;

printf("Liran check start ECPP:\n");
printf("\nValue of DN1: %d: ", DN1);
printf("\nValue of NCLASS_TR: %d: ", NCLASS_TR);
printf("\nValue of ndraft_max: %d: ", ndraft_max);
printf("\nValue of plumetype: %d: ", plumetype);
printf("\nValue of nzstag: %d: ", nzstag);
printf("\nValue of NCLASS_CL: %d: ", NCLASS_CL);
printf("\nValue of NCLASS_PR: %d: ", NCLASS_PR);
printf("\nValue of nens: %d: ", nens);
dm_device.register_and_allocate<int>("crm_cnt",     "number of crm timestep count",  {nens},{"nens"});
dm_device.register_and_allocate<real>("updraftbase", "<description>", {nens}, {"nens"});
dm_device.register_and_allocate<real>("updrafttop", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("dndrafttop", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("dndraftbase", "<description>", {nens}, {"nens"});
auto crm_cnt      = dm_device.get<int,1>("crm_cnt");
auto updraftbase  = dm_device.get<real,1>("updraftbase");
auto updrafttop   = dm_device.get<real,1>("updrafttop");
auto dndrafttop   = dm_device.get<real,1>("dndrafttop");
auto dndraftbase  = dm_device.get<real,1>("dndraftbase");

itavg1 = 0;
itavg2 = 0;
ntavg1 = 0;
ntavg2 = 0;
printf("%s %.2f\n", "\nLiran check ecpp_itavg1 init:", itavg1);
printf("%s %.2f\n", "\nLiran check ecpp_itavg2 init:", itavg2);
printf("\nLiran check start updraftbase calculation\n");
for (int icrm = 0; icrm < nens; ++icrm) {
    crm_cnt (icrm) = 0;
    updraftbase(icrm) = 0;
    updrafttop(icrm) = nz - 1;
    dndrafttop(icrm) = nz - 1;
    dndraftbase(icrm) = 0;
}
 printf("Liran check start ECPP allocation here\n");
// 4D vector allocations

dm_device.register_and_allocate<real>("qlsink_bf" , "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("prain"     , "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qcloud_bf" , "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qcloudsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qcloud_bfsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qrainsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qicesum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qsnowsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qgraupsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qlsinksum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("precrsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("precsolidsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("precallsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("altsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("rhsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("cf3dsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("ecppwwsum1", "<description>", {nzstag,ny,nx,nens}, {"nzstag","y","x","nens"});  //nzstag
dm_device.register_and_allocate<real>("ecppwwsqsum1", "<description>", {nzstag,ny,nx,nens}, {"nzstag","y","x","nens"});  //nzstag
dm_device.register_and_allocate<real>("tkesgssum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qlsink_bfsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("prainsum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("qvssum1", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("liran_test4d", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
dm_device.register_and_allocate<real>("liran_test4davg", "<description>", {nz,ny,nx,nens}, {"z","y","x","nens"});
// 2D vectors
dm_device.register_and_allocate<real>("xkhvsum", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wwqui_cen_sum", "<description>", {nz,nens}, {"z","nens"});\
dm_device.register_and_allocate<real>("liran_test2d", "<description>", {nz,nens}, {"z","nens"});
//dm_device.register_and_allocate<real>("wwqui_bnd_sum", "<description>", {nz,nens}, {"z","nens"}); //nz+1
dm_device.register_and_allocate<real>("wwqui_cloudy_cen_sum", "<description>", {nz,nens}, {"z","nens"});
//dm_device.register_and_allocate<real>("wwqui_cloudy_bnd_sum", "<description>", {nz,nens}, {"z","nens"});//nz+1
//dm_device.register_and_allocate<real>("wup_thresh", "<description>", {nz,nens}, {"z","nens"});//nz+1
//dm_device.register_and_allocate<real>("wdown_thresh", "<description>", {nz,nens}, {"z","nens"});//nz+1
printf("Liran check ECPP allocation done the first part\n");

auto qlsink_bf           = dm_device.get<real,4>("qlsink_bf");
auto prain               = dm_device.get<real,4>("prain");
auto qcloud_bf           = dm_device.get<real,4>("qcloud_bf");
auto qcloudsum1          = dm_device.get<real,4>("qcloudsum1");
auto qcloud_bfsum1       = dm_device.get<real,4>("qcloud_bfsum1");
auto qrainsum1           = dm_device.get<real,4>("qrainsum1");
auto qicesum1            = dm_device.get<real,4>("qicesum1");
auto qsnowsum1           = dm_device.get<real,4>("qsnowsum1");
auto qgraupsum1          = dm_device.get<real,4>("qgraupsum1");
auto qlsinksum1          = dm_device.get<real,4>("qlsinksum1");
auto precrsum1           = dm_device.get<real,4>("precrsum1");
auto precsolidsum1       = dm_device.get<real,4>("precsolidsum1");
auto precallsum1         = dm_device.get<real,4>("precallsum1");
auto altsum1             = dm_device.get<real,4>("altsum1");
auto rhsum1              = dm_device.get<real,4>("rhsum1");
auto cf3dsum1            = dm_device.get<real,4>("cf3dsum1");
auto ecppwwsum1          = dm_device.get<real,4>("ecppwwsum1");
auto ecppwwsqsum1        = dm_device.get<real,4>("ecppwwsqsum1");
auto tkesgssum1          = dm_device.get<real,4>("tkesgssum1");
auto qlsink_bfsum1       = dm_device.get<real,4>("qlsink_bfsum1");
auto prainsum1           = dm_device.get<real,4>("prainsum1");
auto qvssum1             = dm_device.get<real,4>("qvssum1");
auto liran_test4d        = dm_device.get<real,4>("liran_test4d");
auto liran_test4davg     = dm_device.get<real,4>("liran_test4davg");

auto xkhvsum              = dm_device.get<real,2>("xkhvsum");
auto wwqui_cen_sum        = dm_device.get<real,2>("wwqui_cen_sum");
auto liran_test2d         = dm_device.get<real,2>("liran_test2d");
//auto wwqui_bnd_sum        = dm_device.get<real,2>("wwqui_bnd_sum");            // Note: Adjust dimensionality if needed
auto wwqui_cloudy_cen_sum = dm_device.get<real,2>("wwqui_cloudy_cen_sum");
//auto wwqui_cloudy_bnd_sum = dm_device.get<real,2>("wwqui_cloudy_bnd_sum");     // Note: Adjust dimensionality if needed
//auto wup_thresh           = dm_device.get<real,2>("wup_thresh");               // Note: Adjust dimensionality if needed
//auto wdown_thresh         = dm_device.get<real,2>("wdown_thresh");             // Note: Adjust dimensionality if needed

printf("Liran check start ECPP init 0 here\n");
// Initialization of 4D variables
parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int iz, int iy, int ix, int iens) {
  qlsink_bf(iz,iy,ix,iens)          = 0;
  prain(iz,iy,ix,iens)              = 0;
  qcloud_bf(iz,iy,ix,iens)          = 0;
  qcloudsum1(iz,iy,ix,iens)         = 0;
  qcloud_bfsum1(iz,iy,ix,iens)      = 0;
  qrainsum1(iz,iy,ix,iens)          = 0;
  qicesum1(iz,iy,ix,iens)           = 0;
  qsnowsum1(iz,iy,ix,iens)          = 0;
  qgraupsum1(iz,iy,ix,iens)         = 0;
  qlsinksum1(iz,iy,ix,iens)         = 0;
  precrsum1(iz,iy,ix,iens)          = 0;
  precsolidsum1(iz,iy,ix,iens)      = 0;
  precallsum1(iz,iy,ix,iens)        = 0;
  altsum1(iz,iy,ix,iens)            = 0;
  rhsum1(iz,iy,ix,iens)             = 0;
  cf3dsum1(iz,iy,ix,iens)           = 0;
  ecppwwsum1(iz,iy,ix,iens)         = 0;
  ecppwwsqsum1(iz,iy,ix,iens)       = 0;
  tkesgssum1(iz,iy,ix,iens)         = 0;
  qlsink_bfsum1(iz,iy,ix,iens)      = 0;
  prainsum1(iz,iy,ix,iens)          = 0;
  qvssum1(iz,iy,ix,iens)            = 0;
  liran_test4d(iz,iy,ix,iens)       = 1;
  liran_test4davg(iz,iy,ix,iens)    = 0;
});
printf("Liran check start ECPP init 1 here\n");
// Initialization of 2D variables
parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int iz, int iens) {
  xkhvsum(iz,iens)                 = 0;
  wwqui_cen_sum(iz,iens)           = 0;
  liran_test2d(iz,iens)            = 1;
  //wwqui_bnd_sum(iz,iens)           = 0;  // Note: This might need adjustment for nz+1 dimension
  wwqui_cloudy_cen_sum(iz,iens)    = 0;
  //wwqui_cloudy_bnd_sum(iz,iens)    = 0;  // Note: This might need adjustment for nz+1 dimension
  //wup_thresh(iz,iens)              = 0;  // Note: This might need adjustment for nz+1 dimension
  //wdown_thresh(iz,iens)            = 0;  // Note: This might need adjustment for nz+1 dimension
});

printf("Liran check start ECPP init end here\n");

}


// =========================================================================|ecpp_crm_stat|=========
// allocate and initialize variables
inline void ecpp_crm_stat( pam::PamCoupler &coupler , int nstep) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");  // Note that nz   = crm_nz
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  auto itavg1     = coupler.get_option<int>("ecpp_itavg1");
  auto itavg2     = coupler.get_option<int>("ecpp_itavg2");
  auto ntavg1     = coupler.get_option<int>("ecpp_ntavg1");
  auto ntavg2     = coupler.get_option<int>("ecpp_ntavg2");
  auto gcm_dt     = coupler.get_option<real>("gcm_dt");
  auto crm_dt     = coupler.get_option<real>("crm_dt");
  printf("%s %.2f\n", "Liran check gcm_dt:",gcm_dt);
  printf("%s %.2f\n", "Liran check crm_dt:", crm_dt);
  //------------------------------------------------------------------------------------------------
  // get variables allocated in ecpp_crm_init
  auto crm_cnt      = dm_device.get<int,1>("crm_cnt");
  auto qcloudsum1 = dm_device.get<real, 4>("qcloudsum1");
  auto qcloud_bfsum1 = dm_device.get<real, 4>("qcloud_bfsum1");
  auto qrainsum1 = dm_device.get<real, 4>("qrainsum1");
  auto qicesum1 = dm_device.get<real, 4>("qicesum1");
  auto qsnowsum1 = dm_device.get<real, 4>("qsnowsum1");
  auto qgraupsum1 = dm_device.get<real, 4>("qgraupsum1");
  auto qlsinksum1 = dm_device.get<real, 4>("qlsinksum1");
  auto precrsum1 = dm_device.get<real, 4>("precrsum1");
  auto precsolidsum1 = dm_device.get<real, 4>("precsolidsum1");
  auto precallsum1 = dm_device.get<real, 4>("precallsum1");
  auto altsum1 = dm_device.get<real, 4>("altsum1");
  auto rhsum1 = dm_device.get<real, 4>("rhsum1");
  auto cf3dsum1 = dm_device.get<real, 4>("cf3dsum1");
  auto ecppwwsum1 = dm_device.get<real, 4>("ecppwwsum1");
  auto ecppwwsqsum1 = dm_device.get<real, 4>("ecppwwsqsum1");
  auto tkesgssum1 = dm_device.get<real, 4>("tkesgssum1");
  auto qlsink_bfsum1 = dm_device.get<real, 4>("qlsink_bfsum1");
  auto prainsum1 = dm_device.get<real, 4>("prainsum1");
  auto qvssum1 = dm_device.get<real, 4>("qvssum1");
  auto liran_test4d = dm_device.get<real, 4>("liran_test4d");
  auto liran_test4davg = dm_device.get<real, 4>("liran_test4davg");

  // Get values from PAM cloud fields
  auto host_state_shoc_tk       = dm_host.get<real,4>("state_shoc_tk");
  auto host_state_shoc_tkh      = dm_host.get<real,4>("state_shoc_tkh");
  auto host_state_qv            = dm_host.get<real,4>("state_qv");
  auto host_state_qc            = dm_host.get<real,4>("state_qc");
  auto host_state_qr            = dm_host.get<real,4>("state_qr");
  auto host_state_qi            = dm_host.get<real,4>("state_qi");
  auto state_qv                 = dm_host.get<real const,4>("state_qv").createDeviceCopy();
  auto state_qc                 = dm_host.get<real const,4>("state_qc").createDeviceCopy();
  auto state_qr                 = dm_host.get<real const,4>("state_qr").createDeviceCopy();
  auto state_qi                 = dm_host.get<real const,4>("state_qi").createDeviceCopy();
  auto crm_temp                 = dm_device.get<real,4>("temp");
  auto qvloud                   = dm_device.get<real,4>("water_vapor");
  auto qcloud                   = dm_device.get<real,4>("cloud_water");
  auto qrloud                   = dm_device.get<real,4>("rain");
  auto qiloud                   = dm_device.get<real,4>("ice");
  auto qirloud                  = dm_device.get<real,4>("ice_rime");
  //------------------------------------------------------------------------------------------------
  // Define variables used by subroutine categorization_stats
  real4d cloudmixr("cloudmixr",nz,ny,nx,nens);
  real4d cloudmixr_total("cloudmixr_total",nz,ny,nx,nens);
  real4d precmixr_total("precmixr_total",nz,ny,nx,nens);
  //real4d rhoair("rhoair",nz+1); //layer-averaged air density
  // mhwang
  // high thresholds are used to classify transport classes (following Xu et al., 2002, Q.J.R.M.S.
  real cloudthresh_trans = 1e-5; //Cloud mixing ratio beyond which cell is "cloudy" to classify transport classes (kg/kg)   +++mhwang
  // the maxium of cloudthres_trans and 0.01*qvs is used to classify transport class
  real precthresh_trans  = 1e-4; //Preciptation mixing ratio beyond which cell is raining to classify transport classes (kg/kg)  !+++mwhang

  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    crm_cnt(iens) = crm_cnt(iens) + 1;
  });
  // Some how if I move line 358 to 365 above before dm_device.get calls, the value ntavg1 will change. 
  real ecpp_ntavg1_ss = std::min(600.0, gcm_dt); // lesser of 10 minutes or the GCM timestep
  real ecpp_ntavg2_ss = gcm_dt;               // level-2 averaging period is GCM timestep
  printf("Liran check start ecpp_crm_stat 00\n");
  // Ensure ntavg2_ss is a multiple of ntavg1_ss
  ecpp_ntavg1_ss = (int)(ecpp_ntavg2_ss / (ecpp_ntavg2_ss / ecpp_ntavg1_ss));
  ntavg1 = (int)(ecpp_ntavg1_ss / crm_dt);
  ntavg2 = (int)(ecpp_ntavg2_ss / crm_dt);
  printf("%s %.2f\n", "Liran check ecpp_ntavg1_ss 00:", ecpp_ntavg1_ss);
  printf("%s %.2f\n", "Liran check crm_dt 00:", crm_dt);
  printf("%s %.2f\n", "Liran check ecpp_ntavg1_ss / crm_dt)", (ecpp_ntavg1_ss / crm_dt));
  printf("%s %.2f\n", "Liran check int(ecpp_ntavg1_ss / crm_dt))", (int)(ecpp_ntavg1_ss / crm_dt));
  //printf("%s %d %d\n", "Liran check ecpp_itavg1 1:", itavg1, ntavg1);
  //printf("%s %d %d\n", "Liran check ecpp_itavg2 1:", itavg2, ntavg2);

  // Set level-1 and level-2 averaging periods for ECPP
  
  // Calculate number of steps assuming dt evenly divides ntavg[12]_ss
/* 
!------------------------------------------------------------------------
! Main code section...
!------------------------------------------------------------------------
*/
// We could either use itavg1 or crm_cnt. Walter: The later is better for using parallel_for. 
  itavg1 = nstep+1; // The nstep begins from 0
  itavg2 = nstep+1;

  double T_test = 283.14;
  double esat_test = 0;
  double polysvp(double T, int TYPE);


// Increment the 3-D running sums for averaging period 1.
// Increments 3-D running sums for the variables averaged every
// ntavg1_mm minutes.  

printf("Liran check start ECPP ecpp_crm_stat 01\n");
esat_test = esatw_crm(T_test);
printf("%s %.2f\n", "Liran check evp:", esat_test);
parallel_for( "update sums",SimpleBounds<4>(nz, ny, nx, nens),
  YAKL_LAMBDA (int k, int j, int i, int icrm) {
    yakl::atomicAdd(liran_test4d(k,j,i,icrm) , liran_test4d(k,j,i,icrm));
    yakl::atomicAdd(qcloudsum1(k,j,i,icrm) , qcloud(k,j,i,icrm));
    yakl::atomicAdd(qrainsum1(k,j,i,icrm) , qrloud(k,j,i,icrm));
    yakl::atomicAdd(qicesum1(k,j,i,icrm) , qiloud(k,j,i,icrm));
    yakl::atomicAdd(prainsum1(k,j,i,icrm) , qrloud(k,j,i,icrm));
    yakl::atomicAdd(qsnowsum1(k,j,i,icrm) , qirloud(k,j,i,icrm));
    

/*
    yakl::atomicAdd(qcloudsum1(k,j,i,icrm) , qcloud(k,j,i,icrm));
    yakl::atomicAdd(qrainsum1(k,j,i,icrm) , qrloud(k,j,i,icrm));
    yakl::atomicAdd(qicesum1(k,j,i,icrm) , qiloud(k,j,i,icrm));
    yakl::atomicAdd(liran_test4d(k,j,i,icrm) , liran_test4d(k,j,i,icrm));
    yakl::atomicAdd(ecppwwsum1(k,j,i,icrm) , ecppwwsum1(k,j,i,icrm));
    yakl::atomicAdd(ecppwwsqsum1(k,j,i,icrm) , ecppwwsqsum1(k,j,i,icrm));
    yakl::atomicAdd(rhsum1(k,j,i,icrm) , rhsum1(k,j,i,icrm));
    yakl::atomicAdd(cf3dsum1(k,j,i,icrm) , cf3dsum1(k,j,i,icrm));
    yakl::atomicAdd(tkesgssum1(k,j,i,icrm) , tkesgssum1(k,j,i,icrm));
    yakl::atomicAdd(qvssum1(k,j,i,icrm) , qvssum1(k,j,i,icrm));
    yakl::atomicAdd(prainsum1(k,j,i,icrm) , prainsum1(k,j,i,icrm));
    yakl::atomicAdd(precallsum1(k,j,i,icrm) , precallsum1(k,j,i,icrm));
    */


});

/*
      qcloud_bfsum1(:,:,:,icrm) = qcloud_bfsum1(:,:,:,icrm) + qcloud_bf(:,:,:,icrm)
      qsnowsum1    (:,:,:,icrm) = qsnowsum1    (:,:,:,icrm) + qsnow(:,:,:,icrm)
      qgraupsum1   (:,:,:,icrm) = qgraupsum1   (:,:,:,icrm) + qgraup(:,:,:,icrm)
      qlsinksum1   (:,:,:,icrm) = qlsinksum1   (:,:,:,icrm) + qlsink(:,:,:,icrm)*qcloud(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg
      precrsum1    (:,:,:,icrm) = precrsum1    (:,:,:,icrm) + precr(:,:,:,icrm)
      precsolidsum1(:,:,:,icrm) = precsolidsum1(:,:,:,icrm) + precsolid(:,:,:,icrm)
      altsum1      (:,:,:,icrm) = altsum1      (:,:,:,icrm) + alt(:,:,:,icrm
      qlsink_bfsum1(:,:,:,icrm) = qlsink_bfsum1(:,:,:,icrm) + qlsink_bf(:,:,:,icrm)*qcloud_bf(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg

*/
printf("%s %.2f\n", "Liran check liran_test4d:", liran_test4d(10,10,10,10));
printf("%s %d\n", "Liran check itavg1:", itavg1);
printf("%s %d\n", "Liran check ntavg1:", ntavg1);
parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int icrm) {
  if (ntavg1 != 0 && crm_cnt(icrm) % ntavg1 == 0) {
      // itavg1 is divisible by ntavg1
    printf("%s %d %.2f \n", "Liran check liran_test4davg 0:", crm_cnt(icrm),liran_test4d(k,j,i,icrm));
    liran_test4d(k,j,i,icrm) = liran_test4d(k,j,i,icrm)/ntavg1;
    qcloudsum1(k,j,i,icrm)   = qcloudsum1(k,j,i,icrm)  /ntavg1;
    qrainsum1(k,j,i,icrm)    = qrainsum1(k,j,i,icrm)   /ntavg1;
    qicesum1(k,j,i,icrm)     = qicesum1(k,j,i,icrm)    /ntavg1;
    qsnowsum1(k,j,i,icrm)    = qsnowsum1(k,j,i,icrm)   /ntavg1;

    printf("%s %d %.2f \n", "Liran check liran_test4davg 1:", crm_cnt(icrm),liran_test4d(k,j,i,icrm));
  } else {
      // itavg1 is not divisible by ntavg1
    printf("itavg1 is not divisible by ntavg1\n");
  }
});
//printf("%s %d\n", "Liran check itavg1 div ntavg1:", itavg1 % ntavg1);
// Check if we have reached the end of the level 1 time averaging period.


// Increment the running sums for the level two variables that are not
// already incremented. Consolidate from 3-D to 1-D columns.

//------------------------------------------------------------------------------------------------
/*
Start of subroutine categorization_stats(
Transport classification is based on total condensate (cloudmixr_total), and
cloudy (liquid) and clear (non-liquid) classification is based on liquid water,
because wet deposition, aqueous chemistry, and droplet activaton, all are for liquid clouds.
Minghuai Wang, 2010-04
*/

parallel_for( "update sums",SimpleBounds<4>(nz, ny, nx, nens),
  YAKL_LAMBDA (int k, int j, int i, int icrm) {
    cloudmixr(k,j,i,icrm) = qcloudsum1(k,j,i,icrm);
    cloudmixr_total(k,j,i,icrm) = qcloudsum1(k,j,i,icrm) + qicesum1(k,j,i,icrm);
    // total hydrometer (rain, snow, and graupel)
    precmixr_total(k,j,i,icrm) = qrainsum1(k,j,i,icrm)+qsnowsum1(k,j,i,icrm); //+qsnow+qgraup
});





printf("Liran check start ECPP ecpp_crm_stat 02\n");
}








