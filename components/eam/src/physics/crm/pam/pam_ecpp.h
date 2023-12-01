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
  auto nzi         = coupler.get_option<int>("crm_nzi");  // Note that nz   = crm_nz
  auto nx           = coupler.get_option<int>("crm_nx");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");

  auto NCLASS_CL     = coupler.get_option<int>("ecpp_NCLASS_CL");
  auto ndraft_max     = coupler.get_option<int>("ecpp_ndraft_max");
  auto NCLASS_PR     = coupler.get_option<int>("ecpp_NCLASS_PR");
  printf("Liran check start ECPP init 2\n");
  auto gcm_dt       = coupler.get_option<real>("gcm_dt");
  auto crm_dt       = coupler.get_option<real>("crm_dt");
  
  int kbase, ktop;
  int m;

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

  real cloudthresh_trans = 1e-5;
  real precthresh_trans = 1e-4;

  int areaavgtype = 1;
  bool allcomb = false;

  int nupdraft = 0;
  int ndndraft = 0;
  ndraft_max = 0;

  int nupdraft_max = 0;  // Note: This variable is initialized but not used in the given code
  int ndndraft_max = 0;  // Note: This variable is initialized but not used in the given code

  // variables should inside ecppvars.h
  int DN1 = 0; // !First index of downward classes
  int NCLASS_TR = 0; // !Num. of transport classes
  int ncc_in       = 2; // Nnumber of clear/cloudy sub-calsses
  int nprcp_in     = 2; // Number of non-precipitating/precipitating sub-classes.
  NCLASS_CL = ncc_in; // Number of cloud classes
  NCLASS_PR = nprcp_in; // Number of precipitaion classes
   
 
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
printf("\nValue of nupdraft: %d: ", nupdraft);
printf("\nValue of ndndraft: %d: ", ndndraft);
printf("\nValue of nzstag: %d: ", nzstag);
printf("\nValue of NCLASS_CL: %d: ", NCLASS_CL);
printf("\nValue of NCLASS_PR: %d: ", NCLASS_PR);
printf("\nValue of nens: %d: ", nens);

coupler.set_option<int>("ecpp_NCLASS_CL",NCLASS_CL);
coupler.set_option<int>("ecpp_NCLASS_PR",NCLASS_PR);
coupler.set_option<int>("ecpp_ndraft_max",ndraft_max);

dm_device.register_and_allocate<int>("crm_cnt"     , "number of crm timestep count",  {nens},{"nens"});
dm_device.register_and_allocate<int>("ndown"       , "number of down count",  {nens},{"nens"});
dm_device.register_and_allocate<int>("nup"         , "number of nup count",  {nens},{"nens"});
dm_device.register_and_allocate<int>("kup_top"     , "maximum kup"  ,  {nens},{"nens"});
dm_device.register_and_allocate<int>("kdown_top"   , "maximum kdown",  {nens},{"nens"});
dm_device.register_and_allocate<real>("wdown_bar", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("wup_bar", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("wdown_stddev", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("wup_stddev", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("wup_rms", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("wdown_rms", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("updraftbase", "<description>",  {nens}, {"nens"});
dm_device.register_and_allocate<real>("updrafttop", "<description>" ,  {nens}, {"nens"});
dm_device.register_and_allocate<real>("dndrafttop", "<description>" ,  {nens}, {"nens"});
dm_device.register_and_allocate<real>("dndraftbase", "<description>",  {nens}, {"nens"});


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
dm_device.register_and_allocate<int>("nup_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_bar_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<int>("ndown_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_bar_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_stddev_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_stddev_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_rms_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_rms_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_rms_ksmo", "<description>", {nzi,nens}, {"zp1","nens"});
dm_device.register_and_allocate<real>("wdown_rms_ksmo", "<description>", {nzi,nens}, {"zp1","nens"});
printf("Liran check ECPP allocation done the first part\n");

auto crm_cnt         = dm_device.get<int,1>("crm_cnt");
auto ndown           = dm_device.get<int,1>("ndown");
auto nup             = dm_device.get<int,1>("nup");
auto kup_top         = dm_device.get<int,1>("kup_top");
auto kdown_top       = dm_device.get<int,1>("kdown_top");
auto wdown_bar       = dm_device.get<real,1>("wdown_bar");
auto wdown_stddev    = dm_device.get<real,1>("wdown_stddev");
auto wup_stddev      = dm_device.get<real,1>("wup_stddev");
auto wup_rms         = dm_device.get<real,1>("wup_rms");
auto wdown_rms       = dm_device.get<real,1>("wdown_rms");
auto wup_bar         = dm_device.get<real,1>("wup_bar");
auto updraftbase     = dm_device.get<real,1>("updraftbase");
auto updrafttop      = dm_device.get<real,1>("updrafttop");
auto dndrafttop      = dm_device.get<real,1>("dndrafttop");
auto dndraftbase     = dm_device.get<real,1>("dndraftbase");

printf("\nLiran check start updraftbase calculation\n");
for (int icrm = 0; icrm < nens; ++icrm) {
    crm_cnt (icrm) = 0;
    ndown (icrm) = 0;
    nup (icrm) = 0;
    kup_top(icrm) = 0;
    kdown_top(icrm) = 0;
    wup_bar(icrm) = 0.0;
    wdown_bar(icrm) = 0.0;
    wup_stddev(icrm) = 0.0;
    wdown_stddev(icrm) = 0.0;
    wup_rms(icrm) = 0.0;
    wdown_rms(icrm) = 0.0;
    updraftbase(icrm) = 0;
    updrafttop(icrm) = nz - 1;
    dndrafttop(icrm) = nz - 1;
    dndraftbase(icrm) = 0;
}
 printf("Liran check start ECPP allocation here\n");

auto qlsink_bf            = dm_device.get<real,4>("qlsink_bf");
auto prain                = dm_device.get<real,4>("prain");
auto qcloud_bf            = dm_device.get<real,4>("qcloud_bf");
auto qcloudsum1           = dm_device.get<real,4>("qcloudsum1");
auto qcloud_bfsum1        = dm_device.get<real,4>("qcloud_bfsum1");
auto qrainsum1            = dm_device.get<real,4>("qrainsum1");
auto qicesum1             = dm_device.get<real,4>("qicesum1");
auto qsnowsum1            = dm_device.get<real,4>("qsnowsum1");
auto qgraupsum1           = dm_device.get<real,4>("qgraupsum1");
auto qlsinksum1           = dm_device.get<real,4>("qlsinksum1");
auto precrsum1            = dm_device.get<real,4>("precrsum1");
auto precsolidsum1        = dm_device.get<real,4>("precsolidsum1");
auto precallsum1          = dm_device.get<real,4>("precallsum1");
auto altsum1              = dm_device.get<real,4>("altsum1");
auto rhsum1               = dm_device.get<real,4>("rhsum1");
auto cf3dsum1             = dm_device.get<real,4>("cf3dsum1");
auto ecppwwsum1           = dm_device.get<real,4>("ecppwwsum1");
auto ecppwwsqsum1         = dm_device.get<real,4>("ecppwwsqsum1");
auto tkesgssum1           = dm_device.get<real,4>("tkesgssum1");
auto qlsink_bfsum1        = dm_device.get<real,4>("qlsink_bfsum1");
auto prainsum1            = dm_device.get<real,4>("prainsum1");
auto qvssum1              = dm_device.get<real,4>("qvssum1");
auto liran_test4d         = dm_device.get<real,4>("liran_test4d");
auto liran_test4davg      = dm_device.get<real,4>("liran_test4davg");
auto xkhvsum              = dm_device.get<real,2>("xkhvsum");
auto nup_k                = dm_device.get<int,2>("nup_k");
auto ndown_k              = dm_device.get<int,2>("ndown_k");
auto wup_bar_k            = dm_device.get<real,2>("wup_bar_k");
auto wdown_bar_k          = dm_device.get<real,2>("wdown_bar_k");
auto wup_stddev_k         = dm_device.get<real,2>("wup_stddev_k");
auto wdown_stddev_k       = dm_device.get<real,2>("wdown_stddev_k");
auto wup_rms_k            = dm_device.get<real,2>("wup_rms_k");
auto wdown_rms_k          = dm_device.get<real,2>("wdown_rms_k");
auto wup_rms_ksmo         = dm_device.get<real,2>("wup_rms_ksmo");
auto wdown_rms_ksmo       = dm_device.get<real,2>("wdown_rms_ksmo");
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
  nup_k(iz,iens)                   = 0;
  wup_bar_k(iz,iens)               = 0;
  ndown_k(iz,iens)                 = 0;
  wdown_bar_k(iz,iens)             = 0;
  wdown_stddev_k(iz,iens)          = 0;
  wup_stddev_k(iz,iens)            = 0;
  wup_rms_k(iz,iens)               = 0;
  wdown_rms_k(iz,iens)             = 0;
  wup_rms_ksmo(iz,iens)            = 0;
  wdown_rms_ksmo(iz,iens)          = 0;
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
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nz              = coupler.get_option<int>("crm_nz");    // Note that nz   = crm_nz
  auto nzi            = coupler.get_option<int>("crm_nzi");  // Note that nzi  = crm_nz+1
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  auto ntavg1          = coupler.get_option<int>("ecpp_ntavg1");
  auto ntavg2          = coupler.get_option<int>("ecpp_ntavg2");
  auto mode_updnthresh = coupler.get_option<int>("mode_updnthresh");
  auto plumetype       = coupler.get_option<int>("plumetype");
  auto NCLASS_CL       = coupler.get_option<int>("ecpp_NCLASS_CL");
  auto ndraft_max      = coupler.get_option<int>("ecpp_ndraft_max");
  auto NCLASS_PR       = coupler.get_option<int>("ecpp_NCLASS_PR");
  auto gcm_dt          = coupler.get_option<real>("gcm_dt");
  auto crm_dt          = coupler.get_option<real>("crm_dt");
  
  printf("%s %.2f\n", "Liran check gcm_dt:",gcm_dt);
  printf("%s %.2f\n", "Liran check crm_dt:", crm_dt);
  printf("Liran check start ECPP stage2:\n");
  printf("\nValue of ndraft_max 2: %d: ", ndraft_max);
  printf("\nValue of NCLASS_CL 2: %d: ", NCLASS_CL);
  printf("\nValue of NCLASS_PR 2: %d: ", NCLASS_PR);
  //------------------------------------------------------------------------------------------------
  // get variables allocated in ecpp_crm_init
  auto crm_cnt                  = dm_device.get<int,1>("crm_cnt");
  auto ndown                    = dm_device.get<int,1>("ndown");
  auto nup                      = dm_device.get<int,1>("nup");
  auto kup_top                  = dm_device.get<int,1>("kup_top");
  auto kdown_top                = dm_device.get<int,1>("kdown_top");
  auto nup_k                    = dm_device.get<int, 2>("nup_k");
  auto ndown_k                  = dm_device.get<int, 2>("ndown_k");
  auto wdown_bar                = dm_device.get<real,1>("wdown_bar");
  auto wdown_stddev             = dm_device.get<real,1>("wdown_stddev");
  auto wup_stddev               = dm_device.get<real,1>("wup_stddev");
  auto wup_rms                  = dm_device.get<real,1>("wup_rms");
  auto wdown_rms                = dm_device.get<real,1>("wdown_rms");
  auto wup_bar                  = dm_device.get<real,1>("wup_bar");
  auto qcloudsum1               = dm_device.get<real, 4>("qcloudsum1");
  auto qcloud_bfsum1            = dm_device.get<real, 4>("qcloud_bfsum1");
  auto qrainsum1                = dm_device.get<real, 4>("qrainsum1");
  auto qicesum1                 = dm_device.get<real, 4>("qicesum1");
  auto qsnowsum1                = dm_device.get<real, 4>("qsnowsum1");
  auto qgraupsum1               = dm_device.get<real, 4>("qgraupsum1");
  auto qlsinksum1               = dm_device.get<real, 4>("qlsinksum1");
  auto precrsum1                = dm_device.get<real, 4>("precrsum1");
  auto precsolidsum1            = dm_device.get<real, 4>("precsolidsum1");
  auto precallsum1              = dm_device.get<real, 4>("precallsum1");
  auto altsum1                  = dm_device.get<real, 4>("altsum1");
  auto rhsum1                   = dm_device.get<real, 4>("rhsum1");
  auto cf3dsum1                 = dm_device.get<real, 4>("cf3dsum1");
  auto ecppwwsum1               = dm_device.get<real, 4>("ecppwwsum1");
  auto ecppwwsqsum1             = dm_device.get<real, 4>("ecppwwsqsum1");
  auto tkesgssum1               = dm_device.get<real, 4>("tkesgssum1");
  auto qlsink_bfsum1            = dm_device.get<real, 4>("qlsink_bfsum1");
  auto prainsum1                = dm_device.get<real, 4>("prainsum1");
  auto qvssum1                  = dm_device.get<real, 4>("qvssum1");
  auto wup_bar_k                = dm_device.get<real, 2>("wup_bar_k");
  auto wdown_bar_k              = dm_device.get<real, 2>("wdown_bar_k");
  auto wdown_stddev_k           = dm_device.get<real, 2>("wdown_stddev_k");
  auto wup_stddev_k             = dm_device.get<real, 2>("wup_stddev_k");
  auto wup_rms_k                = dm_device.get<real, 2>("wup_rms_k");
  auto wdown_rms_k              = dm_device.get<real, 2>("wdown_rms_k");
  auto wup_rms_ksmo             = dm_device.get<real, 2>("wup_rms_ksmo");
  auto wdown_rms_ksmo           = dm_device.get<real, 2>("wdown_rms_ksmo");
  auto liran_test4d             = dm_device.get<real, 4>("liran_test4d");
  auto liran_test4davg          = dm_device.get<real, 4>("liran_test4davg");

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
  auto ref_pres                 = dm_device.get<real,2>("ref_pres");
  auto crm_wvel                 = dm_device.get<real,4>("wvel");
  auto updraftbase              = dm_device.get<real,1>("updraftbase");
  auto updrafttop               = dm_device.get<real,1>("updrafttop");
  auto dndrafttop               = dm_device.get<real,1>("dndrafttop");
  auto dndraftbase              = dm_device.get<real,1>("dndraftbase");

  //------------------------------------------------------------------------------------------------
  // Define variables used by subroutine categorization_stats
  real7d mask_bnd("mask_bnd",nzi,ny,nx,nens,NCLASS_CL,ndraft_max,NCLASS_PR);
  real7d mask_cen("mask_cen",nzi,ny,nx,nens,NCLASS_CL,ndraft_max,NCLASS_PR);
  real4d cloudmixr("cloudmixr",nz,ny,nx,nens);
  real4d cloudmixr_total("cloudmixr_total",nz,ny,nx,nens);
  real4d precmixr_total("precmixr_total",nz,ny,nx,nens);
  real4d qvs("qvs",nz,ny,nx,nens);
  real4d alt("alt",nz,ny,nx,nens);
  real3d cloudtop("cloudtop"   ,ny,nx,nens);
  real3d cloudtop_upaa("cloudtop_upaa",ny,nx,nens);
  real3d cloudtop_downaa("cloudtop_downaa",ny,nx,nens);
  real3d cloudtop_upbb("cloudtop_upbb",ny,nx,nens);
  real3d cloudtop_downbb("cloudtop_downbb",ny,nx,nens);
  real3d wup_thresh_k("wup_thresh_k",nz,2,nens);
  real3d wdown_thresh_k("wdown_thresh_k",nz,2,nens);
  real2d rhoair("rhoair",nzi,nens); //layer-averaged air density
  real2d tmpveca("rhoair",nz,nens);
  real2d tmpvecb("rhoair",nz,nens);
  real2d wup_thresh("wup_thresh",nz,nens);
  real2d wdown_thresh("wdown_thresh",nz,nens);

  real upthresh = 1.0;
  real downthresh = 1.0;
  real upthresh2 = 0.5;
  real downthresh2 = 0.5;
  real cloudthresh = 1e-6;
  real prcpthresh = 1e-6;
  // mhwang
  // high thresholds are used to classify transport classes (following Xu et al., 2002, Q.J.R.M.S.
  real cloudthresh_trans = 1e-5; //Cloud mixing ratio beyond which cell is "cloudy" to classify transport classes (kg/kg)   +++mhwang
  // the maxium of cloudthres_trans and 0.01*qvs is used to classify transport class
  real precthresh_trans  = 1e-4; //Preciptation mixing ratio beyond which cell is raining to classify transport classes (kg/kg)  !+++mwhang
  int runcount     = 0;
  int ijdel_upaa   = 0;
  int ijdel_downaa = 0;
  int ijdel_upbb   = 0;
  int ijdel_downbb = 0;

  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    crm_cnt(iens) = crm_cnt(iens) + 1;
  });
  runcount = crm_cnt(0);
  printf("%s %d\n", "Liran check runcount:", runcount);
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

  // Set level-1 and level-2 averaging periods for ECPP
  
  // Calculate number of steps assuming dt evenly divides ntavg[12]_ss
/* 
!------------------------------------------------------------------------
! Main code section...
!------------------------------------------------------------------------
*/

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

    real EVS = esatw_crm(crm_temp(k,j,i,icrm)); //   ! saturation water vapor pressure (PA)
    qvs(k,j,i,icrm) = .622*EVS/(ref_pres(icrm,k)*100.-EVS); //  ! pres(icrm,kk) with unit of hPa
    alt(k,j,i,icrm) =  287.0*crm_temp(k,j,i,icrm)/(100.*ref_pres(icrm,k));
    yakl::atomicAdd(altsum1(k,j,i,icrm) , alt(k,j,i,icrm));
    yakl::atomicAdd(liran_test4d(k,j,i,icrm) , liran_test4d(k,j,i,icrm));
    yakl::atomicAdd(qcloudsum1(k,j,i,icrm) , qcloud(k,j,i,icrm));
    yakl::atomicAdd(qrainsum1(k,j,i,icrm) , qrloud(k,j,i,icrm));
    yakl::atomicAdd(qicesum1(k,j,i,icrm) , qiloud(k,j,i,icrm));
    yakl::atomicAdd(prainsum1(k,j,i,icrm) , qrloud(k,j,i,icrm));
    yakl::atomicAdd(qsnowsum1(k,j,i,icrm) , qsnowsum1(k,j,i,icrm)); // This is ZERO!! for now
    yakl::atomicAdd(ecppwwsum1(k,j,i,icrm) , crm_wvel(k,j,i,icrm));
    

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
if (ntavg1 != 0 && runcount % ntavg1 == 0) {
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (ntavg1 != 0 && crm_cnt(icrm) % ntavg1 == 0) {
        // itavg1 is divisible by ntavg1
      printf("%s %d %.2f \n", "Liran check liran_test4davg 0:", crm_cnt(icrm),liran_test4d(k,j,i,icrm));
      liran_test4d(k,j,i,icrm) = liran_test4d(k,j,i,icrm)/ntavg1;
      qcloudsum1(k,j,i,icrm)   = qcloudsum1(k,j,i,icrm)  /ntavg1;
      qrainsum1(k,j,i,icrm)    = qrainsum1(k,j,i,icrm)   /ntavg1;
      qicesum1(k,j,i,icrm)     = qicesum1(k,j,i,icrm)    /ntavg1;
      qsnowsum1(k,j,i,icrm)    = qsnowsum1(k,j,i,icrm)   /ntavg1;
      altsum1(k,j,i,icrm)      = altsum1(k,j,i,icrm)     /ntavg1;
      ecppwwsum1(k,j,i,icrm)   = ecppwwsum1(k,j,i,icrm)  /ntavg1;
      printf("%s %d %.2f \n", "Liran check liran_test4davg 1:", crm_cnt(icrm),liran_test4d(k,j,i,icrm));
    } else {
        // itavg1 is not divisible by ntavg1
      printf("%s %d\n", "\nitavg1 is not divisible by ntavg1", crm_cnt(icrm));
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
  printf("Liran check categorization_stats\n");
  parallel_for( "update sums",SimpleBounds<4>(nz, ny, nx, nens),
    YAKL_LAMBDA (int k, int j, int i, int icrm) {
      cloudmixr(k,j,i,icrm) = qcloudsum1(k,j,i,icrm);
      cloudmixr_total(k,j,i,icrm) = qcloudsum1(k,j,i,icrm) + qicesum1(k,j,i,icrm);
      // total hydrometer (rain, snow, and graupel)
      precmixr_total(k,j,i,icrm) = qrainsum1(k,j,i,icrm)+qsnowsum1(k,j,i,icrm); //+qsnow+qgraup
  });
  int nxy = nx*ny;
  parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int icrm) {
    for (int k_gcm=0; k_gcm<nzi; k_gcm++) {
      //int l = plev-(k+1);
      int k_crm= (gcm_nlev+1)-1-k_gcm;
      /*
      ! Get cloud top height
      ! Cloud top height is used to determine whether there is updraft/downdraft. No updraft and
      ! downdraft is allowed above the condensate level (both liquid and ice).
      */
      cloudtop(j,i,icrm) = 1; // !Default to bottom level if no cloud in column.
      if (cloudmixr_total(k_crm,j,i,icrm) >= cloudthresh_trans) {
        cloudtop(j,i,icrm) = k_crm; 
      }
    }
    for (int k_crm=0; k_crm<nzi; k_crm++) {
      int km0 = std::min(nz,k_crm);
      int km1 = std::max(1,k_crm-1);
      rhoair(k_crm,icrm) = rhoair(k_crm,icrm)+0.5*(1.0/altsum1(km1,j,i,icrm) + 1.0/altsum1(km0,j,i,icrm))/nxy;
    }
  });

  printf("Liran check categorization_stats 2\n");
   //------------------------------------------------------------------------------------------------
    //parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
      //crm_cnt(iens) = crm_cnt(iens) + 1;
    //});

  /*
    !------------------------------------------------------------------------
    subroutine determine_transport_thresh( &
      nx, ny, nz, &
      mode_updnthresh, upthresh, downthresh, &
      upthresh2, downthresh2, cloudthresh, &
      !     ctime, &
      ww, rhoair, &
      wdown_thresh_k, wup_thresh_k         &
      , cloudtop                           &
      , wup_rms_k, wup_bar_k, wup_stddev_k  &
      , wdown_rms_k, wdown_bar_k, wdown_stddev_k  &
      , kup_top, kdown_top)
      !
      ! Deterines the velocity thresholds used to indicate whether a cell's
      ! motion is up, down, or quiescent. This is down for two threshold values
      ! in each direction by level. A dozen options are available on how this
      ! is done as documented below and at the top of postproc_wrfout.
      !
      ! William.Gustafosn@pnl.gov; 11-Sep-2008
      ! Modified: William.Gustafosn@pnl.gov; 14-Apr-2009
      !------------------------------------------------------------------------
  */

  /*
      ! Calc cloudtop_upaa(i,j) = max( cloudtop(i-del:i+del,j-del:j+del) )
      ! and similar for cloudtop_upbb, cloudtop_downaa/bb
      ! (assume periodic BC here)
  */

  printf("Liran check determine_transport_thresh 1\n");

  // if ((mode_updnthresh == 12) .or. (mode_updnthresh == 13)) then This is ignored
  int ijdel = std::max({ijdel_upaa, ijdel_upbb, ijdel_downaa, ijdel_downbb});

  printf("\nValue of mode_updnthresh: %d: ", mode_updnthresh);
  // Value of mode_updnthresh: 16
  printf("\nValue of ijdel: %d: ", ijdel);
  // Value of ijdel: 0

  //if (ijdel > 0) then remove line 972 to 1008

  /*
      ! new coding here and below
      !   cloudtop_up/downaa - only grid cells with k<=cloudtop_up/downaa
      !                        are used for calc of wup_rms and wdn_rms
      !   cloudtop_up/downbb - only grid cells with k<=cloudtop_up/downbb
      !                        can be classified as up/downdraft
  */

  // if ((mode_updnthresh == 12) .or. (mode_updnthresh == 13)) then remove line 1015 to 1021
  /*
      ! mode_updnthresh /= 12,13 corresponds to pre 11-jan-2008 versions of preprocessor
      !   where only grid cells with k <= cloudtop(i,j) are used for calc of wup/dn_rms,
      !   but any grid cells can be up/dn [even those with k >> cloudtop(i,j)]
  */

  parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int icrm) {
    cloudtop_upaa(j,i,icrm)   = cloudtop(j,i,icrm);
    cloudtop_downaa(j,i,icrm) = cloudtop(j,i,icrm);
    cloudtop_upbb(j,i,icrm)   = nz;
    cloudtop_downbb(j,i,icrm) = nz;
  });

  /*
      ! Get standard deviation of up and down vertical velocity below the
      ! cloud tops. For now, each cell is treated equally. We may want to
      ! consider weighting each cell by its volume or mass.
      !
      ! Get the mean values first for wup and wdown
  */

  parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int icrm) {
    
    for (int k_crm=0; k_crm<cloudtop_upaa(j,i,icrm); k_crm++) {
  /*
            !It is dimmensionally ok since w is dimmed nz+1
            !We intentially ignore when w==0 as to not bias one direction
            !over the other for the count. This differs from the Ferret code which
            !assigns w=0 to up values.
  */
      if (ecppwwsum1(k_crm,j,i,icrm) > 0.0) {
        nup(icrm) = nup(icrm) + 1;
        wup_bar(icrm) = wup_bar(icrm) + ecppwwsum1(k_crm,j,i,icrm);
        nup_k(k_crm,icrm) = nup_k(k_crm,icrm) + 1;
        wup_bar_k(k_crm,icrm) = wup_bar_k(k_crm,icrm) + ecppwwsum1(k_crm,j,i,icrm);
        kup_top(icrm)   = std::max(kup_top(icrm), k_crm);
      }
    }

    for (int k_crm=0; k_crm<cloudtop_downaa(j,i,icrm); k_crm++) {
      if (ecppwwsum1(k_crm,j,i,icrm) < 0.0) {
        ndown(icrm) = ndown(icrm) + 1;
        wdown_bar(icrm) = wdown_bar(icrm) + ecppwwsum1(k_crm,j,i,icrm);
        ndown_k(k_crm,icrm) = ndown_k(k_crm,icrm) + 1;
        wdown_bar_k(k_crm,icrm) = wdown_bar_k(k_crm,icrm) + ecppwwsum1(k_crm,j,i,icrm);
        kdown_top(icrm)   = std::max(kdown_top(icrm), k_crm);
      }
    }

  });

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    if (nup(icrm) > 0.0) {
      wup_bar(icrm)   = wup_bar(icrm) / nup(icrm);
    }
    if (ndown(icrm) > 0.0) {
      wdown_bar(icrm)   = wdown_bar(icrm) / ndown(icrm);
    }
  });

  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k_crm, int icrm) {
    if (nup_k(k_crm,icrm) > 0.0) {
      wup_bar_k(k_crm,icrm)   = wup_bar_k(k_crm,icrm) / nup_k(k_crm,icrm);
    }
    if (ndown_k(k_crm,icrm) > 0.0) {
      wdown_bar_k(k_crm,icrm)   = wdown_bar_k(k_crm,icrm) / ndown_k(k_crm,icrm);
    }
  });

  // !Now, we can get the std. dev. of wup and wdown.
  parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int icrm) {
    
    for (int k_crm=0; k_crm<cloudtop_upaa(j,i,icrm); k_crm++) {
  /*
            !We intentionally ignore when w==0 as to not bias one direction
            !over the other.
  */
      if (ecppwwsum1(k_crm,j,i,icrm) > 0.0) {
        wup_stddev(icrm) = wup_stddev(icrm) + (ecppwwsum1(k_crm,j,i,icrm)-wup_bar(icrm))*(ecppwwsum1(k_crm,j,i,icrm)-wup_bar(icrm));
        wup_stddev_k(k_crm,icrm) = wup_stddev_k(k_crm,icrm) + (ecppwwsum1(k_crm,j,i,icrm)-wup_bar_k(k_crm,icrm))*(ecppwwsum1(k_crm,j,i,icrm)-wup_bar_k(k_crm,icrm));
      }
    }

    for (int k_crm=0; k_crm<cloudtop_downaa(j,i,icrm); k_crm++) {
      if (ecppwwsum1(k_crm,j,i,icrm) < 0.0) {
        wdown_stddev(icrm) = wdown_stddev(icrm) + (ecppwwsum1(k_crm,j,i,icrm)-wdown_bar(icrm))*(ecppwwsum1(k_crm,j,i,icrm)-wdown_bar(icrm));
        wdown_stddev_k(k_crm,icrm) = wdown_stddev_k(k_crm,icrm) + (ecppwwsum1(k_crm,j,i,icrm)-wdown_bar_k(k_crm,icrm))*(ecppwwsum1(k_crm,j,i,icrm)-wdown_bar_k(k_crm,icrm));
      }
    }
  });

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    if (nup(icrm) > 0.0) {
      wup_stddev(icrm)   = wup_stddev(icrm) / nup(icrm);
    }
    if (ndown(icrm) > 0.0) {
      wdown_stddev(icrm)   = wdown_stddev(icrm) / ndown(icrm);
    }
  });

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    wup_rms(icrm)   =  std::sqrt(wup_bar(icrm)*wup_bar(icrm)+wup_stddev(icrm)*wup_stddev(icrm));
    wdown_rms(icrm)   = std::sqrt(wdown_bar(icrm)*wdown_bar(icrm)+wdown_stddev(icrm)*wdown_stddev(icrm));
  });

  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k_crm, int icrm) {
    if (nup_k(k_crm,icrm) > 0.0) {
      wup_stddev_k(k_crm,icrm)   = wup_stddev_k(k_crm,icrm) / nup_k(k_crm,icrm);
    }
    if (ndown_k(k_crm,icrm) > 0.0) {
      wdown_stddev_k(k_crm,icrm)   = wdown_stddev_k(k_crm,icrm) / ndown_k(k_crm,icrm);
    }
    wup_rms_k(k_crm,icrm) = std::sqrt( wup_bar_k(k_crm,icrm)*wup_bar_k(k_crm,icrm) + wup_stddev_k(k_crm,icrm)*wup_stddev_k(k_crm,icrm) );
    wdown_rms_k(k_crm,icrm) = std::sqrt( wdown_bar_k(k_crm,icrm)*wdown_bar_k(k_crm,icrm) + wdown_stddev_k(k_crm,icrm)*wdown_stddev_k(k_crm,icrm) );
  });

  // ! calculated smoothed (3-point) wup/down_rms
  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k_crm, int icrm) {
    tmpveca(k_crm,icrm) = wup_rms_k(k_crm,icrm);
    tmpvecb(k_crm,icrm) = wdown_rms_k(k_crm,icrm);
  });


  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    for (int k_crm=1; k_crm<nz; k_crm++) {
      wup_rms_ksmo(k_crm,icrm) = 0.0;
      wdown_rms_ksmo(k_crm,icrm) = 0.0;
      real tmpsuma = 0.0;
      for (int k_crm2=k_crm-1; k_crm<=k_crm+1; k_crm++) {
        wup_rms_ksmo(k_crm,icrm) = wup_rms_ksmo(k_crm,icrm) + tmpveca(k_crm2,icrm);
        wdown_rms_ksmo(k_crm,icrm) = wdown_rms_ksmo(k_crm,icrm) + tmpvecb(k_crm2,icrm);
        tmpsuma = tmpsuma + 1.0;
      }
      tmpsuma = std::max(tmpsuma,1.0);
      wup_rms_ksmo(k_crm,icrm) = wup_rms_ksmo(k_crm,icrm)/tmpsuma;
      wdown_rms_ksmo(k_crm,icrm) = wdown_rms_ksmo(k_crm,icrm)/tmpsuma;
    }
    wup_rms_ksmo(0,icrm) = wup_rms_ksmo(1,icrm);
    wdown_rms_ksmo(0,icrm) = wdown_rms_ksmo(1,icrm);
    wup_rms_ksmo(nzi,icrm) = wup_rms_ksmo(nz,icrm);
    wdown_rms_ksmo(nzi,icrm) = wdown_rms_ksmo(nz,icrm);
  });
  printf("Liran check determine_transport_thresh 2\n");
  /*
        ! case 16 & 17 -- added on 10-dec-2009
        !    updraft   and k  > "updraft   center k",  use max( wup_rms_k, wup_rms )
        !    updraft   and k <= "updraft   center k",  use wup_rms_k
        !    downdraft and k  > "downdraft center k",  use max( wdown_rms_k, wdown_rms )
        !    downdraft and k <= "downdraft center k",  use wdown_rms_k
        ! The idea is to have a higher threshold in upper troposphere to
        ! filter out gravity waves motions
  */

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    real tmpsuma = 0.0;
    real tmpw = 0.0;
    real tmpw_minval = 0.10;
    real tmpsumb = 1.0e-30;
    for (int k_crm=0; k_crm<=nz; k_crm++) {
      tmpw = wup_rms_k(k_crm,icrm);
      tmpw = std::max(1.0e-4,tmpw);
      tmpw = tmpw * rhoair(k_crm,icrm);
      tmpsuma = tmpsuma + tmpw*k_crm; 
      tmpsumb = tmpsumb + tmpw;
    }
    int kup_center1 = std::round(tmpsuma / tmpsumb);
    for (int k_crm=0; k_crm<=nz; k_crm++) {
      tmpw = wup_rms_k(k_crm,icrm);
      if (k_crm > kup_center1){
       tmpw = std::max( tmpw, wup_rms(icrm) );
     }
      tmpw = std::max( tmpw, tmpw_minval );
      wup_thresh_k(k_crm,0,icrm) = tmpw*std::abs(upthresh);
      wup_thresh_k(k_crm,1,icrm) = tmpw*std::abs(upthresh2);
    }
  });

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    real tmpsuma = 0.0;
    real tmpw = 0.0;
    real tmpw_minval = 0.10;
    real tmpsumb = 1.0e-30;
    for (int k_crm=0; k_crm<=nz; k_crm++) {
      tmpw = wdown_rms_k(k_crm,icrm);
      tmpw = std::max(1.0e-4,tmpw);
      tmpw = tmpw * rhoair(k_crm,icrm);
      tmpsuma = tmpsuma + tmpw*k_crm; 
      tmpsumb = tmpsumb + tmpw;
    }
    int kup_center2 = std::round(tmpsuma / tmpsumb);
    for (int k_crm=0; k_crm<=nz; k_crm++) {
      tmpw = wdown_rms_k(k_crm,icrm);
      if (k_crm > kup_center2) { 
        tmpw = std::max( tmpw, wup_rms(icrm) );
      }
      tmpw = std::max( tmpw, tmpw_minval );
      wdown_thresh_k(k_crm,0,icrm) = tmpw*std::abs(upthresh);
      wdown_thresh_k(k_crm,1,icrm) = tmpw*std::abs(upthresh2);
    }

  });

  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int icrm) {
    wdown_thresh(k_crm,icrm) = wdown_thresh_k(k_crm,0,icrm);
    wup_thresh(k_crm,icrm) = wup_thresh_k(k_crm,0,icrm);
  });

}


printf("\nLiran check start ECPP ecpp_crm_stat 02\n");
}








