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
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto &dm_host     = coupler.get_data_manager_host_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");  // Note that nz   = crm_nz
  auto nzi         = coupler.get_option<int>("crm_nzi");  // Note that nz   = crm_nz
  auto nx           = coupler.get_option<int>("crm_nx");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");

  auto NCLASS_CL     = coupler.get_option<int>("ecpp_NCLASS_CL");
  auto ndraft_max    = coupler.get_option<int>("ecpp_ndraft_max");
  auto NCLASS_PR     = coupler.get_option<int>("ecpp_NCLASS_PR");
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
// Liran: For other mode_updnthresh options, the code need to change!!!!
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

  int CLR = 0;       // Clear sub-class
  int CLD = 1;       // Cloudy sub-class
  int PRN = 0;       // Not precipitating sub-class
  int PRY = 1;       // Is precipitating sub-class
  int DN1 = 0;       // !First index of downward classes
  int NCLASS_TR = 0; // !Num. of transport classes <value defined at line 122>
  int QUI       = 1; // Quiescent class
  int UP1       = 2; // First index for upward classes
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

//printf("Liran check start ECPP:\n");
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
dm_device.register_and_allocate<int>("crm_level2_cnt"     , "number of crm level 2 average count",  {nens},{"nens"});
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
dm_device.register_and_allocate<real>("wwqui_cen", "<description>", {nz,nens}, {"z","nens"});\
dm_device.register_and_allocate<real>("liran_test2d", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("cldtot2d", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wwqui_bnd", "<description>", {nzi,nens}, {"zi","nens"}); //nz+1
dm_device.register_and_allocate<real>("wwqui_cloudy_cen", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wwqui_cloudy_bnd", "<description>", {nzi,nens}, {"zi","nens"});//nz+1
dm_device.register_and_allocate<int>("nup_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_bar_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<int>("ndown_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_bar_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_stddev_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_stddev_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_rms_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wdown_rms_k", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("wup_rms_ksmo", "<description>", {nzi,nens}, {"zi","nens"});
dm_device.register_and_allocate<real>("wdown_rms_ksmo", "<description>", {nzi,nens}, {"zi","nens"});

dm_device.register_and_allocate<real>("ecpp_cat_wwqui_bar_cen", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_wwqui_cloudy_bar_cen", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_tbeg", "<description>",  {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_wwqui_bar_bnd", "<description>", {nzi,nens}, {"zi","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_wwqui_cloudy_bar_bnd", "<description>", {nzi,nens}, {"zi","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_area_cen_final", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_area_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_rh_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_qcloud_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_qice_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_precsolidcen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_area_bnd_final", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_area_bnd", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});
dm_device.register_and_allocate<real>("ecpp_cat_mass_bnd", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});

dm_device.register_and_allocate<real>("ecpp_sum_wwqui_bar_cen", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_wwqui_cloudy_bar_cen", "<description>", {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_tbeg", "<description>",  {nz,nens}, {"z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_wwqui_bar_bnd", "<description>", {nzi,nens}, {"zi","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_wwqui_cloudy_bar_bnd", "<description>", {nzi,nens}, {"zi","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_area_cen_final", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_area_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_rh_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_qcloud_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_qice_cen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_precsolidcen", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","z","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_area_bnd_final", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_area_bnd", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});
dm_device.register_and_allocate<real>("ecpp_sum_mass_bnd", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","zi","nens"});


auto crm_cnt         = dm_device.get<int,1>("crm_cnt");
auto crm_level2_cnt  = dm_device.get<int,1>("crm_level2_cnt");
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

for (int icrm = 0; icrm < nens; ++icrm) {
    crm_cnt (icrm) = 0;
    crm_level2_cnt (icrm) = 0;
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
 printf("\nValue of crm_cnt: %d: ", crm_cnt(0));
 printf("\nValue of crm_level2_cnt: %d: ", crm_level2_cnt(0));

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
auto wwqui_cen            = dm_device.get<real,2>("wwqui_cen");
auto wwqui_bnd            = dm_device.get<real,2>("wwqui_bnd");
auto liran_test2d         = dm_device.get<real,2>("liran_test2d");
auto cldtot2d             = dm_device.get<real,2>("cldtot2d");
auto wwqui_cloudy_cen     = dm_device.get<real,2>("wwqui_cloudy_cen");
auto wwqui_cloudy_bnd     = dm_device.get<real,2>("wwqui_cloudy_bnd");

auto ecpp_output_wwqui_cen         = dm_host.get<real,2>("ecpp_output_wwqui_cen");
auto ecpp_output_wwqui_cloudy_cen  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_cen");
auto ecpp_output_wwqui_bnd         = dm_host.get<real,2>("ecpp_output_wwqui_bnd");
auto ecpp_output_wwqui_cloudy_bnd  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_bnd");

auto ecpp_sum_wwqui_bar_cen         = dm_device.get<real,2>("ecpp_sum_wwqui_bar_cen");
auto ecpp_sum_wwqui_cloudy_bar_cen  = dm_device.get<real,2>("ecpp_sum_wwqui_cloudy_bar_cen");
auto ecpp_sum_wwqui_bar_bnd         = dm_device.get<real,2>("ecpp_sum_wwqui_bar_bnd");
auto ecpp_sum_wwqui_cloudy_bar_bnd  = dm_device.get<real,2>("ecpp_sum_wwqui_cloudy_bar_bnd");
auto ecpp_sum_tbeg                  = dm_device.get<real,2>("ecpp_sum_tbeg");
auto ecpp_sum_area_cen_final        = dm_device.get<real,5>("ecpp_sum_area_cen_final");
auto ecpp_sum_area_cen              = dm_device.get<real,5>("ecpp_sum_area_cen");
auto ecpp_sum_area_bnd_final        = dm_device.get<real,5>("ecpp_sum_area_bnd_final");
auto ecpp_sum_area_bnd              = dm_device.get<real,5>("ecpp_sum_area_bnd");
auto ecpp_sum_mass_bnd              = dm_device.get<real,5>("ecpp_sum_mass_bnd");
auto ecpp_sum_rh_cen                = dm_device.get<real,5>("ecpp_sum_rh_cen");
auto ecpp_sum_qcloud_cen            = dm_device.get<real,5>("ecpp_sum_qcloud_cen");
auto ecpp_sum_qice_cen              = dm_device.get<real,5>("ecpp_sum_qice_cen");
auto ecpp_sum_precsolidcen          = dm_device.get<real,5>("ecpp_sum_precsolidcen");

auto ecpp_cat_wwqui_bar_cen         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_cen");
auto ecpp_cat_wwqui_bar_bnd         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_bnd");
auto ecpp_cat_wwqui_cloudy_bar_cen  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_cen");
auto ecpp_cat_wwqui_cloudy_bar_bnd  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_bnd");
auto ecpp_cat_tbeg                  = dm_device.get<real,2>("ecpp_cat_tbeg");
auto ecpp_cat_area_cen_final        = dm_device.get<real,5>("ecpp_cat_area_cen_final");
auto ecpp_cat_area_cen              = dm_device.get<real,5>("ecpp_cat_area_cen");
auto ecpp_cat_area_bnd_final        = dm_device.get<real,5>("ecpp_cat_area_bnd_final");
auto ecpp_cat_area_bnd              = dm_device.get<real,5>("ecpp_cat_area_bnd");
auto ecpp_cat_mass_bnd              = dm_device.get<real,5>("ecpp_cat_mass_bnd");
auto ecpp_cat_rh_cen                = dm_device.get<real,5>("ecpp_cat_rh_cen");
auto ecpp_cat_qcloud_cen            = dm_device.get<real,5>("ecpp_cat_qcloud_cen");
auto ecpp_cat_qice_cen              = dm_device.get<real,5>("ecpp_cat_qice_cen");
auto ecpp_cat_precsolidcen          = dm_device.get<real,5>("ecpp_cat_precsolidcen");

//auto wwqui_cloudy_bnd = dm_device.get<real,2>("wwqui_cloudy_bnd");     // Note: Adjust dimensionality if needed
//auto wup_thresh           = dm_device.get<real,2>("wup_thresh");               // Note: Adjust dimensionality if needed
//auto wdown_thresh         = dm_device.get<real,2>("wdown_thresh");             // Note: Adjust dimensionality if needed

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
// Initialization of 2D variables
parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int iz, int iens) {
  xkhvsum(iz,iens)                        = 0;
  nup_k(iz,iens)                          = 0;
  wup_bar_k(iz,iens)                      = 0;
  ndown_k(iz,iens)                        = 0;
  wdown_bar_k(iz,iens)                    = 0;
  wdown_stddev_k(iz,iens)                 = 0;
  wup_stddev_k(iz,iens)                   = 0;
  wup_rms_k(iz,iens)                      = 0;
  wdown_rms_k(iz,iens)                    = 0;
  wwqui_cen(iz,iens)                      = 0;
  liran_test2d(iz,iens)                   = 1;
  cldtot2d(iz,iens)                       = 0.0;
  wwqui_cloudy_cen(iz,iens)               = 0;
  ecpp_output_wwqui_cen(iz,iens)          = 0;
  ecpp_output_wwqui_cloudy_cen(iz,iens)   = 0;
  ecpp_cat_wwqui_bar_cen(iz,iens)         = 0;
  ecpp_cat_wwqui_cloudy_bar_cen(iz,iens)  = 0;
  ecpp_cat_tbeg(iz,iens)                  = 0;
  ecpp_sum_wwqui_bar_cen(iz,iens)         = 0;
  ecpp_sum_wwqui_cloudy_bar_cen(iz,iens)  = 0;
  ecpp_sum_tbeg(iz,iens)                  = 0;
});

parallel_for(SimpleBounds<2>(nzi,nens), YAKL_LAMBDA (int iz, int iens) {
  wwqui_bnd(iz,iens)                      = 0; 
  wwqui_cloudy_bnd(iz,iens)               = 0;
  wup_rms_ksmo(iz,iens)                   = 0;
  wdown_rms_ksmo(iz,iens)                 = 0;
  ecpp_output_wwqui_bnd(iz,iens)          = 0;
  ecpp_output_wwqui_cloudy_bnd(iz,iens)   = 0;
  ecpp_cat_wwqui_bar_bnd(iz,iens)         = 0;
  ecpp_cat_wwqui_cloudy_bar_bnd(iz,iens)  = 0;
  ecpp_sum_wwqui_bar_bnd(iz,iens)         = 0;
  ecpp_sum_wwqui_cloudy_bar_bnd(iz,iens)  = 0;
});

parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
  ecpp_cat_area_cen_final(iPR,iTR,iCL,k,icrm) = 0;
  ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_cat_rh_cen(iPR,iTR,iCL,k,icrm)         = 0;
  ecpp_cat_qcloud_cen(iPR,iTR,iCL,k,icrm)     = 0;
  ecpp_cat_qice_cen(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_cat_precsolidcen(iPR,iTR,iCL,k,icrm)   = 0;
  ecpp_sum_area_cen_final(iPR,iTR,iCL,k,icrm) = 0;
  ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_sum_rh_cen(iPR,iTR,iCL,k,icrm)         = 0;
  ecpp_sum_qcloud_cen(iPR,iTR,iCL,k,icrm)     = 0;
  ecpp_sum_qice_cen(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_sum_precsolidcen(iPR,iTR,iCL,k,icrm)   = 0;
});

parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
  ecpp_cat_area_bnd_final(iPR,iTR,iCL,k,icrm) = 0;
  ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_cat_mass_bnd(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_sum_area_bnd_final(iPR,iTR,iCL,k,icrm) = 0;
  ecpp_sum_area_bnd(iPR,iTR,iCL,k,icrm)       = 0;
  ecpp_sum_mass_bnd(iPR,iTR,iCL,k,icrm)       = 0;
});

}

// =========================================================================|ecpp_crm_stat|=========
// allocate and initialize variables
inline void ecpp_crm_stat( pam::PamCoupler &coupler , int nstep) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");    // Note that nz   = crm_nz
  auto nzi             = coupler.get_option<int>("crm_nzi");  // Note that nzi  = crm_nz+1
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  auto ntavg1          = coupler.get_option<int>("ecpp_ntavg1");
  auto ntavg2          = coupler.get_option<int>("ecpp_ntavg2");
  auto mode_updnthresh = coupler.get_option<int>("mode_updnthresh");
  auto plumetype       = coupler.get_option<int>("plumetype");
  //auto NCLASS_CL       = coupler.get_option<int>("ecpp_NCLASS_CL");
  auto ndraft_max      = coupler.get_option<int>("ecpp_ndraft_max");
  //auto NCLASS_PR       = coupler.get_option<int>("ecpp_NCLASS_PR");
  auto gcm_dt          = coupler.get_option<real>("gcm_dt");
  auto crm_dt          = coupler.get_option<real>("crm_dt");
  int nupdraft = 0;
  int ndndraft = 0;
  int CLR = 0;       // Clear sub-class
  int CLD = 1;       // Cloudy sub-class
  int PRN = 0;       // Not precipitating sub-class
  int PRY = 1;       // Is precipitating sub-class
  int DN1 = 0; // !First index of downward classes
  int QUI       = 1; // Quiescent class
  int UP1       = 2; // First index for upward classes
  nupdraft = 1;
  ndndraft = 1;
  int ncc_in       = 2; // Nnumber of clear/cloudy sub-calsses
  int nprcp_in     = 2; // Number of non-precipitating/precipitating sub-classes.
  int NCLASS_CL = ncc_in; // Number of cloud classes
  int NCLASS_PR = nprcp_in; // Number of precipitaion classes
  int dndn = 0;
  int dnup = 0;
  dndn = ndndraft;
  dnup = nupdraft;

  printf("%s %.2f\n", "Liran check gcm_dt:",gcm_dt);
  printf("%s %.2f\n", "Liran check crm_dt:", crm_dt);
  
  printf("\nValue of ndraft_max 2: %d ", ndraft_max);
  printf("\nValue of dndn: %d: ", dndn);
  printf("\nValue of dnup: %d: ", dnup);
  printf("\nValue of NCLASS_CL: %d: ", NCLASS_CL);
  printf("\nValue of NCLASS_PR: %d: ", NCLASS_PR);

  //------------------------------------------------------------------------------------------------
  // get variables allocated in ecpp_crm_init
  auto crm_cnt                  = dm_device.get<int,1>("crm_cnt");
  auto crm_level2_cnt           = dm_device.get<int,1>("crm_level2_cnt");
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
  auto wwqui_cen                = dm_device.get<real, 2>("wwqui_cen");
  auto wwqui_bnd                = dm_device.get<real, 2>("wwqui_bnd");
  auto liran_test2d             = dm_device.get<real, 2>("liran_test2d");
  auto cldtot2d                 = dm_device.get<real, 2>("cldtot2d");
  auto wwqui_cloudy_cen         = dm_device.get<real, 2>("wwqui_cloudy_cen");
  auto wwqui_cloudy_bnd         = dm_device.get<real, 2>("wwqui_cloudy_bnd");
  auto ecpp_output_wwqui_cen         = dm_host.get<real,2>("ecpp_output_wwqui_cen");
  auto ecpp_output_wwqui_cloudy_cen  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_cen");
  auto ecpp_output_wwqui_bnd         = dm_host.get<real,2>("ecpp_output_wwqui_bnd");
  auto ecpp_output_wwqui_cloudy_bnd  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_bnd");
  auto ecpp_output_tbeg              = dm_host.get<real,2>("ecpp_output_tbeg");
  auto ecpp_output_acen              = dm_host.get<real,5>("ecpp_output_acen");
  auto ecpp_output_abnd              = dm_host.get<real,5>("ecpp_output_abnd");
  auto ecpp_output_acen_tf           = dm_host.get<real,5>("ecpp_output_acen_tf");
  auto ecpp_output_abnd_tf           = dm_host.get<real,5>("ecpp_output_abnd_tf");
  auto ecpp_output_massflxbnd        = dm_host.get<real,5>("ecpp_output_massflxbnd");
  auto ecpp_output_rhcen             = dm_host.get<real,5>("ecpp_output_rhcen");
  auto ecpp_output_qcloudcen         = dm_host.get<real,5>("ecpp_output_qcloudcen");
  auto ecpp_output_qlsinkcen         = dm_host.get<real,5>("ecpp_output_qlsinkcen");
  auto ecpp_output_precrcen          = dm_host.get<real,5>("ecpp_output_precrcen");
  auto ecpp_output_precsolidcen      = dm_host.get<real,5>("ecpp_output_precsolidcen");
  
  auto ecpp_sum_wwqui_bar_cen         = dm_device.get<real,2>("ecpp_sum_wwqui_bar_cen");
  auto ecpp_sum_wwqui_cloudy_bar_cen  = dm_device.get<real,2>("ecpp_sum_wwqui_cloudy_bar_cen");
  auto ecpp_sum_wwqui_bar_bnd         = dm_device.get<real,2>("ecpp_sum_wwqui_bar_bnd");
  auto ecpp_sum_wwqui_cloudy_bar_bnd  = dm_device.get<real,2>("ecpp_sum_wwqui_cloudy_bar_bnd");
  auto ecpp_sum_tbeg                  = dm_device.get<real,2>("ecpp_sum_tbeg");
  auto ecpp_sum_area_cen_final        = dm_device.get<real,5>("ecpp_sum_area_cen_final");
  auto ecpp_sum_area_cen              = dm_device.get<real,5>("ecpp_sum_area_cen");
  auto ecpp_sum_area_bnd_final        = dm_device.get<real,5>("ecpp_sum_area_bnd_final");
  auto ecpp_sum_area_bnd              = dm_device.get<real,5>("ecpp_sum_area_bnd");
  auto ecpp_sum_mass_bnd              = dm_device.get<real,5>("ecpp_sum_mass_bnd");
  auto ecpp_sum_rh_cen                = dm_device.get<real,5>("ecpp_sum_rh_cen");
  auto ecpp_sum_qcloud_cen            = dm_device.get<real,5>("ecpp_sum_qcloud_cen");
  auto ecpp_sum_qice_cen              = dm_device.get<real,5>("ecpp_sum_qice_cen");
  auto ecpp_sum_precsolidcen          = dm_device.get<real,5>("ecpp_sum_precsolidcen");  
  auto ecpp_cat_wwqui_bar_cen         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_cen");
  auto ecpp_cat_wwqui_bar_bnd         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_bnd");
  auto ecpp_cat_wwqui_cloudy_bar_cen  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_cen");
  auto ecpp_cat_wwqui_cloudy_bar_bnd  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_bnd");
  auto ecpp_cat_area_cen_final        = dm_device.get<real,5>("ecpp_cat_area_cen_final");
  auto ecpp_cat_area_cen              = dm_device.get<real,5>("ecpp_cat_area_cen");
  auto ecpp_cat_area_bnd_final        = dm_device.get<real,5>("ecpp_cat_area_bnd_final");
  auto ecpp_cat_area_bnd              = dm_device.get<real,5>("ecpp_cat_area_bnd");
  auto ecpp_cat_mass_bnd              = dm_device.get<real,5>("ecpp_cat_mass_bnd");
  auto ecpp_cat_rh_cen                = dm_device.get<real,5>("ecpp_cat_rh_cen");
  auto ecpp_cat_qcloud_cen            = dm_device.get<real,5>("ecpp_cat_qcloud_cen");
  auto ecpp_cat_qice_cen              = dm_device.get<real,5>("ecpp_cat_qice_cen");
  auto ecpp_cat_precsolidcen          = dm_device.get<real,5>("ecpp_cat_precsolidcen");
  auto ecpp_cat_tbeg                  = dm_device.get<real,2>("ecpp_cat_tbeg");
  // Get values from PAM cloud fields
  auto host_state_shoc_tk       = dm_host.get<real,4>("state_shoc_tk");
  auto host_state_shoc_tkh      = dm_host.get<real,4>("state_shoc_tkh");
  
  auto host_state_qv            = dm_host.get<real,4>("state_qv");
  auto host_state_qc            = dm_host.get<real,4>("state_qc");
  auto host_state_qr            = dm_host.get<real,4>("state_qr");
  auto host_state_qi            = dm_host.get<real,4>("state_qi");
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
  auto cldfrac                  = dm_device.get<real,4>( "cldfrac");
  auto host_state_shoc_tke      = dm_host.get<real,4>("state_shoc_tke");
  //------------------------------------------------------------------------------------------------
  // Define variables used by subroutine categorization_stats
  real7d mask_bnd("mask_bnd",nzi,ny,nx,nens,NCLASS_CL,ndraft_max,NCLASS_PR);
  real7d mask_cen("mask_cen",nzi,ny,nx,nens,NCLASS_CL,ndraft_max,NCLASS_PR);
  real5d area_cen_final("area_cen_final",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d area_cen("area_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d area_bnd_final("area_bnd_final",NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens);
  real5d area_bnd("area_bnd",NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens);
  real5d mass_bnd_final("mass_bnd",NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens);
  real5d mass_bnd("mass_bnd",NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens);
  real5d mass_cen_final("mass_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d mass_cen("mass_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d ent_bnd("ent_bnd",NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens);
  real5d rh_cen("rh_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qcloud_cen("qcloud_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qcloud_bf_cen("qcloud_bf_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qrain_cen("qrain_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qice_cen("qice_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qsnow_cen("qsnow_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qgraup_cen("qgraup_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qlsink_cen("qlsink_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d precr_cen("precr_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d precsolid_cen("precsolid_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d precall_cen("precall_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qlsink_bf_cen("qlsink_bf_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d qlsink_avg_cen("qlsink_avg_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);
  real5d prain_cen("prain_cen",NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens);      
  real2d acldy_cen_tbeg("acldy_cen_tbeg",nz,nens);
  real2d wwqui_bar_cen("wwqui_bar_cen",nz,nens);
  real2d wwqui_bar_bnd("wwqui_bar_bnd",nzi,nens);
  real2d wwqui_cloudy_bar_cen("wwqui_cloudy_bar_cen",nz,nens);
  real2d wwqui_cloudy_bar_bnd("wwqui_cloudy_bar_bnd",nz,nens);
  real4d cloudmixr("cloudmixr",nz,ny,nx,nens);
  real4d cloudmixr_total("cloudmixr_total",nz,ny,nx,nens);
  real4d precmixr_total("precmixr_total",nz,ny,nx,nens);
  real4d qvs("qvs",nz,ny,nx,nens);
  real4d e_vapor("e_vapor",nz,ny,nx,nens);
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
  real2d thresh_factorbb_up("thresh_factorbb_up",nz,nens);
  real2d thresh_factorbb_down("thresh_factorbb_down",nz,nens);
  int2d maskup("maskup", nzi, nupdraft);
  int2d maskdn("maskdn", nzi, nupdraft);
  int1d maskqu("maskqu", nzi);
  int1d maskcld("maskcld", nzi);
  int1d maskclr("maskclr", nzi);
  int1d maskcld_bnd("maskcld_bnd", nzi);
  int1d maskclr_bnd("maskclr_bnd", nzi);
  int1d maskpry("maskpry", nzi);
  int1d maskprn("maskprn", nzi);
  int1d maskpry_bnd("maskpry_bnd", nzi);
  int1d maskprn_bnd("maskprn_bnd", nzi);
// Liran: Threshold values needed are combined here
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
  /* Liran: module_data_ecpp1.F90 line 142
  !! subarea-average vertical mass fluxes (kg/m2/s) smaller than
   !! aw_draft_cut*rho are treated as zero
   !! note that with a*w = 1e-4 m/s, dz over 1 day = 8.6 m which is small
   */
  real aw_draft_cut = 1.0e-4;   // m/s 
  real w_draft_max = 50.0;   // m/s maximum expected updraft
  // fractional areas below afrac_cut are ignored
  real afrac_cut = aw_draft_cut/w_draft_max;
  real afrac_cut_bb  = afrac_cut*0.5;
  real afrac_cut_0p5 = afrac_cut*0.5;
  real afrac_cut_0p2 = afrac_cut*0.2;
  real afrac_cut_0p1 = afrac_cut*0.1;

  int runcount     = 0;
  int ijdel_upaa   = 0;
  int ijdel_downaa = 0;
  int ijdel_upbb   = 0;
  int ijdel_downbb = 0;
  int temp0        = 0;
  int temp1        = 0;

  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
     area_cen_final(iPR,iTR,iCL,k,icrm) = 0.0;
     area_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     mass_cen_final(iPR,iTR,iCL,k,icrm) = 0.0;
     mass_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     rh_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qcloud_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qcloud_bf_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qrain_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qice_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qsnow_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qgraup_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qlsink_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     precr_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     precsolid_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     precall_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qlsink_bf_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     qlsink_avg_cen(iPR,iTR,iCL,k,icrm) = 0.0;
     prain_cen(iPR,iTR,iCL,k,icrm) = 0.0;   
  });

  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
     ent_bnd(iPR,iTR,iCL,k,icrm) = 0.0;
     area_bnd_final(iPR,iTR,iCL,k,icrm) = 0.0;
     area_bnd(iPR,iTR,iCL,k,icrm) = 0.0;
     mass_bnd_final(iPR,iTR,iCL,k,icrm) = 0.0;
     mass_bnd(iPR,iTR,iCL,k,icrm) = 0.0;  
  });

  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    crm_cnt(iens) = crm_cnt(iens) + 1;
  });
  runcount = crm_cnt(0);
  printf("%s %d\n", "Liran check runcount:", runcount);
  // Some how if I move line 358 to 365 above before dm_device.get calls, the value ntavg1 will change. 
  real ecpp_ntavg1_ss = std::min(600.0, gcm_dt); // lesser of 10 minutes or the GCM timestep
  real ecpp_ntavg2_ss = gcm_dt;               // level-2 averaging period is GCM timestep
  // Ensure ntavg2_ss is a multiple of ntavg1_ss
  ecpp_ntavg1_ss = (int)(ecpp_ntavg2_ss / (ecpp_ntavg2_ss / ecpp_ntavg1_ss));
  temp0 = (int)(ecpp_ntavg1_ss);
  temp1 = (int)(crm_dt);
  ntavg1 = temp0/temp1;
  ntavg2 = ntavg1; // Liran: setting level 1 averaging and level 2 averaging the same for now. 
  temp0 = (int)(ecpp_ntavg2_ss);
  ntavg2 =  temp0/ temp1;
  //printf("%s %d\n", "Liran check crm_dt 00:", temp1);
  printf("%s %d\n", "Liran check int(ecpp_ntavg1_ss / crm_dt))", temp0/temp1);
  printf("%s %d\n", "Liran check ntavg1", ntavg1);
  printf("%s %d\n", "Liran check ntavg2", ntavg2);
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

esat_test = esatw_crm(T_test);
printf("%s %.2f\n", "Liran check evp:", esat_test);
parallel_for( "update sums",SimpleBounds<4>(nz, ny, nx, nens),
  YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real EVS = esatw_crm(crm_temp(k,j,i,icrm)); //   ! saturation water vapor pressure (PA)
    qvs(k,j,i,icrm) = .622*EVS/(ref_pres(icrm,k)*100.-EVS); //  ! pres(icrm,kk) with unit of hPa
    alt(k,j,i,icrm) =  287.0*crm_temp(k,j,i,icrm)/(100.*ref_pres(icrm,k));
    real rh_temp = host_state_qv(k,j,i,icrm)/qvs(k,j,i,icrm);
    altsum1(k,j,i,icrm) = altsum1(k,j,i,icrm) + alt(k,j,i,icrm);
    liran_test4davg(k,j,i,icrm) = liran_test4davg(k,j,i,icrm) + liran_test4d(k,j,i,icrm);
    qcloudsum1(k,j,i,icrm) = qcloudsum1(k,j,i,icrm) + qcloud(k,j,i,icrm);
    qrainsum1(k,j,i,icrm) = qrainsum1(k,j,i,icrm) + qrloud(k,j,i,icrm);
    qicesum1(k,j,i,icrm)  = qicesum1(k,j,i,icrm) + qiloud(k,j,i,icrm);
    prainsum1(k,j,i,icrm) = prainsum1(k,j,i,icrm) + qrloud(k,j,i,icrm);
    qsnowsum1(k,j,i,icrm) = qsnowsum1(k,j,i,icrm) + qsnowsum1(k,j,i,icrm); // This is ZERO!! for now
    ecppwwsum1(k,j,i,icrm) = ecppwwsum1(k,j,i,icrm) + crm_wvel(k,j,i,icrm);
    rhsum1(k,j,i,icrm) = rhsum1(k,j,i,icrm) + rh_temp;
});
real r_nx_ny  = 1.0/(nx*ny);  // precompute reciprocal to avoid costly divisions
parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    yakl::atomicAdd( acldy_cen_tbeg (k,iens), cldfrac(k,j,i,iens) * r_nx_ny );
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
printf("%s %.2f\n", "Liran check liran_test4davg:", liran_test4davg(10,10,10,10));
// ----------------------------------------------------------------------
// Check if we have reached the end of the level 1 time averaging period   
// ----------------------------------------------------------------------
if (runcount >=ntavg1 && runcount % ntavg1 == 0) {
  printf("%s %.2f\n", "Liran check liran_test4davg average start:", liran_test4davg(10,10,10,10));
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int icrm) {
        // itavg1 is divisible by ntavg1
      //printf("%s %d %.2f \n", "Liran check liran_test4davg 0:", crm_cnt(icrm),liran_test4d(k,j,i,icrm));
      real temp = 0.0;
      liran_test4davg(k,j,i,icrm) = liran_test4davg(k,j,i,icrm)/ntavg1;
      qcloudsum1(k,j,i,icrm)      = qcloudsum1(k,j,i,icrm)     /ntavg1;
      qrainsum1(k,j,i,icrm)       = qrainsum1(k,j,i,icrm)      /ntavg1;
      qicesum1(k,j,i,icrm)        = qicesum1(k,j,i,icrm)       /ntavg1;
      qsnowsum1(k,j,i,icrm)       = qsnowsum1(k,j,i,icrm)      /ntavg1;
      altsum1(k,j,i,icrm)         = altsum1(k,j,i,icrm)        /ntavg1;
      ecppwwsum1(k,j,i,icrm)      = ecppwwsum1(k,j,i,icrm)     /ntavg1;
      rhsum1(k,j,i,icrm)          = rhsum1(k,j,i,icrm)         /ntavg1;
      //printf("%s %d %d %d %d %d %.2f \n", "Liran check liran_test4davg 1:", crm_cnt(icrm),k,j,i,icrm,liran_test4d(k,j,i,icrm));
  });

  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    acldy_cen_tbeg (k,iens) = acldy_cen_tbeg (k,iens)  /ntavg1;
  });
  printf("%s %.2f\n", "Liran check liran_test4davg average end:", liran_test4davg(10,10,10,10));
  printf("Liran check categorization_stats\n");
  //printf("%s %d\n", "Liran check itavg1 div ntavg1:", itavg1 % ntavg1);
  


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

  int nxy = nx*ny;
  parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int icrm) {
    for (int k_gcm=0; k_gcm<nz; k_gcm++) {
      //int l = plev-(k+1);
      int k_crm= (gcm_nlev+1)-1-k_gcm;
      /*
      ! Get cloud top height
      ! Cloud top height is used to determine whether there is updraft/downdraft. No updraft and
      ! downdraft is allowed above the condensate level (both liquid and ice).
      */
      cloudtop(j,i,icrm) = 1; // !Default to bottom level if no cloud in column.
      // BELOW: 0.01*qvs may be too large at low level.
      // if( cloudmixr_total(i,j,k) >= max(0.01*qvs(i,j,k), cloudthresh_trans) ) then
      if (cloudmixr_total(k_crm,j,i,icrm) >= cloudthresh_trans) {
        cloudtop(j,i,icrm) = k_crm; 
        break; // exit the loop
      }
    }
    for (int k_crm=0; k_crm<nzi; k_crm++) {
      int km0 = std::min(nz,k_crm);
      int km1 = std::max(1,k_crm-1);
      rhoair(k_crm,icrm) = rhoair(k_crm,icrm)+0.5*(1.0/altsum1(km1,j,i,icrm) + 1.0/altsum1(km0,j,i,icrm))/nxy;
    }
  }); // end of parallel_for( SimpleBounds<3>(ny,nx,nens)

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
    cloudtop_upbb(j,i,icrm)   = nz-1;
    cloudtop_downbb(j,i,icrm) = nz-1;
  });

  /*
      ! Get standard deviation of up and down vertical velocity below the
      ! cloud tops. For now, each cell is treated equally. We may want to
      ! consider weighting each cell by its volume or mass.
      !
      ! Get the mean values first for wup and wdown
  */
  printf("Liran check categorization_stats 3\n");
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

  }); // end of parallel_for( SimpleBounds<3>(ny,nx,nens)
  printf("Liran check categorization_stats 4\n");

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
  printf("Liran check categorization_stats 5\n");
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
  printf("Liran check categorization_stats 6\n");
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
      int kmin = 0;
      kmin = std::max(k_crm-1, 1);
      int kmax = 0;
      kmax = std::min(k_crm+1, nz);
      for (int k_crm2=kmin; k_crm<=kmax; k_crm++) {
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
  printf("Liran check categorization_stats 7\n");
  /*
    ! Get masks to determine (cloud vs. clear) (up vs. down vs. other) categories.
    ! Vertical velocities are checked on the cell vertical interfaces to determine
    ! if they pass the threshold criteria. Clouds below the interface are then
    ! used for updrafts and above the int. for downdrafts. Quiescent (other)
    ! drafts use an average of the cloud above and below the interface to
    ! determine cloudiness.

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
  printf("Liran check categorization_stats 8\n");
  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    real tmpsuma = 0.0;
    real tmpw = 0.0;
    real tmpw_minval = 0.10;
    real tmpsumb = 1.0e-30;
    for (int k_crm=0; k_crm<=nz; k_crm++) {
      tmpw = wdown_rms_k(k_crm,icrm);
      tmpw = std::max(1.0e-4,tmpw);
      tmpw = tmpw * rhoair(k_crm,icrm);  // weighted by mass?
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
      wdown_thresh_k(k_crm,0,icrm) = -tmpw*std::abs(downthresh);
      wdown_thresh_k(k_crm,1,icrm) = -tmpw*std::abs(downthresh2);
    }

  });
  printf("Liran check categorization_stats 9\n");
  // End of call determine_transport_thresh

  // Starting line 608
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int icrm) {
    wdown_thresh(k_crm,icrm) = wdown_thresh_k(k_crm,0,icrm);
    wup_thresh(k_crm,icrm) = wup_thresh_k(k_crm,0,icrm);
  });

  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int icrm) {
    thresh_factorbb_up(k_crm,icrm) = 1.0;
    thresh_factorbb_down(k_crm,icrm) = 1.0;
  });
  
  
  printf("Liran check categorization_stats 10\n");

  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    bool thresh_calc_not_done = true;
    int iter = 0;
    int nzstag = nz + 1;

    while (thresh_calc_not_done) {
      iter = iter + 1;
    
      for (int k_crm = 1; k_crm <= nzstag; ++k_crm) {
          real tmpa, tmpb;
          if (k_crm == 1) {
              tmpa = thresh_factorbb_up(k_crm,icrm); 
              tmpb = thresh_factorbb_down(k_crm,icrm);
          } else if (k_crm == nzstag) {
              tmpa = thresh_factorbb_up(k_crm-1,icrm); 
              tmpb = thresh_factorbb_down(k_crm-1,icrm);
          } else {
              tmpa = std::max(thresh_factorbb_up(k_crm-1,icrm), thresh_factorbb_up(k_crm,icrm));
              tmpb = std::max(thresh_factorbb_down(k_crm-1,icrm), thresh_factorbb_down(k_crm,icrm));
          } // end of if (k_crm == 1) 
          wup_thresh_k(k_crm,0,icrm) = wup_thresh_k(k_crm,0,icrm)*tmpa;
          wup_thresh_k(k_crm,1,icrm) = wup_thresh_k(k_crm,1,icrm)*tmpa;
          wdown_thresh_k(k_crm,0,icrm) = wdown_thresh_k(k_crm,0,icrm)*tmpb;
          wdown_thresh_k(k_crm,1,icrm) = wdown_thresh_k(k_crm,1,icrm)*tmpb;
      } // end of for (int k_crm = 1; k_crm <= nzstag; ++k_crm)


      for (int k_crm = 0; k_crm <= std::max(1, kup_top(icrm)-1); ++k_crm) {
        wup_thresh(k_crm,icrm) = wup_thresh_k(k_crm,0,icrm);
      }
      for (int k_crm = 0; k_crm <= std::max(1, kup_top(icrm)-1); ++k_crm) {
        wdown_thresh(k_crm,icrm) = wdown_thresh_k(k_crm,0,icrm);
      }

      /*
      !
      !  fix a bug in the WRF_ECPP, Minghuai Wang, 2009-12.
      !  set wdown_thresh_k and wup_thresh_k to be an extreme value
      !  above updraft (kup_top) and downdraft top(kdown_top).
      !  This will make sure there is no updraft or downdraft above kup_top and kdown_top
      !
      */

      real wlarge = 1.0e10;   // m/s   Liran note: This seems important to adjust!!

      for (int k_crm = kup_top(icrm); k_crm <= nz; ++k_crm) {
          wup_thresh_k(k_crm, 0,icrm) = wlarge;
          wup_thresh_k(k_crm, 1,icrm) = wlarge;
      }
      for (int k_crm = kdown_top(icrm); k_crm <= nz; ++k_crm) {
          wdown_thresh_k(k_crm, 0,icrm) =  -1. * wlarge;
          wdown_thresh_k(k_crm, 1,icrm) =  -1. * wlarge;
      }


  /* ==========================================================================
  subroutine setup_class_masks( &
      nx, ny, nz, nupdraft, ndndraft, ndraft_max, &
      cloudmixr, cf3d, precall, ww, &
      wdown_thresh_k, wup_thresh_k, &
      cloudthresh, prcpthresh, &
      mask_bnd, mask_cen,  &
      cloudmixr_total, cloudthresh_trans, precthresh_trans, &
      qvs, precmixr_total )
      !
      ! Sets up the masks used for determining quiescent/up/down, clear/cloudy,
      ! and non-precipitatin/precipitating classes.
      !
      ! William.Gustafosn@pnl.gov; 20-Nov-2008
      ! Last modified: William.Gustafson@pnl.gov; 16-Apr-2009

      ! Modification by Minghuai Wang (Minghuai.Wang@pnl.gov), April 23, 2010
      ! use total condensate (liquid+ice),  different condensate and precipitating thresholds
      ! to classify transport classes.
      ! See Xu et al., 2002, Q.J.R.M.S.
      !
  */
  printf("Liran check categorization_stats 11\n");

  /*
      ! Initialize the masks to zero and then we will accumulate values into
      ! them as we identify the various classes.
  */

      parallel_for(SimpleBounds<7>(nzi,ny,nx,nens,NCLASS_CL,ndraft_max,NCLASS_PR), YAKL_LAMBDA (int k,int j,int i,int icrm,int icl,int itr,int ipr) {
        mask_bnd(k,j,i,icrm,icl,itr,ipr)           = 0;
        mask_cen(k,j,i,icrm,icl,itr,ipr)           = 0;
      });

      real cloudthresh_trans_temp=0.0;
      real precthresh_trans_temp=0.0;
      int itr;
      int ipr;
      int icl;
      // Loop over the horizontal dimensions...
      for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
        // Set initial mask values for the vertical cell boundaries...
          for (int k=0; k<nzstag; k++) {
            maskup(k,0) =  0;
            maskdn(k,0) =  0;
            maskqu(k) =  0;
            maskcld(k) =  0;
            maskclr(k) =  0;
            maskcld_bnd(k) =  0;
            maskclr_bnd(k) =  0;
            maskpry(k) =  0;
            maskprn(k) =  0;
            maskpry_bnd(k) =  0;
            maskprn_bnd(k) =  0;
  /*        
            !Transport upward at cell boundaries...
            !We have to take into account the possibility of multiple
            !updraft categories. At this point, we handle only the
            !cases of one or two categories. We do not yet handle the
            !allcomb option.
            !
            ! updraft only exist in cloudy area or precipitating clear area ++++mhwang
  */          
            cloudthresh_trans_temp = cloudthresh_trans;
            if ( (cloudmixr_total(std::max(k - 1, 0), j, i, icrm) + cloudmixr_total(std::min(k, nz-1), j, i, icrm)) * 0.5 > cloudthresh_trans_temp ||
               (precmixr_total(std::max(k - 1, 0),j,i,icrm)  + precmixr_total(std::min(k, nz-1),j,i,icrm)) * 0.5 > precthresh_trans ){
               // Liran Only one threshold
              if (crm_wvel(k,j,i,icrm) > wup_thresh_k(k, 0,icrm)) {
                  maskup(k,0) = 1;
              } 
            } // end of if (((cloudmixr_total(...

/*
        !Transport downward at cell boundaries...
        !
        ! downdraft only exist in cloudy area or precipitating clear area   +++mhwang
*/
            if ( (cloudmixr_total(std::max(k - 1, 0), j, i, icrm) + cloudmixr_total(std::min(k, nz-1), j, i, icrm)) * 0.5 > cloudthresh_trans_temp ||
                (precmixr_total(std::max(k - 1, 0),j,i,icrm)  + precmixr_total(std::min(k, nz-1),j,i,icrm)) * 0.5 > precthresh_trans ){
               // !Only one threshold
              if (crm_wvel(k,j,i,icrm) < wdown_thresh_k(k, 0,icrm)) {
                  maskdn(k,0) = 1;
              } 
            } // end of if (((cloudmixr_total(...
/*
      !Transport quiescent at cell boundaries if neither up or
      !down triggered...
*/
            if ( maskup(k,0) + maskdn(k,0) < 1 ){
              maskqu(k) = 1;
            } 

            // Cloudy or clear at cell boundaries...
            if ((cloudmixr(std::max(k - 1, 0), j, i, icrm)+cloudmixr(std::min(k, nz-1), j, i, icrm))*0.5>cloudthresh){
              maskcld_bnd(k) = 1;
            } else {
              maskclr_bnd(k) = 1;
            }

            // Raining or not at cell boundaries...
            if ((precmixr_total(std::max(k - 1, 0), j, i, icrm)+precmixr_total(std::min(k, nz-1), j, i, icrm))*0.5>prcpthresh){
              maskpry_bnd(k) = 1;
            } else {
              maskprn_bnd(k) = 1;
            }
          } // end of for (int k=0; k<=nzstag; k++) 

          for (int k=0; k<=nz; k++) {
            // Cloudy or clear at cell centers...
            if (cloudmixr_total(k,j,i,icrm)>cloudthresh){
              maskcld(k) = 1;
            } else {
              maskclr(k) = 1;
            }
            // Raining or not at cell centers...
            if (precmixr_total(k,j,i,icrm)>prcpthresh){
              maskpry(k) = 1;
            } else {
              maskprn(k) = 1;
            }
          } // end of for (int k=0; k<=nz; k++) {

/*
    ! Now, use the initial boundary masks by class to generate a combined
    ! mask for the cell boundaries.
*/
          for (int k=0; k<nzstag; k++) {
            // Upward, or at least upward quiescent
            // Liran nupdraft = 1 for now so no need to estimate sum(maskup(k,:))
            if (maskup(k,0)>0 || (maskqu(k)>0 && crm_wvel(k,j,i,icrm)>0)){

              // !Are we are here because of maskup? If so, then we need to
              // !parse the correct updraft category.
              if (maskqu(k) < 1){
                itr = UP1; //+ maskup(k,0)-1;   
                // Liran: To exclude the condition maskup(k,0)=0 [maskdn(k)=1 and crm_wvel(k,j,i,icrm)>0] = [maskqu(k)=1]  
              } else {
                itr = QUI; // Liran: itr = QUI even maskup(k,0)>0
              }
              // For upward motion, determine cloud and precip characteristics
              // based on the cell-center values below the boundary.
              if (k==0){
                icl = CLR;
                ipr = PRN;
              } else {
                //setup_class_masks: bnd cloud up
                if (maskcld(k-1) > 0 && maskclr(k-1) < 1) {
                    icl = CLD;
                } else if (maskclr(k-1) > 0 && maskcld(k-1) < 1) {
                    icl = CLR;
                } 
                //setup_class_masks: bnd prcp up
                if (maskpry(k-1) > 0 && maskprn(k-1) < 1) {
                    ipr = PRY;
                } else if (maskprn(k-1) > 0 && maskpry(k-1) < 1) {
                    ipr = PRN;
                } 
                // call cloud_prcp_check
                //printf("\nsetup_class_masks: bnd prcp up %d %d: ", k-1,icl);
                //printf("\nsetup_class_masks: bnd prcp up %d %d: ", k-1,ipr);
              }
            // Downward, or at least downward quiescent  
            } else if (maskdn(k,0)>0 || (maskqu(k)>0 && crm_wvel(k,j,i,icrm)<0)){  
            
              // Are we here because of maskdn? If so, then we need to
              // parse the correct downdraft category.  
              if (maskqu(k) < 1){
                itr = DN1; //+ maskup(k,0)-1;   
              } else {
                itr = QUI;
              }
              // For downward motion, determine cloud and precip characteristics
              // based on the cell-center values above the boundary.
              if( k==nzstag-1 ) {
                icl = CLR;
                ipr = PRN;
              }else{
                // call cloud_prcp_check
                // setup_class_masks: bnd cloud down
                if (maskcld(k) > 0 && maskclr(k) < 1) {
                    icl = CLD;
                } else if (maskclr(k) > 0 && maskcld(k) < 1) {
                    icl = CLR;
                } 
                // setup_class_masks: bnd prcp down
                if (maskpry(k) > 0 && maskprn(k) < 1) {
                    ipr = PRY;
                } else if (maskprn(k) > 0 && maskpry(k) < 1) {
                    ipr = PRN;
                } 
              }
            // Quiescent with w=0. Use the cell-center values averaged
            // surrounding the boundary for the cloud/prcp states.  
            } else {
              itr = QUI;
              // call cloud_prcp_check
              // setup_class_masks: bnd cloud quiescent
              if (maskcld_bnd(k) > 0 && maskclr_bnd(k) < 1) {
                  icl = CLD;
              } else if (maskclr_bnd(k) > 0 && maskcld_bnd(k) < 1) {
                  icl = CLR;
              } 
              // setup_class_masks: bnd prcp quiescent
              if (maskpry_bnd(k) > 0 && maskprn_bnd(k) < 1) {
                  ipr = PRY;
              } else if (maskprn_bnd(k) > 0 && maskpry_bnd(k) < 1) {
                  ipr = PRN;
              } 
            } // end of if (maskup(k,0)>0 || (maskqu(k)>0 && crm_wvel(k,j,i,icrm)>0))
/*
! +++mhwang Line 1612
    ! Total condensate and different thresholds are used to classify transport classes. So the following change
    ! is not needed anymore. Minghuai Wang, 2010-04-23.
    !
    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have all the class indices determined so now we can set
    !the correct mask location to 1.
    !           mask_bnd(i,j,k,icl,itr,ipr) = 1.
    ! use fractioal cloudiness in SAM
*/
            if( icl==CLR ) {
              mask_bnd(k,j,i,icrm,icl,itr,ipr) = 1;
            }else if (icl==CLD){
              mask_bnd(k,j,i,icrm,CLD,itr,ipr) = (cldfrac(std::max(k - 1, 0), j, i, icrm) + cldfrac(std::min(k, nz-1), j, i, icrm)) * 0.5; 
              mask_bnd(k,j,i,icrm,CLR,itr,ipr) = 1.0 - (cldfrac(std::max(k - 1, 0), j, i, icrm) + cldfrac(std::min(k, nz-1), j, i, icrm)) * 0.5; 
            }
            //printf("\nmask_bnd: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,icl,itr,ipr,mask_bnd(k,j,i,icrm,icl,itr,ipr));
          } // end of for (int k=0; k<=nzstag; k++)
/*
  ! Now, use the initial boundary masks by class to generate a combined
  ! mask for the cell centers. We determine the transport class based on
  ! splitting the cell conceptually in half with the upper boundary
  ! influencing the top half of the cell and the bottom boundary the bottom
  ! half. Each contributes either 0 or 0.5 of the total contribution of the
  ! cell's transport. e.g. if both boundaries are upward, then the cell is
  ! fully an "up" transport cell. If the two boundaries are opposite, then
  ! the cell is weighted half in each direction for the masking.
  !
*/
          for (int k=0; k<nz; k++) {
            // Get the cloud/prcp characteristics at cell center.
            // call cloud_prcp_check
            if (maskcld(k) > 0 && maskclr(k) < 1) {
                icl = CLD;
            } else if (maskclr(k) > 0 && maskcld(k) < 1) {
                icl = CLR;
            } 
            if (maskpry(k) > 0 && maskprn(k) < 1) {
                ipr = PRY;
            } else if (maskprn(k) > 0 && maskpry(k) < 1) {
                ipr = PRN;
            } 
            // Look at the bottom boundary first and determine it's
            // contribution to the cell center transport class.
            if (maskup(k,0)>0){
              itr = UP1 + maskup(k,0)-1; 
            } else if (maskdn(k,0)>0){  
              itr = DN1 + maskdn(k,0)-1;   
            } else if (maskqu(k) > 0){
              itr = QUI;
            } else {
              // call stop
            }
/*
    ! +++mhwang
    ! ! Total condensate and different thresholds are used to classify transport classes. So the following change
    ! is not needed anymore. Minghuai Wang, 2010-04-23.

    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have what we need for the cell bottom classes so increment
    !the center mask for the bottom half...
    !           mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    ! Use fractional cloudiness at SAM
*/
            if (icl==CLR){
              //printf("\nmask_cen 0: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,icl,itr,ipr,mask_cen(k,j,i,icrm,icl,itr,ipr));
              mask_cen(k,j,i,icrm,icl,itr,ipr) =  mask_cen(k,j,i,icrm,icl,itr,ipr) + 0.5; 
              //printf("\nmask_cen 1: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,icl,itr,ipr,mask_cen(k,j,i,icrm,icl,itr,ipr));

            } else if (icl==CLD){  
              mask_cen(k,j,i,icrm,CLD,itr,ipr) =  mask_cen(k,j,i,icrm,CLD,itr,ipr) + cldfrac(k, j, i, icrm) * 0.5; 
              mask_cen(k,j,i,icrm,CLR,itr,ipr) =  mask_cen(k,j,i,icrm,CLR,itr,ipr) + (1.0-cldfrac(k, j, i, icrm))*0.5;  
              //printf("\nmask_cen 20: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,CLD,itr,ipr,mask_cen(k,j,i,icrm,CLD,itr,ipr));
              //printf("\nmask_cen 21: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,CLR,itr,ipr,mask_cen(k,j,i,icrm,CLR,itr,ipr));
            }

            // !Next, look at the top boundary and determine it's
            // !contribution to the cell center transport class.
            // Liran In C++ convert 1 index to 0 index, so k+1 becomes k.
            if (maskup(k,0)>0){
              itr = UP1 + maskup(k,0)-1; 
            } else if (maskdn(k,0)>0){  
              itr = DN1 + maskdn(k,0)-1;   
            } else if (maskqu(k) > 0){
              itr = QUI;
            } else {
              // call stop
            }
/* 
    ! +++mhwang
    ! In the clear, and non-precipitating class, it is classified as quiescent class in the MMF simulation.
    ! If this is classed as updraft or downdraft in mode 16, this would lead to too much upraft and downdraft mass fluxes.
    ! Minghuai Wang, 2010-01-18 (Minghuai.Wang@pnl.gov)
    !           if(icl.eq.CLR .and. ipr.eq.PRN) then
    !             itr = QUI
    !           end if
    !---mhwang

    !We have what we need for the cell top classes so increment
    !the center mask for the top half...
    !           mask_cen(i,j,k,icl,itr,ipr) = mask_cen(i,j,k,icl,itr,ipr) + 0.5
    ! use fractional cloudiness in SAM
*/
            if (icl==CLR){
              //printf("\nmask_cen 03: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,icl,itr,ipr,mask_cen(k,j,i,icrm,icl,itr,ipr));
              mask_cen(k,j,i,icrm,icl,itr,ipr) =  mask_cen(k,j,i,icrm,icl,itr,ipr) + 0.5; 
              //printf("\nmask_cen 3: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,icl,itr,ipr,mask_cen(k,j,i,icrm,icl,itr,ipr));

            } else if (icl==CLD){  
              mask_cen(k,j,i,icrm,CLD,itr,ipr) =  mask_cen(k,j,i,icrm,CLD,itr,ipr) + cldfrac(k, j, i, icrm) * 0.5; 
              mask_cen(k,j,i,icrm,CLR,itr,ipr) =  mask_cen(k,j,i,icrm,CLR,itr,ipr) + (1.0-cldfrac(k, j, i, icrm))*0.5;  
              //printf("\nmask_cen 4: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,CLD,itr,ipr,mask_cen(k,j,i,icrm,CLD,itr,ipr));
              //printf("\nmask_cen 5: %d %d %d %d %d %d %d %.2f: ", k,j,i,icrm,CLR,itr,ipr,mask_cen(k,j,i,icrm,CLR,itr,ipr));
            }

          } // end of for (int k=0; k<=nz; k++) 

        } // end of for (int j=0; y_crm<=ny; y_crm++)
      } // end of for (int i=0; x_crm<=nx; x_crm++)
      // Other loops over k, as in the original code
  
// end of setup_class_masks
/* ==========================================================================
      !
      ! ( code added on 14-dec-2009 to guarantee quiescent class
      !    area > acen_quiesc_minaa )
      ! at each level
      !    calculate total fractional area for quiescent class
      !       using the current level-1 averages
      !    if (acen_quiesc < acen_quiesc_minaa), increase the
      !       thresh_factorbb_up/down(k) by factor of 1.5 or 1.2
      !    (also, if acen_down > acen_up, increase thresh_factorbb_up by less
      !
*/

      real acen_quiesc;
      real acen_up;
      real acen_down;
      real abnd_quiesc;
      real abnd_up;
      real abnd_down;
      real acen_quiesc_minaa;

      thresh_calc_not_done = false;
      // Liran: This code is to calculate the sum
      for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
          for (int i=0; i<nx; i++) {
            for (int iCL=0; iCL<NCLASS_CL; iCL++) {
              for (int iPR=0; iPR<NCLASS_PR; iPR++) {
                acen_quiesc = acen_quiesc + mask_cen(k,j,i,icrm,iCL,QUI,iPR);
                acen_up = acen_up + mask_cen(k,j,i,icrm,iCL,UP1,iPR); 
                abnd_quiesc = abnd_quiesc + mask_bnd(k,j,i,icrm,iCL,QUI,iPR); 
                abnd_up = abnd_up + mask_bnd(k,j,i,icrm,iCL,UP1,iPR); 
                //printf("\nmask_cen8: %d %d %d %d %d %d %d %.8f %.8f : ", k,j,i,icrm,iCL,QUI,iPR,mask_cen(k,j,i,icrm,iCL,QUI,iPR),mask_bnd(k,j,i,icrm,iCL,QUI,iPR));
              } // end of for (int iPR=0
            } // end of for (int iCL=0
          } // end of for (int i=0
        } // end of for (int j=0
        acen_quiesc = std::max(acen_quiesc/nxy, 0.0);
        acen_up = std::max(acen_up/nxy, 0.0);
        acen_down = std::max(1.0-acen_quiesc-acen_up, 0.0);
        abnd_quiesc = std::max(abnd_quiesc/nxy, 0.0);
        abnd_up = std::max(abnd_up/nxy, 0.0);
        abnd_down = std::max( (1.0 - abnd_quiesc - abnd_up),0.0);

        real tmpa;

        if (std::min(acen_quiesc, abnd_quiesc)<acen_quiesc_minaa){
          thresh_calc_not_done = true;
          if (acen_down > acen_up){
            tmpa = acen_up/acen_down; 
          } else if (abnd_down > abnd_up){  
            tmpa = abnd_up/abnd_down;  
          } else {
            tmpa = 1.0;
          }

          if (std::min(acen_quiesc,abnd_quiesc) < 0.5*acen_quiesc_minaa){
            thresh_factorbb_down(k,icrm) = thresh_factorbb_down(k,icrm)*1.5; 
            thresh_factorbb_up(k,icrm) = thresh_factorbb_up(k,icrm)*std::max(1.5*tmpa,1.25);
          } else {
            thresh_factorbb_down(k,icrm) = thresh_factorbb_down(k,icrm)*1.25;
            thresh_factorbb_up(k,icrm) = thresh_factorbb_up(k,icrm)*std::max(1.25*tmpa, 1.125);
          }

          if (iter>5){
            printf("\nLiran iter greater than 5\n");
          }
        }
      } // end of for (int k=0; k<=nz; k++)
    } //while (thresh_calc_not_done
    real mask_tmp;
    int km0;
    int km1;
    int km2;
    real testgt0 = 0.0;
    real wwrho_k = 0.0;
    real wwrho_km1 = 0.0;
    real tempwork = 0.0;

    //printf("\ncheck const: %d %d %d %d %d %d %d: ", nxy,QUI,CLD,CLR,NCLASS_CL,ndraft_max,NCLASS_PR);
    // nxy = 8, QUI = 1, CLD = 2, CLR = 1, NCLASS_CL = 2, ndraft_max = 3, NCLASS_PR = 2
    // iCL,QUI,iPR
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        for (int iCL=0; iCL<NCLASS_CL; iCL++) {
          for (int iTR=0; iTR<ndraft_max; iTR++) {
            for (int iPR=0; iPR<NCLASS_PR; iPR++) {
              for (int k=0; k<nz; k++) {
                // We now have enough information to aggregate the variables into domain
                // averages by class. Do this first for the cell centers...
                //if (iCL==CLR){
                //  printf("\narea_cen0: %d %d %d %d %d %d %d %.8f %.8f : ", k,j,i,icrm,iCL,iTR,iPR,mask_cen(k,j,i,icrm,iCL,iTR,iPR),area_cen(iPR,iTR,iCL,k,icrm));
                //}
                mask_tmp = mask_cen(k,j,i,icrm,iCL,iTR,iPR)/nxy;
                area_cen_final(iPR,iTR,iCL,k,icrm) = area_cen_final(iPR,iTR,iCL,k,icrm) + mask_tmp;
                area_cen(iPR,iTR,iCL,k,icrm) = area_cen(iPR,iTR,iCL,k,icrm) + mask_tmp;
                //if (iCL==CLR){
                //  printf("\narea_cen : %d %d %d %d %d %d %d %.8f %.8f %.8f: ", k,j,i,icrm,iCL,iTR,iPR,mask_tmp,mask_cen(k,j,i,icrm,iCL,iTR,iPR),area_cen(iPR,iTR,iCL,k,icrm));
                //}
                rh_cen(iPR,iTR,iCL,k,icrm) = rh_cen(iPR,iTR,iCL,k,icrm) + rhsum1(k,j,i,icrm) *mask_tmp;
                qcloud_cen(iPR,iTR,iCL,k,icrm) = qcloud_cen(iPR,iTR,iCL,k,icrm) + qcloud(k,j,i,icrm)*mask_tmp;
                qrain_cen(iPR,iTR,iCL,k,icrm) = qrain_cen(iPR,iTR,iCL,k,icrm) + qrloud(k,j,i,icrm)*mask_tmp;
                qice_cen(iPR,iTR,iCL,k,icrm) = qice_cen(iPR,iTR,iCL,k,icrm) + qiloud(k,j,i,icrm)*mask_tmp;
                // This list is not complete! 
                // ! calculate the mean vertical velocity over the quiescent class  +++mhwang
                if(iTR==QUI){
                  wwqui_bar_cen(k,icrm) = wwqui_bar_cen(k,icrm)+(crm_wvel(k,j,i,icrm)+crm_wvel(k+1,j,i,icrm))*0.5*mask_tmp;
                  if(iCL==CLD){
                    wwqui_cloudy_bar_cen(k,icrm)=wwqui_cloudy_bar_cen(k,icrm)+(crm_wvel(k,j,i,icrm)+crm_wvel(k+1,j,i,icrm))*0.5*mask_tmp;
                  } // if(icl==CLD)
                } // if(itr==QUI)
              } // end of for (int k=0
              // Now, we can do a similar aggregation for the cell boundaries.
              for (int k = 0; k < nzstag; ++k) {
                mask_tmp = mask_bnd(k,j,i,icrm,iCL,iTR,iPR)/nxy;
                area_bnd_final(iPR,iTR,iCL,k,icrm) = area_bnd_final(iPR,iTR,iCL,k,icrm) + mask_tmp;
                // NOTE: technically we should interpolate and not do a simple
                //       average to get density at the cell interface
                km0 = std::min(nz,k);   // Liran: should we change k to k-1 in c++?
                km1 = std::max(0,k-1);
                km2 = std::max(0,k-2);

                wwrho_k   = 0.5*(1.0/altsum1(km1,j,i,icrm) + 1.0/altsum1(km0,j,i,icrm))*crm_wvel(k,j,i,icrm);
                wwrho_km1 = 0.5*(1.0/altsum1(km2,j,i,icrm) + 1.0/altsum1(km1,j,i,icrm))*crm_wvel(km1,j,i,icrm);
                testgt0 = std::max(0.0,wwrho_k-wwrho_km1);
                //if (iCL==CLR){
                //  printf("\narea_bnd 0: %d %d %d %d %d %.8f %.8f : ", iPR,iTR,iCL,k,icrm,mask_tmp,area_bnd(iPR,iTR,iCL,k,icrm));
                //}
                area_bnd(iPR,iTR,iCL,k,icrm) = area_bnd(iPR,iTR,iCL,k,icrm) + mask_tmp;
                //if (iCL==CLR){
                //  printf("\narea_bnd 1: %d %d %d %d %d %.8f %.8f : ", iPR,iTR,iCL,k,icrm,mask_tmp,area_bnd(iPR,iTR,iCL,k,icrm));
                //}
                mass_bnd(iPR,iTR,iCL,k,icrm) = mass_bnd(iPR,iTR,iCL,k,icrm) + wwrho_k*mask_tmp;
                ent_bnd(iPR,iTR,iCL,k,icrm) = ent_bnd(iPR,iTR,iCL,k,icrm) + testgt0*mask_tmp;
                // ! calculate the mean vertical velocity over the quiescent class  +++mhwang
                if(iTR==QUI){
                  wwqui_bar_bnd(k,icrm) = wwqui_bar_bnd(k,icrm)+(crm_wvel(k,j,i,icrm))*mask_tmp;
                  if(iCL==CLD){
                    wwqui_cloudy_bar_bnd(k,icrm)=wwqui_cloudy_bar_bnd(k,icrm)+(crm_wvel(k,j,i,icrm))*mask_tmp;
                  } // if(icl==CLD)
                } // if(itr==QUI)
              } // end of for (int k = 1; k <= nzstag; ++k)
            } // end of for (int iPR=0
          } // end of for (int iND=0
        } // end of for (int iCL=0
      } // end of for (int i=0
    } // end of for (int j=0


// ! calcualte vertical velocity variance for quiescent class (total and cloudy)  +++mhwang


    real abnd_up;
    real abnd_down;
    real acen_quiesc_minaa;

    for (int k=0; k<nz; k++) {
      real sum_mask_cen_CL = 0.0;
      real sum_mask_bnd_CL = 0.0;
      real sum_mask_cen_CLD = 0.0;
      real sum_mask_bnd_CLD = 0.0;
      for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
          for (int iPR=0; iPR<NCLASS_PR; iPR++) {
            for (int iCL=0; iCL<NCLASS_CL; iCL++) {
              sum_mask_cen_CL = sum_mask_cen_CL + mask_cen(k,j,i,icrm,iCL,QUI,iPR);
              sum_mask_bnd_CL = sum_mask_bnd_CL + mask_bnd(k,j,i,icrm,iCL,QUI,iPR);
            } // end of or (int iCL=0
            sum_mask_cen_CLD = sum_mask_cen_CLD + mask_cen(k,j,i,icrm,CLD,QUI,iPR);
            sum_mask_bnd_CLD = sum_mask_bnd_CLD + mask_bnd(k,j,i,icrm,CLD,QUI,iPR);
          } // end of for (int iPR=0
        } // end of for (int i=0
      } // end of for (int j=0
    
      if (sum_mask_cen_CL>0.5){
        wwqui_bar_cen(k,icrm) = wwqui_bar_cen(k,icrm) *nxy/sum_mask_cen_CL;
      } else {
        wwqui_bar_cen(k,icrm) = 0.0;
      } // end of if (sum_mask_cen_CL>0.5)

      if (sum_mask_cen_CLD>0.5){
        wwqui_cloudy_bar_cen(k,icrm) = wwqui_cloudy_bar_cen(k,icrm) *nxy/sum_mask_cen_CL;
      } else {
        wwqui_cloudy_bar_cen(k,icrm) = 0.0;
      } // end of if (sum_mask_cen_CLD>0.5)

      if (sum_mask_bnd_CL>0.5){
        wwqui_bar_bnd(k,icrm) = wwqui_bar_bnd(k,icrm) *nxy/sum_mask_bnd_CL;
      } else {
        wwqui_bar_bnd(k,icrm) = 0.0;
      } // end of if (sum_mask_bnd_CL>0.5)

      if (sum_mask_bnd_CLD>0.5){
        wwqui_cloudy_bar_bnd(k,icrm) = wwqui_cloudy_bar_bnd(k,icrm) *nxy/sum_mask_cen_CL;
      } else {
        wwqui_cloudy_bar_bnd(k,icrm) = 0.0;
      } // end of if (sum_mask_bnd_CLD>0.5)

    } // end of for (int k=0


    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        for (int iCL=0; iCL<NCLASS_CL; iCL++) {
          for (int iPR=0; iPR<NCLASS_PR; iPR++) {
            for (int k=0; k<nz; k++) {
              mask_tmp = mask_cen(k,j,i,icrm,iCL,QUI,iPR)/nxy;
/*
              ! calculate the vertical velocity variance over the quiescent class  +++mhwang
              ! wwqui_bar_cen is used in for both all sky and cloudy sky.
              ! when wwqui_cloudy_bar_cen was used for cloudy sky, wwqui_cloudy_cen will be smaller than wwqui_cen.
              !
*/
              tempwork = ((crm_wvel(k,j,i,icrm)+crm_wvel(k+1,j,i,icrm))*0.5-wwqui_bar_cen(k,icrm));
              wwqui_cen(k,icrm) = wwqui_cen(k,icrm)+mask_tmp * tempwork*tempwork + mask_tmp * host_state_shoc_tke(k,j,i,icrm)/3.0;
              if(iCL==CLD){
                tempwork = ((crm_wvel(k,j,i,icrm)+crm_wvel(k+1,j,i,icrm))*0.5-wwqui_bar_cen(k,icrm));
                wwqui_cloudy_cen(k,icrm)=wwqui_cloudy_cen(k,icrm)+mask_tmp * tempwork*tempwork + mask_tmp * host_state_shoc_tke(k,j,i,icrm)/3.0;
              } // if(icl==CLD)
            } // end of for (int k=0

            // ! Now, we can do a similar aggregation for the cell boundaries.
            for (int k = 0; k < nzstag; ++k) {
              mask_tmp = mask_bnd(k,j,i,icrm,iCL,QUI,iPR)/nxy;
              // !NOTE: technically we should interpolate and not do a simple
              //     average to get density at the cell interface
              km0 = std::min(nz,k);   // Liran: should we change k to k-1 in c++?
              km1 = std::max(0,k-1);
/*
              ! calculate the mean vertical velocity over the quiescent class  +++mhwang
              ! wwqui_bar_bnd is used in both all sky and cloudy sky.
              ! when wwqui_cloudy_bar_bnd was used for cloudy sky, wwqui_cloudy_bnd will be smaller than wwqui_bnd.
              !
*/
              tempwork = (crm_wvel(k,j,i,icrm)-wwqui_bar_cen(k,icrm));
              wwqui_bnd(k,icrm) = wwqui_bnd(k,icrm)+mask_tmp * tempwork*tempwork + mask_tmp * (host_state_shoc_tke(km0,j,i,icrm)+host_state_shoc_tke(km1,j,i,icrm)) * 0.5/3.0;
              if(iCL==CLD){
                tempwork = (crm_wvel(k,j,i,icrm)-wwqui_bar_cen(k,icrm));
                wwqui_cloudy_bnd(k,icrm)=wwqui_cloudy_bnd(k,icrm)+mask_tmp * tempwork*tempwork + mask_tmp * (host_state_shoc_tke(km0,j,i,icrm)+host_state_shoc_tke(km1,j,i,icrm))* 0.5/3.0;
              } // if(icl==CLD)
            } // for (int k = 1; k < nzstag; ++k)  
          } // end of for (int iPR=0
        } // end of for (int iCL=0
      } // end of for (int i=0
    } // end of for (int j=0

    // testing small queiscent fraction +++mhwang
    real temp_area_cen;
    for (int k=0; k<nz; k++) {
      for (int iCL=0; iCL<=NCLASS_CL; iCL++) {
        for (int iPR=0; iPR<=NCLASS_PR; iPR++) {
          temp_area_cen = temp_area_cen + area_cen_final(k,icrm,iCL,0,iPR);
        }
      }
      if (temp_area_cen < 1.0e-3){
        //printf("%s %.2f\n", "ecpp, area_cen_final, quiescent", temp_area_cen);
      }
    } // end of for (int k=0
// ==============================End of call categorization_stats============
  }); // end of parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) 
/*
      ! If we want final area categories based on the last avg1 period in each
      ! avg2 then we need to zero out the running sum just created for the areas
      ! if it is not the last block of time in ntavg2
*/    
  printf("Liran check categorization_stats 12\n");
  parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) {
    for (int k=0; k<nz; k++) {
      ecpp_sum_wwqui_bar_cen(k,icrm)        = ecpp_sum_wwqui_bar_cen(k,icrm)        + wwqui_bar_cen(k,icrm);
      ecpp_sum_wwqui_cloudy_bar_cen(k,icrm) = ecpp_sum_wwqui_cloudy_bar_cen(k,icrm) + wwqui_cloudy_bar_cen(k,icrm);
      ecpp_sum_wwqui_bar_bnd(k,icrm)        = ecpp_sum_wwqui_bar_bnd(k,icrm)        + wwqui_bar_bnd(k,icrm);
      ecpp_sum_wwqui_cloudy_bar_bnd(k,icrm) = ecpp_sum_wwqui_cloudy_bar_bnd(k,icrm) + wwqui_cloudy_bar_bnd(k,icrm);
      ecpp_sum_tbeg(k,icrm)                 =  ecpp_sum_tbeg(k,icrm)                + acldy_cen_tbeg(k,icrm);
    }
/*
    for (int k=0; k<nz; k++) {
      printf("\narea_bnd (iPR=0,iTR=0,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(0,0,0,k,icrm));
      printf("\narea_bnd (iPR=0,iTR=0,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(0,0,1,k,icrm));
      printf("\narea_bnd (iPR=0,iTR=1,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(0,1,0,k,icrm));
      printf("\narea_bnd (iPR=0,iTR=1,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(0,1,1,k,icrm));
      printf("\narea_bnd (iPR=0,iTR=2,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(0,2,0,k,icrm));
      printf("\narea_bnd (iPR=0,iTR=2,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(0,2,1,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=0,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(1,0,0,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=0,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(1,0,1,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=1,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(1,1,0,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=1,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(1,1,1,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=2,iCL=0): %d %d %.2f  : ", k,icrm,area_bnd(1,2,0,k,icrm));
      printf("\narea_bnd (iPR=1,iTR=2,iCL=1): %d %d %.2f  : ", k,icrm,area_bnd(1,2,1,k,icrm));
    }
*/
    for (int iCL=0; iCL<NCLASS_CL; iCL++) {
      for (int iTR=0; iTR<ndraft_max; iTR++) {
        for (int iPR=0; iPR<NCLASS_PR; iPR++) {
          for (int k=0; k<nz; k++) {
            ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm)        =  ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm)       + area_cen(iPR,iTR,iCL,k,icrm);
            ecpp_sum_area_cen_final(iPR,iTR,iCL,k,icrm)  =  ecpp_sum_area_cen_final(iPR,iTR,iCL,k,icrm) + area_cen_final(iPR,iTR,iCL,k,icrm);
            ecpp_sum_rh_cen(iPR,iTR,iCL,k,icrm)          =  ecpp_sum_rh_cen(iPR,iTR,iCL,k,icrm)         + rh_cen(iPR,iTR,iCL,k,icrm);
            ecpp_sum_qcloud_cen(iPR,iTR,iCL,k,icrm)      =  ecpp_sum_qcloud_cen(iPR,iTR,iCL,k,icrm)     + qcloud_cen(iPR,iTR,iCL,k,icrm);
            ecpp_sum_qice_cen(iPR,iTR,iCL,k,icrm)        =  ecpp_sum_qice_cen(iPR,iTR,iCL,k,icrm)       + qice_cen(iPR,iTR,iCL,k,icrm);
            ecpp_sum_precsolidcen(iPR,iTR,iCL,k,icrm)    =  ecpp_sum_precsolidcen(iPR,iTR,iCL,k,icrm)   + qrain_cen(iPR,iTR,iCL,k,icrm);
          }
        }
      }
    }

    for (int iCL=0; iCL<NCLASS_CL; iCL++) {
      for (int iTR=0; iTR<ndraft_max; iTR++) {
        for (int iPR=0; iPR<NCLASS_PR; iPR++) {
          for (int k=0; k<nzi; k++) {
            ecpp_sum_area_bnd(iPR,iTR,iCL,k,icrm)        =  ecpp_sum_area_bnd(iPR,iTR,iCL,k,icrm)       + area_bnd(iPR,iTR,iCL,k,icrm);
            ecpp_sum_area_bnd_final(iPR,iTR,iCL,k,icrm)  =  ecpp_sum_area_bnd_final(iPR,iTR,iCL,k,icrm) + area_bnd_final(iPR,iTR,iCL,k,icrm);
            ecpp_sum_mass_bnd(iPR,iTR,iCL,k,icrm)        =  ecpp_sum_mass_bnd(iPR,iTR,iCL,k,icrm)       + mass_bnd(iPR,iTR,iCL,k,icrm);
          }
        }
      }
    }
    crm_level2_cnt(icrm) = crm_level2_cnt(icrm) + 1;
  }); //end of parallel_for( SimpleBounds<1>(nens) , YAKL_LAMBDA (int icrm) 
  //printf("\nLevel 1: Value of crm_level2_cnt: %d: ", crm_level2_cnt(0));

  // Done with time level one averages so zero them out for next period.
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int iz, int iy, int ix, int iens) {
    liran_test4davg(iz,iy,ix,iens)    = 0;
    qcloudsum1(iz,iy,ix,iens)    = 0;
    qrainsum1(iz,iy,ix,iens)     = 0;
    qicesum1(iz,iy,ix,iens)      = 0;
    qsnowsum1(iz,iy,ix,iens)     = 0;
    altsum1(iz,iy,ix,iens)       = 0;
    ecppwwsum1(iz,iy,ix,iens)    = 0;
    rhsum1(iz,iy,ix,iens)        = 0;
  });



printf("\nLiran check start ECPP ecpp_crm_stat 04\n");
} // if (runcount >=ntavg1 && runcount % ntavg1 == 0)   --------------end of level 1 averaging-----------------------

// ============================== End of time level one averaging period ======

// Check if we have reached the end of a level 2 averaging period.
// Liran: Again, we haven't distinguish level 2 averaging from level 1? 

// ! Check if we have reached the end of a level 2 averaging period.
if (runcount >=ntavg2 && runcount % ntavg2 == 0){
    printf("Liran check categorization_stats 15\n");
/*
      ! Turn the running sums into averages. ncnt1 in this case is the number
      ! of calls to categorization_stats during the level 2 averaging period,
      ! which increment the bnd/cen arrays.
*/
  //printf("\nLiran check start level2 averaging 00\n");
  //printf("\nLevel2 Value of crm_level2_cnt: %d: ", crm_level2_cnt(0));
  
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int icrm) {
    ecpp_cat_wwqui_bar_cen(k_crm,icrm)         = ecpp_sum_wwqui_bar_cen(k_crm,icrm)     /crm_level2_cnt(icrm);
    ecpp_cat_wwqui_cloudy_bar_cen(k_crm,icrm)  = ecpp_sum_wwqui_cloudy_bar_cen(k_crm,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_tbeg(k_crm,icrm)                  = ecpp_sum_tbeg(k_crm,icrm) /crm_level2_cnt(icrm);
  });
  //printf("\nLiran check start level2 averaging 01\n");
  parallel_for( SimpleBounds<2>(nzi,nens) , YAKL_LAMBDA (int k_crm, int icrm) {
    ecpp_cat_wwqui_bar_bnd(k_crm,icrm)             = ecpp_sum_wwqui_bar_bnd(k_crm,icrm)       /crm_level2_cnt(icrm);
    ecpp_cat_wwqui_cloudy_bar_bnd(k_crm,icrm)      = ecpp_sum_wwqui_cloudy_bar_bnd(k_crm,icrm)/crm_level2_cnt(icrm);
  });
  printf("Liran check categorization_stats 16\n");
  //printf("\nLiran check start level2 averaging 02\n");
  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm)         = ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_area_cen_final(iPR,iTR,iCL,k,icrm)   = ecpp_sum_area_cen_final(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_rh_cen(iPR,iTR,iCL,k,icrm)           = ecpp_sum_rh_cen(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_qcloud_cen(iPR,iTR,iCL,k,icrm)       = ecpp_sum_qcloud_cen(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_qice_cen(iPR,iTR,iCL,k,icrm)         = ecpp_sum_qice_cen(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_precsolidcen(iPR,iTR,iCL,k,icrm)     = ecpp_sum_precsolidcen(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
  });
  printf("\nLiran check start level2 averaging 03\n");
  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm)         = ecpp_sum_area_bnd(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_area_bnd_final(iPR,iTR,iCL,k,icrm)   = ecpp_sum_area_bnd_final(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
    ecpp_cat_mass_bnd(iPR,iTR,iCL,k,icrm)         = ecpp_sum_mass_bnd(iPR,iTR,iCL,k,icrm)/crm_level2_cnt(icrm);
  });
printf("\nLiran check start level2 averaging %d \n");
/*
! get in-cloud value for rh, qcloud, qrain, qice, qsnow, qgraup,
! percr, precsolid, and precall. (qlsink is already in-cloud values)
*/
  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    if (ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm) >afrac_cut){
      //ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm)            = ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm)/ecpp_sum_area_cen(iPR,iTR,iCL,k,icrm) ;
      //ecpp_cat_area_cen_final(iPR,iTR,iCL,k,icrm)      = ecpp_cat_area_cen_final(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      ecpp_cat_rh_cen(iPR,iTR,iCL,k,icrm)              = ecpp_cat_rh_cen(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      ecpp_cat_qcloud_cen(iPR,iTR,iCL,k,icrm)          = ecpp_cat_qcloud_cen(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      ecpp_cat_qice_cen(iPR,iTR,iCL,k,icrm)            = ecpp_cat_qice_cen(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      ecpp_cat_precsolidcen(iPR,iTR,iCL,k,icrm)        = ecpp_cat_precsolidcen(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
    } else {
      //ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm)            = 0.0;
      //ecpp_cat_area_cen_final(iPR,iTR,iCL,k,icrm)      = 0.0;
      ecpp_cat_rh_cen(iPR,iTR,iCL,k,icrm)              = 0.0;
      ecpp_cat_qcloud_cen(iPR,iTR,iCL,k,icrm)          = 0.0;
      ecpp_cat_qice_cen(iPR,iTR,iCL,k,icrm)            = 0.0;
      ecpp_cat_precsolidcen(iPR,iTR,iCL,k,icrm)        = 0.0;    }  
  });
  printf("\nLiran check start level2 averaging 05\n");
  parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    if (ecpp_sum_area_bnd(iPR,iTR,iCL,k,icrm) >afrac_cut){
      //ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm)            = ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      //ecpp_cat_area_bnd_final(iPR,iTR,iCL,k,icrm)      = ecpp_cat_area_bnd_final(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
      ecpp_cat_mass_bnd(iPR,iTR,iCL,k,icrm)            = ecpp_cat_mass_bnd(iPR,iTR,iCL,k,icrm)/ecpp_cat_area_cen(iPR,iTR,iCL,k,icrm) ;
    } else {
      //ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm)            = 0.0;
      //ecpp_cat_area_bnd_final(iPR,iTR,iCL,k,icrm)      = 0.0;
      ecpp_cat_mass_bnd(iPR,iTR,iCL,k,icrm)            = 0.0;
    }  
  });


  //parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,nzi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
  //   ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm) = iPR*100000.0+iTR*10000.0+iCL*1000.0+k*10.0+icrm;
  //   printf("\nnecpp_cat_area_bnd: %d %d %d %d %d %.2f  : ", iPR,iTR,iCL,k,icrm,ecpp_cat_area_bnd(iPR,iTR,iCL,k,icrm));
  //});

  printf("\nLiran check end level2 averaging\n");

} // end of if (runcount >=ntavg2 && runcount % ntavg2 == 0)

}

inline void pam_ecpp_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  //printf("\nLiran check pam_ecpp_copy_to_host 0\n");
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  auto nens            = coupler.get_option<int>("ncrms");    // Note that nz   = crm_nz
  auto nz              = coupler.get_option<int>("crm_nz");
  auto NCLASS_CL     = coupler.get_option<int>("ecpp_NCLASS_CL");
  auto ndraft_max    = coupler.get_option<int>("ecpp_ndraft_max");
  auto NCLASS_PR     = coupler.get_option<int>("ecpp_NCLASS_PR");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");
  auto gcm_nlevi    = coupler.get_option<int>("gcm_nlevi");
  dm_device.register_and_allocate<real>("ecpp_output_acen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_acen_tf_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_abnd_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlevi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_nlevi","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_abnd_tf_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlevi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_nlevi","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_massflxbnd_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlevi,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_nlevi","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_rhcen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_qcloudcen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_qlsinkcen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_precrcen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_precsolidcen_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ecpp_output_tbeg_temp", "<description>", {NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens}, {"NCLASS_PR","ndraft_max","NCLASS_CL","gcm_lev","nens"});

  auto ecpp_output_acen_temp         = dm_device.get<real,5>("ecpp_output_acen_temp");
  auto ecpp_output_acen_tf_temp      = dm_device.get<real,5>("ecpp_output_acen_tf_temp");
  auto ecpp_output_abnd_temp         = dm_device.get<real,5>("ecpp_output_abnd_temp");
  auto ecpp_output_abnd_tf_temp      = dm_device.get<real,5>("ecpp_output_abnd_tf_temp");
  auto ecpp_output_massflxbnd_temp   = dm_device.get<real,5>("ecpp_output_massflxbnd_temp");
  auto ecpp_output_rhcen_temp        = dm_device.get<real,5>("ecpp_output_rhcen_temp");
  auto ecpp_output_qcloudcen_temp    = dm_device.get<real,5>("ecpp_output_qcloudcen_temp");
  auto ecpp_output_qlsinkcen_temp    = dm_device.get<real,5>("ecpp_output_qlsinkcen_temp");
  auto ecpp_output_precrcen_temp     = dm_device.get<real,5>("ecpp_output_precrcen_temp");
  auto ecpp_output_precsolidcen_temp = dm_device.get<real,5>("ecpp_output_precsolidcen_temp");
  auto ecpp_output_tbeg_temp         = dm_device.get<real,5>("ecpp_output_tbeg_temp");

  auto ecpp_output_wwqui_cen         = dm_host.get<real,2>("ecpp_output_wwqui_cen");
  auto ecpp_output_wwqui_cloudy_cen  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_cen");
  auto ecpp_output_wwqui_bnd         = dm_host.get<real,2>("ecpp_output_wwqui_bnd");
  auto ecpp_output_wwqui_cloudy_bnd  = dm_host.get<real,2>("ecpp_output_wwqui_cloudy_bnd");
  auto ecpp_output_acen              = dm_host.get<real,5>("ecpp_output_acen");
  auto ecpp_output_abnd              = dm_host.get<real,5>("ecpp_output_abnd");
  auto ecpp_output_acen_tf           = dm_host.get<real,5>("ecpp_output_acen_tf");
  auto ecpp_output_abnd_tf           = dm_host.get<real,5>("ecpp_output_abnd_tf");
  auto ecpp_output_massflxbnd        = dm_host.get<real,5>("ecpp_output_massflxbnd");
  auto ecpp_output_rhcen             = dm_host.get<real,5>("ecpp_output_rhcen");
  auto ecpp_output_qcloudcen         = dm_host.get<real,5>("ecpp_output_qcloudcen");
  auto ecpp_output_qlsinkcen         = dm_host.get<real,5>("ecpp_output_qlsinkcen");
  auto ecpp_output_precrcen          = dm_host.get<real,5>("ecpp_output_precrcen");
  auto ecpp_output_precsolidcen      = dm_host.get<real,5>("ecpp_output_precsolidcen");
  auto ecpp_output_tbeg              = dm_host.get<real,2>("ecpp_output_tbeg");
  auto ecpp_cat_wwqui_bar_cen         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_cen");
  auto ecpp_cat_wwqui_bar_bnd         = dm_device.get<real,2>("ecpp_cat_wwqui_bar_bnd");
  auto ecpp_cat_wwqui_cloudy_bar_cen  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_cen");
  auto ecpp_cat_wwqui_cloudy_bar_bnd  = dm_device.get<real,2>("ecpp_cat_wwqui_cloudy_bar_bnd");
  auto ecpp_cat_area_cen_final        = dm_device.get<real,5>("ecpp_cat_area_cen_final");
  auto ecpp_cat_area_cen              = dm_device.get<real,5>("ecpp_cat_area_cen");
  auto ecpp_cat_area_bnd_final        = dm_device.get<real,5>("ecpp_cat_area_bnd_final");
  auto ecpp_cat_area_bnd          = dm_device.get<real,5>("ecpp_cat_area_bnd");
  auto ecpp_cat_mass_bnd          = dm_device.get<real,5>("ecpp_cat_mass_bnd");
  auto ecpp_cat_rh_cen            = dm_device.get<real,5>("ecpp_cat_rh_cen");
  auto ecpp_cat_qcloud_cen        = dm_device.get<real,5>("ecpp_cat_qcloud_cen");
  auto ecpp_cat_qice_cen          = dm_device.get<real,5>("ecpp_cat_qice_cen");
  auto ecpp_cat_precsolidcen      = dm_device.get<real,5>("ecpp_cat_precsolidcen");
  auto ecpp_cat_tbeg              = dm_device.get<real,2>("ecpp_cat_tbeg");
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
    printf("\nLiran check pam_ecpp_copy_to_host 1\n");

    parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
      ecpp_output_acen_temp(iPR,iTR,iCL,k,icrm)         = 0;
      ecpp_output_acen_tf_temp(iPR,iTR,iCL,k,icrm)      = 0;
      ecpp_output_rhcen_temp(iPR,iTR,iCL,k,icrm)        = 0;
      ecpp_output_qcloudcen_temp(iPR,iTR,iCL,k,icrm)    = 0;
      ecpp_output_qlsinkcen_temp(iPR,iTR,iCL,k,icrm)    = 0;
      ecpp_output_precrcen_temp(iPR,iTR,iCL,k,icrm)     = 0;
      ecpp_output_precsolidcen_temp(iPR,iTR,iCL,k,icrm) = 0;
    });
    parallel_for(SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlevi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
      ecpp_output_abnd_temp(iPR,iTR,iCL,k,icrm)         = 0;
      ecpp_output_abnd_tf_temp(iPR,iTR,iCL,k,icrm)      = 0;
      ecpp_output_massflxbnd_temp(iPR,iTR,iCL,k,icrm)   = 0;
    });
    printf("\nLiran check pam_ecpp_copy_to_host 2\n");
    parallel_for("set output forcing tendencies", SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlev,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k_gcm,int icrm) {
        int k_crm = gcm_nlev-k_gcm;
        int CLR = 0;       // Clear sub-class
        int PRN = 0;       // Not precipitating sub-class
        int QUI       = 1; // Quiescent class
        if (k_crm<nz) {
          ecpp_output_acen_temp (iPR,iTR,iCL,k_gcm,icrm)         = ecpp_cat_area_cen (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_acen_tf_temp (iPR,iTR,iCL,k_gcm,icrm)      = ecpp_cat_area_cen_final (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_rhcen_temp (iPR,iTR,iCL,k_gcm,icrm)        = ecpp_cat_rh_cen (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_qcloudcen_temp (iPR,iTR,iCL,k_gcm,icrm)    = ecpp_cat_qcloud_cen (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_qlsinkcen_temp (iPR,iTR,iCL,k_gcm,icrm)    = ecpp_cat_qice_cen (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_precrcen_temp (iPR,iTR,iCL,k_gcm,icrm)     = ecpp_cat_precsolidcen (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_precsolidcen_temp (iPR,iTR,iCL,k_gcm,icrm) = ecpp_cat_precsolidcen (iPR,iTR,iCL,k_crm,icrm);
          
        } else {
          // Liran: above crm levels are assumed to be clear, non-precipitation, and clear condition
          ecpp_output_acen_temp (PRN,QUI,CLR,k_gcm,icrm)         = 1.;
          ecpp_output_acen_tf_temp (PRN,QUI,CLR,k_gcm,icrm)      = 1.;
          ecpp_output_rhcen_temp (PRN,QUI,CLR,k_gcm,icrm)        = 1.;
          ecpp_output_qcloudcen_temp (PRN,QUI,CLR,k_gcm,icrm)    = 1.; 
          ecpp_output_qlsinkcen_temp (PRN,QUI,CLR,k_gcm,icrm)    = 1.;
          ecpp_output_precrcen_temp (PRN,QUI,CLR,k_gcm,icrm)     = 1.;
          ecpp_output_precsolidcen_temp (PRN,QUI,CLR,k_gcm,icrm) = 1.;
        }
    });
    printf("\nLiran check pam_ecpp_copy_to_host 3\n");
    parallel_for("set output forcing tendencies", SimpleBounds<5>(NCLASS_PR,ndraft_max,NCLASS_CL,gcm_nlevi,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k_gcm,int icrm) {
        int k_crm = gcm_nlev-1-k_gcm;
        int CLR = 0;       // Clear sub-class
        int PRN = 0;       // Not precipitating sub-class
        int QUI       = 1; // Quiescent class
        if (k_crm<nz+1) {
          ecpp_output_abnd_temp (iPR,iTR,iCL,k_gcm,icrm)         = ecpp_cat_area_bnd (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_abnd_tf_temp (iPR,iTR,iCL,k_gcm,icrm)      = ecpp_cat_area_bnd_final (iPR,iTR,iCL,k_crm,icrm);
          ecpp_output_massflxbnd_temp (iPR,iTR,iCL,k_gcm,icrm)   = ecpp_cat_mass_bnd (iPR,iTR,iCL,k_crm,icrm);
        } else {
          // Liran: above crm levels are assumed to be clear, non-precipitation, and clear condition
          ecpp_output_abnd_temp (PRN,QUI,CLR,k_gcm,icrm)         = 1.;
          ecpp_output_abnd_tf_temp (PRN,QUI,CLR,k_gcm,icrm)      = 1.;
          ecpp_output_massflxbnd_temp (PRN,QUI,CLR,k_gcm,icrm)   = 1.;
        }
    });
  printf("\nLiran check pam_ecpp_copy_to_host 4\n");
  // Copy the CRM ECPP state to host arrays
  ecpp_cat_wwqui_bar_cen                 .deep_copy_to( ecpp_output_wwqui_cen               );
  ecpp_cat_wwqui_cloudy_bar_cen          .deep_copy_to( ecpp_output_wwqui_cloudy_cen        );
  ecpp_cat_wwqui_bar_bnd                 .deep_copy_to( ecpp_output_wwqui_bnd               );
  ecpp_cat_wwqui_cloudy_bar_bnd          .deep_copy_to( ecpp_output_wwqui_cloudy_bnd        );
  ecpp_output_acen_temp                  .deep_copy_to( ecpp_output_acen                    );
  ecpp_output_abnd_temp                  .deep_copy_to( ecpp_output_abnd                    );
  ecpp_output_acen_tf_temp               .deep_copy_to( ecpp_output_acen_tf                 );
  ecpp_output_abnd_tf_temp               .deep_copy_to( ecpp_output_abnd_tf                 );
  ecpp_output_massflxbnd_temp            .deep_copy_to( ecpp_output_massflxbnd              );
  ecpp_output_rhcen_temp                 .deep_copy_to( ecpp_output_rhcen                   );
  ecpp_output_qcloudcen_temp             .deep_copy_to( ecpp_output_qcloudcen               );
  ecpp_output_qlsinkcen_temp             .deep_copy_to( ecpp_output_qlsinkcen               );
  ecpp_output_precrcen_temp              .deep_copy_to( ecpp_output_precrcen                );
  ecpp_output_precsolidcen_temp          .deep_copy_to( ecpp_output_precsolidcen            );
  ecpp_cat_tbeg                          .deep_copy_to( ecpp_output_tbeg                    );
  //printf("\nLiran check pam_ecpp_copy_to_host 5\n");
  yakl::fence();
  //------------------------------------------------------------------------------------------------
}
