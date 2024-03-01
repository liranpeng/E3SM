#pragma once

#include "pam_coupler.h"

inline void pam_ecpp_vars_register( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  //------------------------------------------------------------------------------------------------
  nclass_trx    = coupler.get_option<int>("ecpp_nclass_trx");
  nclass_cld    = coupler.get_option<int>("ecpp_nclass_cld");
  nclass_prc    = coupler.get_option<int>("ecpp_nclass_prc");
  //------------------------------------------------------------------------------------------------
  dm_device.register_and_allocate<int> ("ecpp_L1_cnt",   "# of level 2 avg count",     {nens}, {"nens"});
  dm_device.register_and_allocate<int> ("ecpp_L2_cnt",   "# of level 2 avg count",     {nens}, {"nens"});
  dm_device.register_and_allocate<int> ("ndown",         "# of down count",            {nens}, {"nens"});
  dm_device.register_and_allocate<int> ("nup",           "# of nup count",             {nens}, {"nens"});
  dm_device.register_and_allocate<int> ("kup_top",       "maximum kup",                {nens}, {"nens"});
  dm_device.register_and_allocate<int> ("kdown_top",     "maximum kdown",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wdown_bar",     "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wup_bar",       "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wdown_stddev",  "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wup_stddev",    "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wup_rms",       "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("wdown_rms",     "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("updraftbase",   "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("updrafttop",    "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("dndrafttop",    "<description>",              {nens}, {"nens"});
  dm_device.register_and_allocate<real>("dndraftbase",   "<description>",              {nens}, {"nens"});
  //------------------------------------------------------------------------------------------------
  // Level 1 running sums
  dm_device.register_and_allocate<real>("ecpp_L1_sum_cld","<description>", {nz,        nens}, {"z",          "nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_rho","<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_qc", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_qr", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_qi", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_qs", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_rh", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_sum_ww", "<description>", {nz+1,ny,nx,nens}, {"zp1","y","x","nens"});
  // Level 1 averages
  dm_device.register_and_allocate<real>("ecpp_L1_avg_cld","<description>", {nz,        nens}, {"z",          "nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_rho","<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_qc", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_qr", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_qi", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_qs", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_rh", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecpp_L1_avg_ww", "<description>", {nz+1,ny,nx,nens}, {"zp1","y","x","nens"});
  //------------------------------------------------------------------------------------------------
  // Level 2 running sums
  dm_device.register_and_allocate<real>("ecpp_L2_sum_tk", "<description>", {nz,  nens}, {"z",  "nens"}); // eddy viscosity m2/s
  //------------------------------------------------------------------------------------------------
  // 4D ECPP quantities
  dm_device.register_and_allocate<real>("qlsink_bf" ,    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("prain"     ,    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qcloud_bf" ,    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qcloud_bfsum1", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qgraupsum1",    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qlsinksum1",    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("precrsum1",     "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("precsolidsum1", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("precallsum1",   "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("cf3dsum1",      "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("tkesgssum1",    "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qlsink_bfsum1", "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("qvssum1",       "<description>", {nz,  ny,nx,nens}, {"z",  "y","x","nens"});
  dm_device.register_and_allocate<real>("ecppwwsqsum1",  "<description>", {nz+1,ny,nx,nens}, {"zp1","y","x","nens"});
  //------------------------------------------------------------------------------------------------
  // 2D ECPP quantities
  dm_device.register_and_allocate<real>("cldtot2d",         "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wwqui_all_cen",    "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wwqui_all_bnd",    "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("wwqui_cld_cen",    "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wwqui_cld_bnd",    "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<int> ("nup_k",            "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wup_bar_k",        "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<int> ("ndown_k",          "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wdown_bar_k",      "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wdown_stddev_k",   "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wup_stddev_k",     "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wup_rms_k",        "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wdown_rms_k",      "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("wup_rms_ksmo",     "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("wdown_rms_ksmo",   "<description>", {nz+1,nens}, {"zp1","nens"});
  //------------------------------------------------------------------------------------------------
  // variables categorized by the various transport and cloud categories
  dm_device.register_and_allocate<real>("ecpp_cat_wwqui_bar_cen",     "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_wwqui_cld_bar_cen", "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_tbeg",              "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_wwqui_bar_bnd",     "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_wwqui_cld_bar_bnd", "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_area_cen_final",    "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_area_cen",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_rh_cen",            "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_qcloud_cen",        "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_qice_cen",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_precsolidcen",      "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_area_bnd_final",    "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_area_bnd",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_cat_mass_bnd",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  //------------------------------------------------------------------------------------------------
  // aggregated variables used for categorization
  dm_device.register_and_allocate<real>("ecpp_sum_wwqui_bar_cen",     "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_wwqui_cld_bar_cen", "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_tbeg",              "<description>", {nz,  nens}, {"z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_wwqui_bar_bnd",     "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_wwqui_cld_bar_bnd", "<description>", {nz+1,nens}, {"zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_area_cen_final",    "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_area_cen",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_rh_cen",            "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_qcloud_cen",        "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_qice_cen",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_precsolidcen",      "<description>", {nclass_prc,nclass_trx,nclass_cld,nz,  nens}, {"nclass_prc","nclass_trx","nclass_cld","z",  "nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_area_bnd_final",    "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_area_bnd",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  dm_device.register_and_allocate<real>("ecpp_sum_mass_bnd",          "<description>", {nclass_prc,nclass_trx,nclass_cld,nz+1,nens}, {"nclass_prc","nclass_trx","nclass_cld","zp1","nens"});
  //------------------------------------------------------------------------------------------------
}

// initialize values after they are registered and allocated
inline void pam_ecpp_vars_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  //------------------------------------------------------------------------------------------------
  nclass_trx    = coupler.get_option<int>("ecpp_nclass_trx");
  nclass_cld    = coupler.get_option<int>("ecpp_nclass_cld");
  nclass_prc    = coupler.get_option<int>("ecpp_nclass_prc");
  //------------------------------------------------------------------------------------------------
  // initialize 1D quantities
  auto L1_cnt          = dm_device.get<int, 1>("ecpp_L1_cnt");
  auto L2_cnt          = dm_device.get<int, 1>("ecpp_L2_cnt");
  auto kup_top         = dm_device.get<int, 1>("kup_top");
  auto kdown_top       = dm_device.get<int, 1>("kdown_top");
  auto ndown           = dm_device.get<int, 1>("ndown");
  auto nup             = dm_device.get<int, 1>("nup");
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
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    L1_cnt        (icrm) = 0;
    L2_cnt        (icrm) = 0;
    kup_top       (icrm) = 0;
    kdown_top     (icrm) = 0;
    ndown         (icrm) = 0;
    nup           (icrm) = 0;
    updraftbase   (icrm) = 0;
    updrafttop    (icrm) = nz - 1;
    dndrafttop    (icrm) = nz - 1;
    dndraftbase   (icrm) = 0;
    wup_bar       (icrm) = 0.0;
    wdown_bar     (icrm) = 0.0;
    wup_stddev    (icrm) = 0.0;
    wdown_stddev  (icrm) = 0.0;
    wup_rms       (icrm) = 0.0;
    wdown_rms     (icrm) = 0.0;
  });
  //------------------------------------------------------------------------------------------------
  // initialize 2D quantities
  auto xkhvsum                    = dm_device.get<real,2>("xkhvsum");
  auto nup_k                      = dm_device.get<int, 2>("nup_k");
  auto ndown_k                    = dm_device.get<int, 2>("ndown_k");
  auto wup_bar_k                  = dm_device.get<real,2>("wup_bar_k");
  auto wdown_bar_k                = dm_device.get<real,2>("wdown_bar_k");
  auto wup_stddev_k               = dm_device.get<real,2>("wup_stddev_k");
  auto wdown_stddev_k             = dm_device.get<real,2>("wdown_stddev_k");
  auto wup_rms_k                  = dm_device.get<real,2>("wup_rms_k");
  auto wdown_rms_k                = dm_device.get<real,2>("wdown_rms_k");
  auto wup_rms_ksmo               = dm_device.get<real,2>("wup_rms_ksmo");
  auto wdown_rms_ksmo             = dm_device.get<real,2>("wdown_rms_ksmo");
  auto wwqui_all_cen              = dm_device.get<real,2>("wwqui_all_cen");
  auto wwqui_all_bnd              = dm_device.get<real,2>("wwqui_all_bnd");
  auto wwqui_cld_cen              = dm_device.get<real,2>("wwqui_cld_cen");
  auto wwqui_cld_bnd              = dm_device.get<real,2>("wwqui_cld_bnd");
  auto cldtot2d                   = dm_device.get<real,2>("cldtot2d");
  auto ecpp_sum_wwqui_bar_cen     = dm_device.get<real,2>("ecpp_sum_wwqui_bar_cen");
  auto ecpp_sum_wwqui_cld_bar_cen = dm_device.get<real,2>("ecpp_sum_wwqui_cld_bar_cen");
  auto ecpp_sum_wwqui_bar_bnd     = dm_device.get<real,2>("ecpp_sum_wwqui_bar_bnd");
  auto ecpp_sum_wwqui_cld_bar_bnd = dm_device.get<real,2>("ecpp_sum_wwqui_cld_bar_bnd");
  auto ecpp_sum_tbeg              = dm_device.get<real,2>("ecpp_sum_tbeg");
  auto ecpp_cat_wwqui_bar_cen     = dm_device.get<real,2>("ecpp_cat_wwqui_bar_cen");
  auto ecpp_cat_wwqui_bar_bnd     = dm_device.get<real,2>("ecpp_cat_wwqui_bar_bnd");
  auto ecpp_cat_wwqui_cld_bar_cen = dm_device.get<real,2>("ecpp_cat_wwqui_cld_bar_cen");
  auto ecpp_cat_wwqui_cld_bar_bnd = dm_device.get<real,2>("ecpp_cat_wwqui_cld_bar_bnd");
  auto ecpp_cat_tbeg              = dm_device.get<real,2>("ecpp_cat_tbeg");
  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int iz, int iens) {
    xkhvsum                        (iz,iens) = 0;
    nup_k                          (iz,iens) = 0;
    wup_bar_k                      (iz,iens) = 0;
    ndown_k                        (iz,iens) = 0;
    wdown_bar_k                    (iz,iens) = 0;
    wdown_stddev_k                 (iz,iens) = 0;
    wup_stddev_k                   (iz,iens) = 0;
    wup_rms_k                      (iz,iens) = 0;
    wdown_rms_k                    (iz,iens) = 0;
    wwqui_all_cen                  (iz,iens) = 0;
    wwqui_cld_cen                  (iz,iens) = 0;
    cldtot2d                       (iz,iens) = 0;
    ecpp_cat_wwqui_bar_cen         (iz,iens) = 0;
    ecpp_cat_wwqui_cld_bar_cen     (iz,iens) = 0;
    ecpp_cat_tbeg                  (iz,iens) = 0;
    ecpp_sum_wwqui_bar_cen         (iz,iens) = 0;
    ecpp_sum_wwqui_cld_bar_cen     (iz,iens) = 0;
    ecpp_sum_tbeg                  (iz,iens) = 0;
  });
  parallel_for(SimpleBounds<2>(nz+1,nens), YAKL_LAMBDA (int iz, int iens) {
    wwqui_all_bnd                  (iz,iens)= 0; 
    wwqui_cld_bnd                  (iz,iens)= 0;
    wup_rms_ksmo                   (iz,iens)= 0;
    wdown_rms_ksmo                 (iz,iens)= 0;
    ecpp_cat_wwqui_bar_bnd         (iz,iens)= 0;
    ecpp_cat_wwqui_cld_bar_bnd     (iz,iens)= 0;
    ecpp_sum_wwqui_bar_bnd         (iz,iens)= 0;
    ecpp_sum_wwqui_cld_bar_bnd     (iz,iens)= 0;
  });
  //------------------------------------------------------------------------------------------------
  // initialize 4D quantities
  auto qlsink_bf            = dm_device.get<real,4>("qlsink_bf");
  auto prain                = dm_device.get<real,4>("prain");
  auto qcloud_bf            = dm_device.get<real,4>("qcloud_bf");
  auto qcloud_bfsum1        = dm_device.get<real,4>("qcloud_bfsum1");
  auto qgraupsum1           = dm_device.get<real,4>("qgraupsum1");
  auto qlsinksum1           = dm_device.get<real,4>("qlsinksum1");
  auto precrsum1            = dm_device.get<real,4>("precrsum1");
  auto precsolidsum1        = dm_device.get<real,4>("precsolidsum1");
  auto precallsum1          = dm_device.get<real,4>("precallsum1");
  auto cf3dsum1             = dm_device.get<real,4>("cf3dsum1");
  auto ecppwwsqsum1         = dm_device.get<real,4>("ecppwwsqsum1");
  auto tkesgssum1           = dm_device.get<real,4>("tkesgssum1");
  auto qlsink_bfsum1        = dm_device.get<real,4>("qlsink_bfsum1");
  auto qvssum1              = dm_device.get<real,4>("qvssum1");
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int iz, int iy, int ix, int iens) {
    qlsink_bf     (iz,iy,ix,iens) = 0;
    prain         (iz,iy,ix,iens) = 0;
    qcloud_bf     (iz,iy,ix,iens) = 0;
    qcloud_bfsum1 (iz,iy,ix,iens) = 0;
    qgraupsum1    (iz,iy,ix,iens) = 0;
    qlsinksum1    (iz,iy,ix,iens) = 0;
    precrsum1     (iz,iy,ix,iens) = 0;
    precsolidsum1 (iz,iy,ix,iens) = 0;
    precallsum1   (iz,iy,ix,iens) = 0;
    cf3dsum1      (iz,iy,ix,iens) = 0;
    ecppwwsqsum1  (iz,iy,ix,iens) = 0;
    tkesgssum1    (iz,iy,ix,iens) = 0;
    qlsink_bfsum1 (iz,iy,ix,iens) = 0;
    qvssum1       (iz,iy,ix,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // initialize 5D quantities
  auto ecpp_sum_area_cen_final   = dm_device.get<real,5>("ecpp_sum_area_cen_final");
  auto ecpp_sum_area_cen         = dm_device.get<real,5>("ecpp_sum_area_cen");
  auto ecpp_sum_area_bnd_final   = dm_device.get<real,5>("ecpp_sum_area_bnd_final");
  auto ecpp_sum_area_bnd         = dm_device.get<real,5>("ecpp_sum_area_bnd");
  auto ecpp_sum_mass_bnd         = dm_device.get<real,5>("ecpp_sum_mass_bnd");
  auto ecpp_sum_rh_cen           = dm_device.get<real,5>("ecpp_sum_rh_cen");
  auto ecpp_sum_qcloud_cen       = dm_device.get<real,5>("ecpp_sum_qcloud_cen");
  auto ecpp_sum_qice_cen         = dm_device.get<real,5>("ecpp_sum_qice_cen");
  auto ecpp_sum_precsolidcen     = dm_device.get<real,5>("ecpp_sum_precsolidcen");
  auto ecpp_cat_area_cen_final   = dm_device.get<real,5>("ecpp_cat_area_cen_final");
  auto ecpp_cat_area_cen         = dm_device.get<real,5>("ecpp_cat_area_cen");
  auto ecpp_cat_area_bnd_final   = dm_device.get<real,5>("ecpp_cat_area_bnd_final");
  auto ecpp_cat_area_bnd         = dm_device.get<real,5>("ecpp_cat_area_bnd");
  auto ecpp_cat_mass_bnd         = dm_device.get<real,5>("ecpp_cat_mass_bnd");
  auto ecpp_cat_rh_cen           = dm_device.get<real,5>("ecpp_cat_rh_cen");
  auto ecpp_cat_qcloud_cen       = dm_device.get<real,5>("ecpp_cat_qcloud_cen");
  auto ecpp_cat_qice_cen         = dm_device.get<real,5>("ecpp_cat_qice_cen");
  auto ecpp_cat_precsolidcen     = dm_device.get<real,5>("ecpp_cat_precsolidcen");
  parallel_for(SimpleBounds<5>(nclass_prc,nclass_trx,nclass_cld,nz,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    ecpp_cat_area_cen_final   (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_area_cen         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_rh_cen           (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_qcloud_cen       (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_qice_cen         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_precsolidcen     (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_area_cen_final   (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_area_cen         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_rh_cen           (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_qcloud_cen       (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_qice_cen         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_precsolidcen     (iPR,iTR,iCL,k,icrm) = 0;
  });
  parallel_for(SimpleBounds<5>(nclass_prc,nclass_trx,nclass_cld,nz+1,nens), YAKL_LAMBDA (int iPR,int iTR,int iCL,int k,int icrm) {
    ecpp_cat_area_bnd_final   (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_area_bnd         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_cat_mass_bnd         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_area_bnd_final   (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_area_bnd         (iPR,iTR,iCL,k,icrm) = 0;
    ecpp_sum_mass_bnd         (iPR,iTR,iCL,k,icrm) = 0;
  });
  //------------------------------------------------------------------------------------------------
}