#pragma once

#include "pam_coupler.h"


// zero out L1 running sums
inline void pam_ecpp_L1_zero_sums( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  auto L1_sum_cld      = dm_device.get<real,2>("ecpp_L1_sum_cld");
  auto L1_sum_rho      = dm_device.get<real,4>("ecpp_L1_sum_rho");
  auto L1_sum_qc       = dm_device.get<real,4>("ecpp_L1_sum_qc");
  auto L1_sum_qr       = dm_device.get<real,4>("ecpp_L1_sum_qr");
  auto L1_sum_qi       = dm_device.get<real,4>("ecpp_L1_sum_qi");
  auto L1_sum_qs       = dm_device.get<real,4>("ecpp_L1_sum_qs");
  auto L1_sum_ww       = dm_device.get<real,4>("ecpp_L1_sum_ww");
  auto L1_sum_rh       = dm_device.get<real,4>("ecpp_L1_sum_rh");
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<2>(nz, nens), YAKL_LAMBDA (int k, int iens) {
    L1_sum_cld(k,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    L1_sum_rho(k,j,i,iens) = 0;
    L1_sum_qc (k,j,i,iens) = 0;
    L1_sum_qr (k,j,i,iens) = 0;
    L1_sum_qi (k,j,i,iens) = 0;
    L1_sum_qs (k,j,i,iens) = 0;
    L1_sum_ww (k,j,i,iens) = 0;
    L1_sum_rh (k,j,i,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
}


// update L1 running sums
inline void pam_ecpp_L1_update_sums( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  auto ref_pres        = dm_device.get<real const,2>("ref_pres"   );
  auto temp            = dm_device.get<real const,4>("temp"       );
  auto rho_d           = dm_device.get<real const,4>("density_dry");
  auto rho_v           = dm_device.get<real const,4>("water_vapor");
  auto rho_l           = dm_device.get<real const,4>("cloud_water");
  auto rho_i           = dm_device.get<real const,4>("ice"        );
  auto rho_r           = dm_device.get<real const,4>("rain"       );
  auto wvel            = dm_device.get<real const,4>("wvel"       );
  auto cldfrac         = dm_device.get<real const,4>("cldfrac"    );
  //------------------------------------------------------------------------------------------------
  // real2d acldy_cen_tbeg      ("acldy_cen_tbeg",      nz,  nens);
  auto L1_cnt          = dm_device.get<int, 1>("ecpp_L1_cnt");
  auto L1_sum_cld      = dm_device.get<real,2>("ecpp_L1_sum_cld");
  auto L1_sum_rho      = dm_device.get<real,4>("ecpp_L1_sum_rho");
  auto L1_sum_qc       = dm_device.get<real,4>("ecpp_L1_sum_qc");
  auto L1_sum_qr       = dm_device.get<real,4>("ecpp_L1_sum_qr");
  auto L1_sum_qi       = dm_device.get<real,4>("ecpp_L1_sum_qi");
  auto L1_sum_qs       = dm_device.get<real,4>("ecpp_L1_sum_qs");
  auto L1_sum_ww       = dm_device.get<real,4>("ecpp_L1_sum_ww");
  auto L1_sum_rh       = dm_device.get<real,4>("ecpp_L1_sum_rh");
  //------------------------------------------------------------------------------------------------
  // Increment the running sums for the averaging period
  real r_nx_ny  = 1.0/(nx*ny);
  parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
      // calculate saturation water vapor pressure (Pa)
      real evs = esatw_crm(temp(k,j,i,iens)); 
      real qvs = 0.622*EVS / ( ref_pres(iens,k)*100. - evs ); //  ! pres(icrm,kk) with unit of hPa
      real alt = 287.0 * temp(k,j,i,icrm) / (100.*ref_pres(icrm,k));
      real rho_t = rho_d(k,j,i,iens) + rho_v(k,j,i,iens)
      L1_sum_rho(k,j,i,iens) += rho_t
      L1_sum_qc (k,j,i,iens) += rho_c(k,j,i,iens) / rho_t;
      L1_sum_qr (k,j,i,iens) += rho_r(k,j,i,iens) / rho_t;
      L1_sum_qi (k,j,i,iens) += rho_i(k,j,i,iens) / rho_t;
      L1_sum_qs (k,j,i,iens) += 0; // P3 does not have a snow category
      L1_sum_ww (k,j,i,iens) += crm_wvel(k,j,i,iens);
      L1_sum_rh (k,j,i,iens) += ( rho_v(k,j,i,iens) / rho_t ) / qvs(k,j,i,iens);
      yakl::atomicAdd( ecpp_L1_sum_cld(k,iens), cldfrac(k,j,i,iens) * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // Do we need to take care of these other variables from the old SAM code?
  // qcloud_bfsum1(:,:,:,icrm) = qcloud_bfsum1(:,:,:,icrm) + qcloud_bf(:,:,:,icrm)
  // L1_sum_qs    (:,:,:,icrm) = L1_sum_qs    (:,:,:,icrm) + qsnow(:,:,:,icrm)
  // qgraupsum1   (:,:,:,icrm) = qgraupsum1   (:,:,:,icrm) + qgraup(:,:,:,icrm)
  // qlsinksum1   (:,:,:,icrm) = qlsinksum1   (:,:,:,icrm) + qlsink(:,:,:,icrm)*qcloud(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg
  // precrsum1    (:,:,:,icrm) = precrsum1    (:,:,:,icrm) + precr(:,:,:,icrm)
  // precsolidsum1(:,:,:,icrm) = precsolidsum1(:,:,:,icrm) + precsolid(:,:,:,icrm)
  // altsum1      (:,:,:,icrm) = altsum1      (:,:,:,icrm) + alt(:,:,:,icrm
  // qlsink_bfsum1(:,:,:,icrm) = qlsink_bfsum1(:,:,:,icrm) + qlsink_bf(:,:,:,icrm)*qcloud_bf(:,:,:,icrm)  ! Note this is converted back in rsum2ToAvg
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) { L1_cnt += 1; });
  //------------------------------------------------------------------------------------------------
}


// convert L1 running sums to averages
inline void pam_ecpp_L1_sum_to_avg( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  auto L1_cnt          = dm_device.get<int, 1>("ecpp_L1_cnt");
  auto L1_sum_cld      = dm_device.get<real,2>("ecpp_L1_sum_cld");
  auto L1_sum_rho      = dm_device.get<real,4>("ecpp_L1_sum_rho");
  auto L1_sum_qc       = dm_device.get<real,4>("ecpp_L1_sum_qc");
  auto L1_sum_qr       = dm_device.get<real,4>("ecpp_L1_sum_qr");
  auto L1_sum_qi       = dm_device.get<real,4>("ecpp_L1_sum_qi");
  auto L1_sum_qs       = dm_device.get<real,4>("ecpp_L1_sum_qs");
  auto L1_sum_ww       = dm_device.get<real,4>("ecpp_L1_sum_ww");
  auto L1_sum_rh       = dm_device.get<real,4>("ecpp_L1_sum_rh");
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<2>(nz, nens), YAKL_LAMBDA (int k, int iens) {
    L1_avg_cld(k,iens) = L1_sum_cld(k,iens) / L1_cnt;
  });
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    L1_avg_rho(k,j,i,iens) = L1_sum_rho(k,j,i,iens) / L1_cnt(iens);
    L1_avg_qc (k,j,i,iens) = L1_sum_qc (k,j,i,iens) / L1_cnt(iens);
    L1_avg_qr (k,j,i,iens) = L1_sum_qr (k,j,i,iens) / L1_cnt(iens);
    L1_avg_qi (k,j,i,iens) = L1_sum_qi (k,j,i,iens) / L1_cnt(iens);
    L1_avg_qs (k,j,i,iens) = L1_sum_qs (k,j,i,iens) / L1_cnt(iens);
    L1_avg_ww (k,j,i,iens) = L1_sum_ww (k,j,i,iens) / L1_cnt(iens);
    L1_avg_rh (k,j,i,iens) = L1_sum_rh (k,j,i,iens) / L1_cnt(iens);
  });
  //------------------------------------------------------------------------------------------------
}


// zero out L1 running sums
inline void pam_ecpp_zero_L2_sums( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  // auto L2_sum_      = dm_device.get<real,4>("ecpp_L2_sum_");
  //------------------------------------------------------------------------------------------------
  // parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
  //   L2_sum_(k,j,i,iens) = 0;
  // });
  //------------------------------------------------------------------------------------------------
}


// update L2 running sums
inline void pam_ecpp_update_L2_sums( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  auto shoc_tk         = dm_device.get<real const,4>("tk");
  //------------------------------------------------------------------------------------------------
  auto L2_cnt          = dm_device.get<int, 1>("ecpp_L2_cnt");
  auto L2_sum_tk       = dm_device.get<real,2>("ecpp_L2_sum_tk");
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    L2_sum_tk(k,j,i,iens) += shoc_tk(k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) { L2_cnt += 1; });
  //------------------------------------------------------------------------------------------------
}


// convert L1 running sums to averages
inline void pam_ecpp_L1_sum_to_avg( pam::PamCoupler &coupler) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device      = coupler.get_data_manager_device_readwrite();
  auto &dm_host        = coupler.get_data_manager_host_readwrite();
  auto nens            = coupler.get_option<int>("ncrms");
  auto nx              = coupler.get_option<int>("crm_nx");
  auto ny              = coupler.get_option<int>("crm_ny");
  auto nz              = coupler.get_option<int>("crm_nz");
  auto gcm_nlev        = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  // auto L2_cnt          = dm_device.get<int, 1>("ecpp_L2_cnt");
  // auto L2_sum_tk       = dm_device.get<real,2>("ecpp_L2_sum_tk");
  //------------------------------------------------------------------------------------------------
  // parallel_for(SimpleBounds<4>(nz, ny, nx, nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
  //   L2_avg_tk(k,j,i,iens) = L2_sum_tk(k,j,i,iens) / L2_cnt(iens);
  // });
  //------------------------------------------------------------------------------------------------
}
