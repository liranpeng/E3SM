#pragma once

#include "pam_coupler.h"
#include "Dycore.h"
#include <vector>
// #include "ecppvars.h"

// allocate and initialize variables
inline void ecpp_crm_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nzm        = coupler.get_option<int>("crm_nz");  // Note that nzm   = crm_nz
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");

  int kbase, ktop;
  int m;
  int nup, ndn, icrm;
  std::string msg;

  int nxstag = nx + 1;
  int nystag = ny + 1;
  int nzstag = nzm + 1;

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

  double upthresh = 1.0;
  double downthresh = 1.0;
  double upthresh2 = 0.5;
  double downthresh2 = 0.5;
  double cloudthresh = 1e-6;
  double prcpthresh = 1e-6;

  double cloudthresh_trans = 1e-5;
  double precthresh_trans = 1e-4;

  int areaavgtype = 1;
  int plumetype = 1;
  bool allcomb = false;

  int nupdraft = 0;
  int ndndraft = 0;
  int ndraft_max = 0;

  int nupdraft_max = 0;  // Note: This variable is initialized but not used in the given code
  int ndndraft_max = 0;  // Note: This variable is initialized but not used in the given code

  int itavg1 = 0;
  int itavg2 = 0; // level-1 and level-2 counters

  // variables should inside ecppvars.h
  int DN1 = 0; // !First index of downward classes
  int NCLASS_TR = 0; // !Num. of transport classes
  int ncc_in       = 2; // Nnumber of clear/cloudy sub-calsses
  int nprcp_in     = 2; // Number of non-precipitating/precipitating sub-classes.
  int NCLASS_CL = ncc_in; // Number of cloud classes
  int NCLASS_PR = nprcp_in; // Number of precipitaion classes

  std::vector<std::vector<int>> updraftbase;
  std::vector<std::vector<int>> updrafttop;
  std::vector<std::vector<int>> dndrafttop;
  std::vector<std::vector<int>> dndraftbase;

  printf("Liran check start ECPP init\n");

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

switch (plumetype) {
    case 1:  // single plume
        nupdraft = 1;
        ndndraft = 1;
        break;
    case 2:  // core and weak plumes
        nupdraft = 2;
        ndndraft = 2;
        break;
    case 3:
        for (int kbase = 1; kbase < nzm; ++kbase) {
            if (allcomb) {  // all possible tops
                nupdraft += nzm - kbase;
            } else {        // one top per base
                nupdraft += 1;
            }
        }
        for (int ktop = nzm; ktop >= 2; --ktop) {
            if (allcomb) {  // all possible bases
                ndndraft += ktop - 1;
            } else {        // one base per top
                ndndraft += 1;
            }
        }
        break;
    // Optionally handle other cases or a default if necessary
}


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

updraftbase.resize(nupdraft_max, std::vector<int>(nens));
updrafttop.resize(nupdraft_max, std::vector<int>(nens));


for (int icrm = 0; icrm < nens; ++icrm) {
    switch (plumetype) {
        case 1:  // single plume
            updraftbase[0][icrm] = 0;
            updrafttop[0][icrm] = nzm - 1;
            dndrafttop[0][icrm] = nzm - 1;
            dndraftbase[0][icrm] = 0;
            break;

        case 2:
            for (int i = 0; i < 2; ++i) {
                updraftbase[i][icrm] = 0;
                updrafttop[i][icrm] = nzm - 1;
                dndrafttop[i][icrm] = nzm - 1;
                dndraftbase[i][icrm] = 0;
            }
            break;

        case 3:
            m = 0;
            for (int kbase = 0; kbase < nzm - 1; ++kbase) {
                if (allcomb) { // loop over all possible tops.
                    for (int ktop = kbase + 1; ktop < nzm; ++ktop) {
                        updraftbase[m][icrm] = kbase;
                        updrafttop[m][icrm] = ktop;
                        m++;
                    }
                } else { // only one top per base
                    updraftbase[m][icrm] = kbase;
                    updrafttop[m][icrm] = nzm - 1;
                    m++;
                }
            }

            m = 0;
            for (int ktop = nzm - 1; ktop >= 1; --ktop) {
                if (allcomb) { // loop over all possible bases.
                    for (int kbase = ktop - 1; kbase >= 0; --kbase) {
                        dndrafttop[m][icrm] = ktop;
                        dndraftbase[m][icrm] = kbase;
                        m++;
                    }
                } else { // only one base per top
                    dndrafttop[m][icrm] = ktop;
                    dndraftbase[m][icrm] = 0;
                    m++;
                }
            }
            break;
    }
}

// 4D vector allocations
std::vector<std::vector<std::vector<std::vector<double>>>> qlsink(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> precr(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> precsolid(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> rh(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qvs(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));

std::vector<std::vector<std::vector<std::vector<double>>>> qlsink_bf(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> prain(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qcloud_bf(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));

std::vector<std::vector<std::vector<std::vector<double>>>> qcloudsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qcloud_bfsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qrainsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qicesum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qsnowsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qgraupsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qlsinksum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> precrsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> precsolidsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> precallsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> altsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> rhsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> cf3dsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> wwsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzstag, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> wwsqsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzstag, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> tkesgssum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qlsink_bfsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> prainsum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));
std::vector<std::vector<std::vector<std::vector<double>>>> qvssum1(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(nzm, std::vector<double>(nens))));

std::vector<std::vector<double>> xkhvsum(nzm, std::vector<double>(nens));

std::vector<std::vector<double>> wwqui_cen_sum(nzm, std::vector<double>(nens));
std::vector<std::vector<double>> wwqui_bnd_sum(nzm+1, std::vector<double>(nens));
std::vector<std::vector<double>> wwqui_cloudy_cen_sum(nzm, std::vector<double>(nens));
std::vector<std::vector<double>> wwqui_cloudy_bnd_sum(nzm+1, std::vector<double>(nens));

std::vector<std::vector<double>> wup_thresh(nzm+1, std::vector<double>(nens));
std::vector<std::vector<double>> wdown_thresh(nzm+1, std::vector<double>(nens));


std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    area_bnd_final(nzstag, std::vector<std::vector<std::vector<std::vector<double>>>>
        (NCLASS_CL, std::vector<std::vector<std::vector<double>>>
            (ndraft_max, std::vector<std::vector<double>>
                (NCLASS_PR, std::vector<double>(nens)))));

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    area_bnd_sum = area_bnd_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    area_cen_final(nzm, std::vector<std::vector<std::vector<std::vector<double>>>>
        (NCLASS_CL, std::vector<std::vector<std::vector<double>>>
            (ndraft_max, std::vector<std::vector<double>>
                (NCLASS_PR, std::vector<double>(nens)))));

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    area_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    mass_bnd_final = area_bnd_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    mass_bnd_sum = area_bnd_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    mass_cen_final = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    mass_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    ent_bnd_sum = area_bnd_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    rh_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qcloud_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qcloud_bf_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qrain_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qice_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qsnow_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qgraup_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qlsink_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    precr_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    precsolid_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    precall_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qlsink_bf_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    qlsink_avg_cen_sum = area_cen_final;

std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
    prain_cen_sum = area_cen_final;


/*

All the variables listed here are already initialized to zero. 
In C++, when you use the std::vector constructor that takes a size and 
a default value (like you've done), it initializes all the elements 
with that default value. In your case, the default value is another 
vector, which itself is constructed with a size and a default value, 
and so on. For the innermost vector, the default value is of type 
double, which is initialized to zero.

No need to call zero_out_sums1

*/



  printf("Liran check start ECPP init end here\n");

}

