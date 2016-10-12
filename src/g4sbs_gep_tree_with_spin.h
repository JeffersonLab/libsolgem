//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan  7 11:54:23 2016 by ROOT version 5.34/32
// from TTree T/Geant4 SBS Simulation
// found on file: gep_spin_transport_Sx.root
//////////////////////////////////////////////////////////

#ifndef __G4SBS_GEP_TREE_WITH_SPIN_H
#define __G4SBS_GEP_TREE_WITH_SPIN_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
// Fixed size dimensions of array or collections stored in the TTree if any.

class g4sbs_gep_tree_with_spin {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   bool            fEcalBox;// needed to turn on/off the reading of the H(E)CAL_box data

   // Declaration of leaf types
   Double_t        ev_count;
   Double_t        ev_rate;
   Double_t        ev_solang;
   Double_t        ev_sigma;
   Double_t        ev_W2;
   Double_t        ev_xbj;
   Double_t        ev_Q2;
   Double_t        ev_th;
   Double_t        ev_ph;
   Double_t        ev_Aperp;
   Double_t        ev_Apar;
   Double_t        ev_Pt;
   Double_t        ev_Pl;
   Double_t        ev_vx;
   Double_t        ev_vy;
   Double_t        ev_vz;
   Double_t        ev_ep;
   Double_t        ev_np;
   Double_t        ev_epx;
   Double_t        ev_epy;
   Double_t        ev_epz;
   Double_t        ev_npx;
   Double_t        ev_npy;
   Double_t        ev_npz;
   Double_t        ev_nth;
   Double_t        ev_nph;
   Double_t        ev_pmperp;
   Double_t        ev_pmpar;
   Double_t        ev_pmparsm;
   Double_t        ev_z;
   Double_t        ev_phperp;
   Double_t        ev_phih;
   Double_t        ev_MX2;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   Double_t        gen_thbb;
   Double_t        gen_thsbs;
   Double_t        gen_dbb;
   Double_t        gen_dsbs;
   Double_t        gen_dhcal;
   Double_t        gen_drich;
   Double_t        gen_dsbstrkr;
   Double_t        gen_Ebeam;
   
   Int_t           Earm_CDET_hit_nhits;
   std::vector<int>     *Earm_CDET_hit_PMT;
   std::vector<int>     *Earm_CDET_hit_row;
   std::vector<int>     *Earm_CDET_hit_col;
   std::vector<int>     *Earm_CDET_hit_plane;
   std::vector<double>  *Earm_CDET_hit_xcell;
   std::vector<double>  *Earm_CDET_hit_ycell;
   std::vector<double>  *Earm_CDET_hit_zcell;
   std::vector<double>  *Earm_CDET_hit_xgcell;
   std::vector<double>  *Earm_CDET_hit_ygcell;
   std::vector<double>  *Earm_CDET_hit_zgcell;
   std::vector<int>     *Earm_CDET_hit_NumPhotoelectrons;
   std::vector<double>  *Earm_CDET_hit_Time_avg;
   std::vector<double>  *Earm_CDET_hit_Time_rms;
   std::vector<double>  *Earm_CDET_hit_Time_min;
   std::vector<double>  *Earm_CDET_hit_Time_max;
   
   Int_t           Earm_CDET_Scint_hit_nhits;
   std::vector<int>     *Earm_CDET_Scint_hit_row;
   std::vector<int>     *Earm_CDET_Scint_hit_col;
   std::vector<int>     *Earm_CDET_Scint_hit_cell;
   std::vector<int>     *Earm_CDET_Scint_hit_plane;
   std::vector<double>  *Earm_CDET_Scint_hit_xcell;
   std::vector<double>  *Earm_CDET_Scint_hit_ycell;
   std::vector<double>  *Earm_CDET_Scint_hit_zcell;
   std::vector<double>  *Earm_CDET_Scint_hit_xcellg;
   std::vector<double>  *Earm_CDET_Scint_hit_ycellg;
   std::vector<double>  *Earm_CDET_Scint_hit_zcellg;
   std::vector<double>  *Earm_CDET_Scint_hit_xhit;
   std::vector<double>  *Earm_CDET_Scint_hit_yhit;
   std::vector<double>  *Earm_CDET_Scint_hit_zhit;
   std::vector<double>  *Earm_CDET_Scint_hit_sumedep;
   std::vector<double>  *Earm_CDET_Scint_hit_tavg;
   std::vector<double>  *Earm_CDET_Scint_hit_trms;
   std::vector<double>  *Earm_CDET_Scint_hit_tmin;
   std::vector<double>  *Earm_CDET_Scint_hit_tmax;
   
   Int_t           Earm_ECAL_hit_nhits;
   std::vector<int>     *Earm_ECAL_hit_PMT;
   std::vector<int>     *Earm_ECAL_hit_row;
   std::vector<int>     *Earm_ECAL_hit_col;
   std::vector<int>     *Earm_ECAL_hit_plane;
   std::vector<double>  *Earm_ECAL_hit_xcell;
   std::vector<double>  *Earm_ECAL_hit_ycell;
   std::vector<double>  *Earm_ECAL_hit_zcell;
   std::vector<double>  *Earm_ECAL_hit_xgcell;
   std::vector<double>  *Earm_ECAL_hit_ygcell;
   std::vector<double>  *Earm_ECAL_hit_zgcell;
   std::vector<int>     *Earm_ECAL_hit_NumPhotoelectrons;
   std::vector<double>  *Earm_ECAL_hit_Time_avg;
   std::vector<double>  *Earm_ECAL_hit_Time_rms;
   std::vector<double>  *Earm_ECAL_hit_Time_min;
   std::vector<double>  *Earm_ECAL_hit_Time_max;
 
   Int_t           Earm_ECAL_box_hit_nhits;
   std::vector<int>     *Earm_ECAL_box_hit_row;
   std::vector<int>     *Earm_ECAL_box_hit_col;
   std::vector<int>     *Earm_ECAL_box_hit_cell;
   std::vector<int>     *Earm_ECAL_box_hit_plane;
   std::vector<double>  *Earm_ECAL_box_hit_xcell;
   std::vector<double>  *Earm_ECAL_box_hit_ycell;
   std::vector<double>  *Earm_ECAL_box_hit_zcell;
   std::vector<double>  *Earm_ECAL_box_hit_xcellg;
   std::vector<double>  *Earm_ECAL_box_hit_ycellg;
   std::vector<double>  *Earm_ECAL_box_hit_zcellg;
   std::vector<double>  *Earm_ECAL_box_hit_xhit;
   std::vector<double>  *Earm_ECAL_box_hit_yhit;
   std::vector<double>  *Earm_ECAL_box_hit_zhit;
   std::vector<double>  *Earm_ECAL_box_hit_sumedep;
   std::vector<double>  *Earm_ECAL_box_hit_tavg;
   std::vector<double>  *Earm_ECAL_box_hit_trms;
   std::vector<double>  *Earm_ECAL_box_hit_tmin;
   std::vector<double>  *Earm_ECAL_box_hit_tmax;
   
   Int_t           Earm_ECalTF1_hit_nhits;
   std::vector<int>     *Earm_ECalTF1_hit_row;
   std::vector<int>     *Earm_ECalTF1_hit_col;
   std::vector<int>     *Earm_ECalTF1_hit_cell;
   std::vector<int>     *Earm_ECalTF1_hit_plane;
   std::vector<double>  *Earm_ECalTF1_hit_xcell;
   std::vector<double>  *Earm_ECalTF1_hit_ycell;
   std::vector<double>  *Earm_ECalTF1_hit_zcell;
   std::vector<double>  *Earm_ECalTF1_hit_xcellg;
   std::vector<double>  *Earm_ECalTF1_hit_ycellg;
   std::vector<double>  *Earm_ECalTF1_hit_zcellg;
   std::vector<double>  *Earm_ECalTF1_hit_xhit;
   std::vector<double>  *Earm_ECalTF1_hit_yhit;
   std::vector<double>  *Earm_ECalTF1_hit_zhit;
   std::vector<double>  *Earm_ECalTF1_hit_sumedep;
   std::vector<double>  *Earm_ECalTF1_hit_tavg;
   std::vector<double>  *Earm_ECalTF1_hit_trms;
   std::vector<double>  *Earm_ECalTF1_hit_tmin;
   std::vector<double>  *Earm_ECalTF1_hit_tmax;
   
   Int_t           Harm_FPP1_hit_nhits;
   std::vector<int>     *Harm_FPP1_hit_plane;
   std::vector<int>     *Harm_FPP1_hit_strip;
   std::vector<double>  *Harm_FPP1_hit_x;
   std::vector<double>  *Harm_FPP1_hit_y;
   std::vector<double>  *Harm_FPP1_hit_z;
   std::vector<double>  *Harm_FPP1_hit_polx;
   std::vector<double>  *Harm_FPP1_hit_poly;
   std::vector<double>  *Harm_FPP1_hit_polz;
   std::vector<double>  *Harm_FPP1_hit_t;
   std::vector<double>  *Harm_FPP1_hit_trms;
   std::vector<double>  *Harm_FPP1_hit_tmin;
   std::vector<double>  *Harm_FPP1_hit_tmax;
   std::vector<double>  *Harm_FPP1_hit_tx;
   std::vector<double>  *Harm_FPP1_hit_ty;
   std::vector<double>  *Harm_FPP1_hit_txp;
   std::vector<double>  *Harm_FPP1_hit_typ;
   std::vector<int>     *Harm_FPP1_hit_trid;
   std::vector<int>     *Harm_FPP1_hit_mid;
   std::vector<int>     *Harm_FPP1_hit_pid;
   std::vector<double>  *Harm_FPP1_hit_vx;
   std::vector<double>  *Harm_FPP1_hit_vy;
   std::vector<double>  *Harm_FPP1_hit_vz;
   std::vector<double>  *Harm_FPP1_hit_p;
   std::vector<double>  *Harm_FPP1_hit_edep;
   std::vector<double>  *Harm_FPP1_hit_beta;
   
   Int_t           Harm_FPP1_Track_ntracks;
   std::vector<int>     *Harm_FPP1_Track_TID;
   std::vector<int>     *Harm_FPP1_Track_PID;
   std::vector<int>     *Harm_FPP1_Track_MID;
   std::vector<int>     *Harm_FPP1_Track_NumHits;
   std::vector<int>     *Harm_FPP1_Track_NumPlanes;
   std::vector<int>     *Harm_FPP1_Track_NDF;
   std::vector<double>  *Harm_FPP1_Track_Chi2fit;
   std::vector<double>  *Harm_FPP1_Track_Chi2true;
   std::vector<double>  *Harm_FPP1_Track_X;
   std::vector<double>  *Harm_FPP1_Track_Y;
   std::vector<double>  *Harm_FPP1_Track_Xp;
   std::vector<double>  *Harm_FPP1_Track_Yp;
   std::vector<double>  *Harm_FPP1_Track_T;
   std::vector<double>  *Harm_FPP1_Track_P;
   std::vector<double>  *Harm_FPP1_Track_Sx;
   std::vector<double>  *Harm_FPP1_Track_Sy;
   std::vector<double>  *Harm_FPP1_Track_Sz;
   std::vector<double>  *Harm_FPP1_Track_Xfit;
   std::vector<double>  *Harm_FPP1_Track_Yfit;
   std::vector<double>  *Harm_FPP1_Track_Xpfit;
   std::vector<double>  *Harm_FPP1_Track_Ypfit;
   
   Int_t           Harm_FPP2_hit_nhits;
   std::vector<int>     *Harm_FPP2_hit_plane;
   std::vector<int>     *Harm_FPP2_hit_strip;
   std::vector<double>  *Harm_FPP2_hit_x;
   std::vector<double>  *Harm_FPP2_hit_y;
   std::vector<double>  *Harm_FPP2_hit_z;
   std::vector<double>  *Harm_FPP2_hit_polx;
   std::vector<double>  *Harm_FPP2_hit_poly;
   std::vector<double>  *Harm_FPP2_hit_polz;
   std::vector<double>  *Harm_FPP2_hit_t;
   std::vector<double>  *Harm_FPP2_hit_trms;
   std::vector<double>  *Harm_FPP2_hit_tmin;
   std::vector<double>  *Harm_FPP2_hit_tmax;
   std::vector<double>  *Harm_FPP2_hit_tx;
   std::vector<double>  *Harm_FPP2_hit_ty;
   std::vector<double>  *Harm_FPP2_hit_txp;
   std::vector<double>  *Harm_FPP2_hit_typ;
   std::vector<int>     *Harm_FPP2_hit_trid;
   std::vector<int>     *Harm_FPP2_hit_mid;
   std::vector<int>     *Harm_FPP2_hit_pid;
   std::vector<double>  *Harm_FPP2_hit_vx;
   std::vector<double>  *Harm_FPP2_hit_vy;
   std::vector<double>  *Harm_FPP2_hit_vz;
   std::vector<double>  *Harm_FPP2_hit_p;
   std::vector<double>  *Harm_FPP2_hit_edep;
   std::vector<double>  *Harm_FPP2_hit_beta;
   
   Int_t           Harm_FPP2_Track_ntracks;
   std::vector<int>     *Harm_FPP2_Track_TID;
   std::vector<int>     *Harm_FPP2_Track_PID;
   std::vector<int>     *Harm_FPP2_Track_MID;
   std::vector<int>     *Harm_FPP2_Track_NumHits;
   std::vector<int>     *Harm_FPP2_Track_NumPlanes;
   std::vector<int>     *Harm_FPP2_Track_NDF;
   std::vector<double>  *Harm_FPP2_Track_Chi2fit;
   std::vector<double>  *Harm_FPP2_Track_Chi2true;
   std::vector<double>  *Harm_FPP2_Track_X;
   std::vector<double>  *Harm_FPP2_Track_Y;
   std::vector<double>  *Harm_FPP2_Track_Xp;
   std::vector<double>  *Harm_FPP2_Track_Yp;
   std::vector<double>  *Harm_FPP2_Track_T;
   std::vector<double>  *Harm_FPP2_Track_P;
   std::vector<double>  *Harm_FPP2_Track_Sx;
   std::vector<double>  *Harm_FPP2_Track_Sy;
   std::vector<double>  *Harm_FPP2_Track_Sz;
   std::vector<double>  *Harm_FPP2_Track_Xfit;
   std::vector<double>  *Harm_FPP2_Track_Yfit;
   std::vector<double>  *Harm_FPP2_Track_Xpfit;
   std::vector<double>  *Harm_FPP2_Track_Ypfit;
   
   Int_t           Harm_FT_hit_nhits;
   std::vector<int>     *Harm_FT_hit_plane;
   std::vector<int>     *Harm_FT_hit_strip;
   std::vector<double>  *Harm_FT_hit_x;
   std::vector<double>  *Harm_FT_hit_y;
   std::vector<double>  *Harm_FT_hit_z;
   std::vector<double>  *Harm_FT_hit_polx;
   std::vector<double>  *Harm_FT_hit_poly;
   std::vector<double>  *Harm_FT_hit_polz;
   std::vector<double>  *Harm_FT_hit_t;
   std::vector<double>  *Harm_FT_hit_trms;
   std::vector<double>  *Harm_FT_hit_tmin;
   std::vector<double>  *Harm_FT_hit_tmax;
   std::vector<double>  *Harm_FT_hit_tx;
   std::vector<double>  *Harm_FT_hit_ty;
   std::vector<double>  *Harm_FT_hit_txp;
   std::vector<double>  *Harm_FT_hit_typ;
   std::vector<int>     *Harm_FT_hit_trid;
   std::vector<int>     *Harm_FT_hit_mid;
   std::vector<int>     *Harm_FT_hit_pid;
   std::vector<double>  *Harm_FT_hit_vx;
   std::vector<double>  *Harm_FT_hit_vy;
   std::vector<double>  *Harm_FT_hit_vz;
   std::vector<double>  *Harm_FT_hit_p;
   std::vector<double>  *Harm_FT_hit_edep;
   std::vector<double>  *Harm_FT_hit_beta;
   
   Int_t           Harm_FT_Track_ntracks;
   std::vector<int>     *Harm_FT_Track_TID;
   std::vector<int>     *Harm_FT_Track_PID;
   std::vector<int>     *Harm_FT_Track_MID;
   std::vector<int>     *Harm_FT_Track_NumHits;
   std::vector<int>     *Harm_FT_Track_NumPlanes;
   std::vector<int>     *Harm_FT_Track_NDF;
   std::vector<double>  *Harm_FT_Track_Chi2fit;
   std::vector<double>  *Harm_FT_Track_Chi2true;
   std::vector<double>  *Harm_FT_Track_X;
   std::vector<double>  *Harm_FT_Track_Y;
   std::vector<double>  *Harm_FT_Track_Xp;
   std::vector<double>  *Harm_FT_Track_Yp;
   std::vector<double>  *Harm_FT_Track_T;
   std::vector<double>  *Harm_FT_Track_P;
   std::vector<double>  *Harm_FT_Track_Sx;
   std::vector<double>  *Harm_FT_Track_Sy;
   std::vector<double>  *Harm_FT_Track_Sz;
   std::vector<double>  *Harm_FT_Track_Xfit;
   std::vector<double>  *Harm_FT_Track_Yfit;
   std::vector<double>  *Harm_FT_Track_Xpfit;
   std::vector<double>  *Harm_FT_Track_Ypfit;
   
   Int_t           Harm_HCAL_box_hit_nhits;
   std::vector<int>     *Harm_HCAL_box_hit_row;
   std::vector<int>     *Harm_HCAL_box_hit_col;
   std::vector<int>     *Harm_HCAL_box_hit_cell;
   std::vector<int>     *Harm_HCAL_box_hit_plane;
   std::vector<double>  *Harm_HCAL_box_hit_xcell;
   std::vector<double>  *Harm_HCAL_box_hit_ycell;
   std::vector<double>  *Harm_HCAL_box_hit_zcell;
   std::vector<double>  *Harm_HCAL_box_hit_xcellg;
   std::vector<double>  *Harm_HCAL_box_hit_ycellg;
   std::vector<double>  *Harm_HCAL_box_hit_zcellg;
   std::vector<double>  *Harm_HCAL_box_hit_xhit;
   std::vector<double>  *Harm_HCAL_box_hit_yhit;
   std::vector<double>  *Harm_HCAL_box_hit_zhit;
   std::vector<double>  *Harm_HCAL_box_hit_sumedep;
   std::vector<double>  *Harm_HCAL_box_hit_tavg;
   std::vector<double>  *Harm_HCAL_box_hit_trms;
   std::vector<double>  *Harm_HCAL_box_hit_tmin;
   std::vector<double>  *Harm_HCAL_box_hit_tmax;
   
   Int_t           Harm_HCal_hit_nhits;
   std::vector<int>     *Harm_HCal_hit_PMT;
   std::vector<int>     *Harm_HCal_hit_row;
   std::vector<int>     *Harm_HCal_hit_col;
   std::vector<int>     *Harm_HCal_hit_plane;
   std::vector<double>  *Harm_HCal_hit_xcell;
   std::vector<double>  *Harm_HCal_hit_ycell;
   std::vector<double>  *Harm_HCal_hit_zcell;
   std::vector<double>  *Harm_HCal_hit_xgcell;
   std::vector<double>  *Harm_HCal_hit_ygcell;
   std::vector<double>  *Harm_HCal_hit_zgcell;
   std::vector<int>     *Harm_HCal_hit_NumPhotoelectrons;
   std::vector<double>  *Harm_HCal_hit_Time_avg;
   std::vector<double>  *Harm_HCal_hit_Time_rms;
   std::vector<double>  *Harm_HCal_hit_Time_min;
   std::vector<double>  *Harm_HCal_hit_Time_max;
   
   Int_t           Harm_HCalScint_hit_nhits;
   std::vector<int>     *Harm_HCalScint_hit_row;
   std::vector<int>     *Harm_HCalScint_hit_col;
   std::vector<int>     *Harm_HCalScint_hit_cell;
   std::vector<int>     *Harm_HCalScint_hit_plane;
   std::vector<double>  *Harm_HCalScint_hit_xcell;
   std::vector<double>  *Harm_HCalScint_hit_ycell;
   std::vector<double>  *Harm_HCalScint_hit_zcell;
   std::vector<double>  *Harm_HCalScint_hit_xcellg;
   std::vector<double>  *Harm_HCalScint_hit_ycellg;
   std::vector<double>  *Harm_HCalScint_hit_zcellg;
   std::vector<double>  *Harm_HCalScint_hit_xhit;
   std::vector<double>  *Harm_HCalScint_hit_yhit;
   std::vector<double>  *Harm_HCalScint_hit_zhit;
   std::vector<double>  *Harm_HCalScint_hit_sumedep;
   std::vector<double>  *Harm_HCalScint_hit_tavg;
   std::vector<double>  *Harm_HCalScint_hit_trms;
   std::vector<double>  *Harm_HCalScint_hit_tmin;
   std::vector<double>  *Harm_HCalScint_hit_tmax;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_gen;   //!
   
   TBranch        *b_Earm_CDET_hit_nhits;   //!
   TBranch        *b_Earm_CDET_hit_PMT;   //!
   TBranch        *b_Earm_CDET_hit_row;   //!
   TBranch        *b_Earm_CDET_hit_col;   //!
   TBranch        *b_Earm_CDET_hit_plane;   //!
   TBranch        *b_Earm_CDET_hit_xcell;   //!
   TBranch        *b_Earm_CDET_hit_ycell;   //!
   TBranch        *b_Earm_CDET_hit_zcell;   //!
   TBranch        *b_Earm_CDET_hit_xgcell;   //!
   TBranch        *b_Earm_CDET_hit_ygcell;   //!
   TBranch        *b_Earm_CDET_hit_zgcell;   //!
   TBranch        *b_Earm_CDET_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_CDET_hit_Time_avg;   //!
   TBranch        *b_Earm_CDET_hit_Time_rms;   //!
   TBranch        *b_Earm_CDET_hit_Time_min;   //!
   TBranch        *b_Earm_CDET_hit_Time_max;   //!
   
   TBranch        *b_Earm_CDET_Scint_hit_nhits;   //!
   TBranch        *b_Earm_CDET_Scint_hit_row;   //!
   TBranch        *b_Earm_CDET_Scint_hit_col;   //!
   TBranch        *b_Earm_CDET_Scint_hit_cell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_plane;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_yhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_sumedep;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tavg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_trms;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmin;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmax;   //!
   
   TBranch        *b_Earm_ECAL_hit_nhits;   //!
   TBranch        *b_Earm_ECAL_hit_PMT;   //!
   TBranch        *b_Earm_ECAL_hit_row;   //!
   TBranch        *b_Earm_ECAL_hit_col;   //!
   TBranch        *b_Earm_ECAL_hit_plane;   //!
   TBranch        *b_Earm_ECAL_hit_xcell;   //!
   TBranch        *b_Earm_ECAL_hit_ycell;   //!
   TBranch        *b_Earm_ECAL_hit_zcell;   //!
   TBranch        *b_Earm_ECAL_hit_xgcell;   //!
   TBranch        *b_Earm_ECAL_hit_ygcell;   //!
   TBranch        *b_Earm_ECAL_hit_zgcell;   //!
   TBranch        *b_Earm_ECAL_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_ECAL_hit_Time_avg;   //!
   TBranch        *b_Earm_ECAL_hit_Time_rms;   //!
   TBranch        *b_Earm_ECAL_hit_Time_min;   //!
   TBranch        *b_Earm_ECAL_hit_Time_max;   //!
   
   TBranch        *b_Earm_ECAL_box_hit_nhits;   //!
   TBranch        *b_Earm_ECAL_box_hit_row;   //!
   TBranch        *b_Earm_ECAL_box_hit_col;   //!
   TBranch        *b_Earm_ECAL_box_hit_cell;   //!
   TBranch        *b_Earm_ECAL_box_hit_plane;   //!
   TBranch        *b_Earm_ECAL_box_hit_xcell;   //!
   TBranch        *b_Earm_ECAL_box_hit_ycell;   //!
   TBranch        *b_Earm_ECAL_box_hit_zcell;   //!
   TBranch        *b_Earm_ECAL_box_hit_xcellg;   //!
   TBranch        *b_Earm_ECAL_box_hit_ycellg;   //!
   TBranch        *b_Earm_ECAL_box_hit_zcellg;   //!
   TBranch        *b_Earm_ECAL_box_hit_xhit;   //!
   TBranch        *b_Earm_ECAL_box_hit_yhit;   //!
   TBranch        *b_Earm_ECAL_box_hit_zhit;   //!
   TBranch        *b_Earm_ECAL_box_hit_sumedep;   //!
   TBranch        *b_Earm_ECAL_box_hit_tavg;   //!
   TBranch        *b_Earm_ECAL_box_hit_trms;   //!
   TBranch        *b_Earm_ECAL_box_hit_tmin;   //!
   TBranch        *b_Earm_ECAL_box_hit_tmax;   //!
   
   TBranch        *b_Earm_ECalTF1_hit_nhits;   //!
   TBranch        *b_Earm_ECalTF1_hit_row;   //!
   TBranch        *b_Earm_ECalTF1_hit_col;   //!
   TBranch        *b_Earm_ECalTF1_hit_cell;   //!
   TBranch        *b_Earm_ECalTF1_hit_plane;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycell;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_xhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_yhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_zhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_sumedep;   //!
   TBranch        *b_Earm_ECalTF1_hit_tavg;   //!
   TBranch        *b_Earm_ECalTF1_hit_trms;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmin;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmax;   //!
   
   TBranch        *b_Harm_FPP1_hit_nhits;   //!
   TBranch        *b_Harm_FPP1_hit_plane;   //!
   TBranch        *b_Harm_FPP1_hit_strip;   //!
   TBranch        *b_Harm_FPP1_hit_x;   //!
   TBranch        *b_Harm_FPP1_hit_y;   //!
   TBranch        *b_Harm_FPP1_hit_z;   //!
   TBranch        *b_Harm_FPP1_hit_polx;   //!
   TBranch        *b_Harm_FPP1_hit_poly;   //!
   TBranch        *b_Harm_FPP1_hit_polz;   //!
   TBranch        *b_Harm_FPP1_hit_t;   //!
   TBranch        *b_Harm_FPP1_hit_trms;   //!
   TBranch        *b_Harm_FPP1_hit_tmin;   //!
   TBranch        *b_Harm_FPP1_hit_tmax;   //!
   TBranch        *b_Harm_FPP1_hit_tx;   //!
   TBranch        *b_Harm_FPP1_hit_ty;   //!
   TBranch        *b_Harm_FPP1_hit_txp;   //!
   TBranch        *b_Harm_FPP1_hit_typ;   //!
   TBranch        *b_Harm_FPP1_hit_trid;   //!
   TBranch        *b_Harm_FPP1_hit_mid;   //!
   TBranch        *b_Harm_FPP1_hit_pid;   //!
   TBranch        *b_Harm_FPP1_hit_vx;   //!
   TBranch        *b_Harm_FPP1_hit_vy;   //!
   TBranch        *b_Harm_FPP1_hit_vz;   //!
   TBranch        *b_Harm_FPP1_hit_p;   //!
   TBranch        *b_Harm_FPP1_hit_edep;   //!
   TBranch        *b_Harm_FPP1_hit_beta;   //!
   
   TBranch        *b_Harm_FPP1_Track_ntracks;   //!
   TBranch        *b_Harm_FPP1_Track_TID;   //!
   TBranch        *b_Harm_FPP1_Track_PID;   //!
   TBranch        *b_Harm_FPP1_Track_MID;   //!
   TBranch        *b_Harm_FPP1_Track_NumHits;   //!
   TBranch        *b_Harm_FPP1_Track_NumPlanes;   //!
   TBranch        *b_Harm_FPP1_Track_NDF;   //!
   TBranch        *b_Harm_FPP1_Track_Chi2fit;   //!
   TBranch        *b_Harm_FPP1_Track_Chi2true;   //!
   TBranch        *b_Harm_FPP1_Track_X;   //!
   TBranch        *b_Harm_FPP1_Track_Y;   //!
   TBranch        *b_Harm_FPP1_Track_Xp;   //!
   TBranch        *b_Harm_FPP1_Track_Yp;   //!
   TBranch        *b_Harm_FPP1_Track_T;   //!
   TBranch        *b_Harm_FPP1_Track_P;   //!
   TBranch        *b_Harm_FPP1_Track_Sx;   //!
   TBranch        *b_Harm_FPP1_Track_Sy;   //!
   TBranch        *b_Harm_FPP1_Track_Sz;   //!
   TBranch        *b_Harm_FPP1_Track_Xfit;   //!
   TBranch        *b_Harm_FPP1_Track_Yfit;   //!
   TBranch        *b_Harm_FPP1_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP1_Track_Ypfit;   //!
   
   TBranch        *b_Harm_FPP2_hit_nhits;   //!
   TBranch        *b_Harm_FPP2_hit_plane;   //!
   TBranch        *b_Harm_FPP2_hit_strip;   //!
   TBranch        *b_Harm_FPP2_hit_x;   //!
   TBranch        *b_Harm_FPP2_hit_y;   //!
   TBranch        *b_Harm_FPP2_hit_z;   //!
   TBranch        *b_Harm_FPP2_hit_polx;   //!
   TBranch        *b_Harm_FPP2_hit_poly;   //!
   TBranch        *b_Harm_FPP2_hit_polz;   //!
   TBranch        *b_Harm_FPP2_hit_t;   //!
   TBranch        *b_Harm_FPP2_hit_trms;   //!
   TBranch        *b_Harm_FPP2_hit_tmin;   //!
   TBranch        *b_Harm_FPP2_hit_tmax;   //!
   TBranch        *b_Harm_FPP2_hit_tx;   //!
   TBranch        *b_Harm_FPP2_hit_ty;   //!
   TBranch        *b_Harm_FPP2_hit_txp;   //!
   TBranch        *b_Harm_FPP2_hit_typ;   //!
   TBranch        *b_Harm_FPP2_hit_trid;   //!
   TBranch        *b_Harm_FPP2_hit_mid;   //!
   TBranch        *b_Harm_FPP2_hit_pid;   //!
   TBranch        *b_Harm_FPP2_hit_vx;   //!
   TBranch        *b_Harm_FPP2_hit_vy;   //!
   TBranch        *b_Harm_FPP2_hit_vz;   //!
   TBranch        *b_Harm_FPP2_hit_p;   //!
   TBranch        *b_Harm_FPP2_hit_edep;   //!
   TBranch        *b_Harm_FPP2_hit_beta;   //!
   
   TBranch        *b_Harm_FPP2_Track_ntracks;   //!
   TBranch        *b_Harm_FPP2_Track_TID;   //!
   TBranch        *b_Harm_FPP2_Track_PID;   //!
   TBranch        *b_Harm_FPP2_Track_MID;   //!
   TBranch        *b_Harm_FPP2_Track_NumHits;   //!
   TBranch        *b_Harm_FPP2_Track_NumPlanes;   //!
   TBranch        *b_Harm_FPP2_Track_NDF;   //!
   TBranch        *b_Harm_FPP2_Track_Chi2fit;   //!
   TBranch        *b_Harm_FPP2_Track_Chi2true;   //!
   TBranch        *b_Harm_FPP2_Track_X;   //!
   TBranch        *b_Harm_FPP2_Track_Y;   //!
   TBranch        *b_Harm_FPP2_Track_Xp;   //!
   TBranch        *b_Harm_FPP2_Track_Yp;   //!
   TBranch        *b_Harm_FPP2_Track_T;   //!
   TBranch        *b_Harm_FPP2_Track_P;   //!
   TBranch        *b_Harm_FPP2_Track_Sx;   //!
   TBranch        *b_Harm_FPP2_Track_Sy;   //!
   TBranch        *b_Harm_FPP2_Track_Sz;   //!
   TBranch        *b_Harm_FPP2_Track_Xfit;   //!
   TBranch        *b_Harm_FPP2_Track_Yfit;   //!
   TBranch        *b_Harm_FPP2_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP2_Track_Ypfit;   //!
   
   TBranch        *b_Harm_FT_hit_nhits;   //!
   TBranch        *b_Harm_FT_hit_plane;   //!
   TBranch        *b_Harm_FT_hit_strip;   //!
   TBranch        *b_Harm_FT_hit_x;   //!
   TBranch        *b_Harm_FT_hit_y;   //!
   TBranch        *b_Harm_FT_hit_z;   //!
   TBranch        *b_Harm_FT_hit_polx;   //!
   TBranch        *b_Harm_FT_hit_poly;   //!
   TBranch        *b_Harm_FT_hit_polz;   //!
   TBranch        *b_Harm_FT_hit_t;   //!
   TBranch        *b_Harm_FT_hit_trms;   //!
   TBranch        *b_Harm_FT_hit_tmin;   //!
   TBranch        *b_Harm_FT_hit_tmax;   //!
   TBranch        *b_Harm_FT_hit_tx;   //!
   TBranch        *b_Harm_FT_hit_ty;   //!
   TBranch        *b_Harm_FT_hit_txp;   //!
   TBranch        *b_Harm_FT_hit_typ;   //!
   TBranch        *b_Harm_FT_hit_trid;   //!
   TBranch        *b_Harm_FT_hit_mid;   //!
   TBranch        *b_Harm_FT_hit_pid;   //!
   TBranch        *b_Harm_FT_hit_vx;   //!
   TBranch        *b_Harm_FT_hit_vy;   //!
   TBranch        *b_Harm_FT_hit_vz;   //!
   TBranch        *b_Harm_FT_hit_p;   //!
   TBranch        *b_Harm_FT_hit_edep;   //!
   TBranch        *b_Harm_FT_hit_beta;   //!
   TBranch        *b_Harm_FT_Track_ntracks;   //!
   TBranch        *b_Harm_FT_Track_TID;   //!
   TBranch        *b_Harm_FT_Track_PID;   //!
   TBranch        *b_Harm_FT_Track_MID;   //!
   TBranch        *b_Harm_FT_Track_NumHits;   //!
   TBranch        *b_Harm_FT_Track_NumPlanes;   //!
   TBranch        *b_Harm_FT_Track_NDF;   //!
   TBranch        *b_Harm_FT_Track_Chi2fit;   //!
   TBranch        *b_Harm_FT_Track_Chi2true;   //!
   TBranch        *b_Harm_FT_Track_X;   //!
   TBranch        *b_Harm_FT_Track_Y;   //!
   TBranch        *b_Harm_FT_Track_Xp;   //!
   TBranch        *b_Harm_FT_Track_Yp;   //!
   TBranch        *b_Harm_FT_Track_T;   //!
   TBranch        *b_Harm_FT_Track_P;   //!
   TBranch        *b_Harm_FT_Track_Sx;   //!
   TBranch        *b_Harm_FT_Track_Sy;   //!
   TBranch        *b_Harm_FT_Track_Sz;   //!
   TBranch        *b_Harm_FT_Track_Xfit;   //!
   TBranch        *b_Harm_FT_Track_Yfit;   //!
   TBranch        *b_Harm_FT_Track_Xpfit;   //!
   TBranch        *b_Harm_FT_Track_Ypfit;   //!
   
   TBranch        *b_Harm_HCAL_box_hit_nhits;   //!
   TBranch        *b_Harm_HCAL_box_hit_row;   //!
   TBranch        *b_Harm_HCAL_box_hit_col;   //!
   TBranch        *b_Harm_HCAL_box_hit_cell;   //!
   TBranch        *b_Harm_HCAL_box_hit_plane;   //!
   TBranch        *b_Harm_HCAL_box_hit_xcell;   //!
   TBranch        *b_Harm_HCAL_box_hit_ycell;   //!
   TBranch        *b_Harm_HCAL_box_hit_zcell;   //!
   TBranch        *b_Harm_HCAL_box_hit_xcellg;   //!
   TBranch        *b_Harm_HCAL_box_hit_ycellg;   //!
   TBranch        *b_Harm_HCAL_box_hit_zcellg;   //!
   TBranch        *b_Harm_HCAL_box_hit_xhit;   //!
   TBranch        *b_Harm_HCAL_box_hit_yhit;   //!
   TBranch        *b_Harm_HCAL_box_hit_zhit;   //!
   TBranch        *b_Harm_HCAL_box_hit_sumedep;   //!
   TBranch        *b_Harm_HCAL_box_hit_tavg;   //!
   TBranch        *b_Harm_HCAL_box_hit_trms;   //!
   TBranch        *b_Harm_HCAL_box_hit_tmin;   //!
   TBranch        *b_Harm_HCAL_box_hit_tmax;   //!
   
   TBranch        *b_Harm_HCal_hit_nhits;   //!
   TBranch        *b_Harm_HCal_hit_PMT;   //!
   TBranch        *b_Harm_HCal_hit_row;   //!
   TBranch        *b_Harm_HCal_hit_col;   //!
   TBranch        *b_Harm_HCal_hit_plane;   //!
   TBranch        *b_Harm_HCal_hit_xcell;   //!
   TBranch        *b_Harm_HCal_hit_ycell;   //!
   TBranch        *b_Harm_HCal_hit_zcell;   //!
   TBranch        *b_Harm_HCal_hit_xgcell;   //!
   TBranch        *b_Harm_HCal_hit_ygcell;   //!
   TBranch        *b_Harm_HCal_hit_zgcell;   //!
   TBranch        *b_Harm_HCal_hit_NumPhotoelectrons;   //!
   TBranch        *b_Harm_HCal_hit_Time_avg;   //!
   TBranch        *b_Harm_HCal_hit_Time_rms;   //!
   TBranch        *b_Harm_HCal_hit_Time_min;   //!
   TBranch        *b_Harm_HCal_hit_Time_max;   //!
   
   TBranch        *b_Harm_HCalScint_hit_nhits;   //!
   TBranch        *b_Harm_HCalScint_hit_row;   //!
   TBranch        *b_Harm_HCalScint_hit_col;   //!
   TBranch        *b_Harm_HCalScint_hit_cell;   //!
   TBranch        *b_Harm_HCalScint_hit_plane;   //!
   TBranch        *b_Harm_HCalScint_hit_xcell;   //!
   TBranch        *b_Harm_HCalScint_hit_ycell;   //!
   TBranch        *b_Harm_HCalScint_hit_zcell;   //!
   TBranch        *b_Harm_HCalScint_hit_xcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_ycellg;   //!
   TBranch        *b_Harm_HCalScint_hit_zcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_xhit;   //!
   TBranch        *b_Harm_HCalScint_hit_yhit;   //!
   TBranch        *b_Harm_HCalScint_hit_zhit;   //!
   TBranch        *b_Harm_HCalScint_hit_sumedep;   //!
   TBranch        *b_Harm_HCalScint_hit_tavg;   //!
   TBranch        *b_Harm_HCalScint_hit_trms;   //!
   TBranch        *b_Harm_HCalScint_hit_tmin;   //!
   TBranch        *b_Harm_HCalScint_hit_tmax;   //!

   g4sbs_gep_tree_with_spin(TTree *tree=0, bool ecalbox = true);
   virtual ~g4sbs_gep_tree_with_spin();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif



