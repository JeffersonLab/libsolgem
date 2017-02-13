#ifndef __CINT__

#include "SBSSpec.h"
#include "TSBSSimDecoder.h"
#include "TSolSimFile.h"

#include "THaInterface.h"
#include "THaGlobals.h"
#include "THaTextvars.h"
#include "THaAnalyzer.h"
#include "THaDetector.h"

#include "TSystem.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"

#include <iostream>

#endif

void ReplayMCdigDemo(const char* filename = "digdemo3", const char* detsuf = "fpp",
		     Int_t nevent = -1, Int_t nseg = 0, bool do_cuts = true )
{
  if( nseg < 0 || nseg > 100 ) {
    cerr << "Invalid number of run segments = " << nseg
	 << ", must be 0-100" << endl;
    return;
  }
  bool do_parts = true;
  if( nseg == 0 ) {
    do_parts = false;
    nseg = 1;
  }
  
  gSystem->Load("libsolgem.so");
  gSystem->Load("libTreeSearch-SBS.so");
  
  TSolDBManager* manager = TSolDBManager::GetInstance();
  manager->LoadGeneralInfo(Form("db_generalinfo_%s.dat", detsuf));
  
  // dde = new TSBSSimDecoder();
  // dde->SetCrateMapName("db_sbssim_cratemap.dat");
  
  THaInterface::SetDecoder( TSBSSimDecoder::Class() );
  
  cout << "Reading FPP: " << endl;
  THaApparatus* SBS_FPP = new SBS::SBSSpec( "sbs_fpp", "SBS GEp Focal Plane Polarimeter", 1 );
  gHaApps->Add( SBS_FPP );
  cout << "Just read FPP." << endl;
  cout << "Reading FT: " << endl;
  THaApparatus* SBS_FT = new SBS::SBSSpec( "sbs_ft", "SBS GEp front tracker", 1 );
  gHaApps->Add( SBS_FT );
  cout << "Just read FT." << endl;
  
  TString db_prefix_1 = SBS_FPP->GetName();
  TString db_prefix_2 = SBS_FT->GetName();
  db_prefix_1 += ".tracker";
  gHaTextvars->Add( "DET", db_prefix_1.Data() );
  gHaTextvars->Add( "APP", SBS_FPP->GetName() );
  db_prefix_2 += ".tracker";
  gHaTextvars->Add( "DET", db_prefix_2.Data() );
  gHaTextvars->Add( "APP", SBS_FT->GetName() );
  
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  TString rootfile(Form("%s_%s", filename, detsuf)), infile0(Form("%s_%s", filename, detsuf));
  TString odeffile("sbssim.odef"), cutfile("sbssim.cuts");
  rootfile.Append("_replayed_new.root");
  analyzer->EnableBenchmarks();
  analyzer->SetOutFile(rootfile);
  analyzer->SetOdefFile(odeffile);
  if( do_cuts ) analyzer->SetCutFile(cutfile);
  analyzer->SetSummaryFile(Form("%s_new.sum", filename));
  analyzer->SetCrateMapFileName("sbssim_cratemap");
  
  //static int Nrun = TMath::Max(nseg,1);
  THaRunBase* run[0];
  TString title0 = "Digitized MC data";
  for( int i=0; i<nseg; ++i ) {
    TString title(title0), infile(infile0);
    if( do_parts ) {
      title.Append(Form(" part %d", i+1));
      infile.Append(Form("_p%d", i+1));
    }
    infile.Append(".root");
    run[i] = new TSolSimFile(infile,title);
  }
  if( nseg == 1 && nevent > 0 )
    run[0]->SetLastEvent(nevent);

  bool fail = true;
  if( analyzer->Init(run[0]) == 0 ) {
    cout << "initialization successful..." << endl;
    THaDetector* tracker1 = SBS_FPP->GetDetector("tracker.1");
    THaDetector* tracker2 = SBS_FT->GetDetector("tracker.1");
    if( tracker2 ) {
      // The SoLID trackers' origin really is the origin of the first plane
      Double_t z1 = tracker1->GetOrigin().Z();
      Double_t z2 = tracker2->GetOrigin().Z();
      cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
      manager->SetZ0(z2);
      //TSBSSimDecoder::SetZ0(z2);
    } else {
      cerr << "ERROR: cannot get tracker detector! z0 may be wrong" << endl;
    }
    // Emulate dummy calorimeter planes
    // TSBSSimDecoder::SetCaloZ(3.20);  // Calo front is at z = 320 cm
    // TSBSSimDecoder::SetCaloRes(0.01);  // Calo resolution 1 cm (sigma)
    // TSBSSimDecoder::EmulateCalorimeter(false);
    manager->EmulateCalorimeter(false);
    
    
    // Process the runs
    Int_t ret = 0, ntotal = 0;
    for( int i=0; i<nseg && ret >= 0; ++i ) {
      //cout << "processing segment i / nseg : " << i << "/" << nseg << endl;
      if( i>0 )
	run[i]->SetDate(run[0]->GetDate());
      //cout << "Start processing " << endl;
      ret = analyzer->Process(run[i]);
      cout << "ret = " << ret << endl;
      if( ret > 0 )
	ntotal += ret;
    }
    fail = (ret < 0 );
    if( fail )
      cerr << "Terminated with analyzer error = " << ret << endl;
    else
      cout << "Analyzed " << ntotal << " events" << endl;
    analyzer->Close();
  }
  
  for( int i=0; i<nseg; ++i ) {
    delete run[i]; run[i] = 0;
  }
  delete analyzer; analyzer = 0;
  gHaApps->Delete();
  
  //TFile* f =
  if( !fail )
    new TFile(rootfile,"r");
} 
