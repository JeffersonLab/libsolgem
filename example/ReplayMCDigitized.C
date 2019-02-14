#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSSimFile.h"
#include "TSBSDBManager.h"
#include "TSBSGEMChamber.h"
#include "TSBSGEMPlane.h"
#include "TSBSSpec.h"
#include "TSBSSimGEMDigitization.h"
#include "TSBSSimDecoder.h"
#include "THaInterface.h"
#include "THaTextvars.h"
#include "THaRunBase.h"
#include "THaAnalyzer.h"
#include "THaAnalysisObject.h"
#include "SBSSpec.h"
#endif
/// #ifndef __CINT__

//#include "SBSSpec.h"
//#include "TSBSSimDecoder.h"
//#include "TSBSSimFile.h"

// #include "THaInterface.h"
// #include "THaGlobals.h"
// #include "THaTextvars.h"
// #include "THaAnalyzer.h"
// #include "THaDetector.h"

// #include "TSystem.h"
// #include "TList.h"
// #include "TString.h"
// #include "TFile.h"
// #include "TVector3.h"

// #include <iostream>

// #endif
using namespace std;

void ReplayMCDigitized(const char* filename = "digitized", 
		       const char* detsuffix = "ft",//detector suffix: 
		       //"bbgem" for BigBite spectrometer (GMn, GEn, SIDIS);
		       //"FT" for Front Tracker spectrometer (GEp);
		       //"FPP" for Focal Plane Polarimeters (GEp).
		       bool bkgd = false,// flag to indicate if digitized file includes background or not.
		       Int_t nevent = -1, // number of events to process
		       Int_t nseg = 0, // number of segments
		       bool do_cuts = true )
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
  
  string bg = "nobkgd";
  if(bkgd)bg = "bkgd";
  
  gSystem->AddDynamicPath("${LIBSBSGEM}");
  gSystem->AddDynamicPath("${TREESEARCH}");
  gSystem->Load("../libsolgem.so");
  gSystem->Load("/home/efuchey/g4work/Tracking/TreeSearch/libTreeSearch.so");
  gSystem->Load("/home/efuchey/g4work/Tracking/TreeSearch/libTreeSearch-SBS.so");
  gSystem->Load("libMinuit");

  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  manager->LoadGeneralInfo(Form("db_generalinfo_%s.dat", detsuffix));
  manager->LoadGeoInfo(Form("g4sbs_%s", detsuffix));

  int Ns = manager->GetNSector();
  const int Nsect = Ns;

  cout << "detector " << detsuffix << ": " << Nsect << " sectors." << endl;

  // dde = new TSBSSimDecoder();
  // dde->SetCrateMapName("db_sbssim_cratemap.dat");

  THaInterface::SetDecoder( TSBSSimDecoder::Class() );
  
  cout << "Reading " << detsuffix << endl;
  THaApparatus* SBS_GEMdet = new SBS::SBSSpec( Form("sbs_%s",detsuffix), "SBS / FT", Nsect );
  gHaApps->Add( SBS_GEMdet );
  cout << "Just read " << detsuffix << endl;

  SBS_GEMdet->Print("DET");

  TString db_prefix = SBS_GEMdet->GetName();
  db_prefix += ".tracker";
  gHaTextvars->Add( "DET", db_prefix.Data() );
  gHaTextvars->Add( "APP", SBS_GEMdet->GetName() );

  THaAnalyzer* analyzer = new THaAnalyzer;

  TString rootfile(Form("%s_%s_%s", filename, detsuffix, bg.c_str())), infile0(Form("%s_%s_%s", filename, detsuffix, bg.c_str()));
  TString odeffile("sbssim.odef"), cutfile(Form("sbs_%ssim.cuts",detsuffix));
  rootfile.Append("_replayed_new.root");
  analyzer->EnableBenchmarks();
  analyzer->SetOutFile(rootfile);
  analyzer->SetOdefFile(odeffile);
  if( do_cuts ) analyzer->SetCutFile(cutfile);
  analyzer->SetSummaryFile(Form("%s_%s_new.sum", filename, detsuffix));
  analyzer->SetCrateMapFileName("sbssim_cratemap");

  //static int Nrun = TMath::Max(nseg,1);
  THaRunBase* run[0];
  TString title0 = "Digitized MC data";
  for( int i=0; i<nseg; ++i ) {
    TString title(title0), infile("rootfiles/"+infile0);
    if( do_parts ) {
      title.Append(Form(" part %d", i+1));
      infile.Append(Form("_p%d", i+1));
    }
    infile.Append(".root");
    run[i] = new TSBSSimFile(infile,title);
  }
  if( nseg == 1 && nevent > 0 )
    run[0]->SetLastEvent(nevent);
  bool fail = true;
  if( analyzer->Init(run[0]) == 0 ) {
    cout << "initialization successful..." << endl;
    THaDetector* tracker[Nsect];
    
    for(int ns = 0; ns < Nsect; ns++){
      //THaDetector* 
      cout << ns << " " << tracker[ns] << endl;
      tracker[ns] = SBS_GEMdet->GetDetector(Form("tracker.%d", ns+1));
      tracker[ns]->Print("");
      if( tracker[ns] ) {
	// The SBS trackers' origin really is the origin of the first plane
	Double_t z0 = tracker[ns]->GetOrigin().Z();
	cout << "z0 = " << z0 << endl;
	manager->SetZ0(z0);
      } else {
	cerr << "ERROR: cannot get tracker detector! z0 may be wrong" << endl;
      }
      // Emulate dummy calorimeter planes
      // TSBSSimDecoder::SetCaloZ(3.20);  // Calo front is at z = 320 cm
      // TSBSSimDecoder::SetCaloRes(0.01);  // Calo resolution 1 cm (sigma)
      // TSBSSimDecoder::EmulateCalorimeter(false);
    }//end loop on trackers
    /*
      THaDetector* tracker = SBS_GEMdet->GetDetector("tracker.1");
      tracker->Print("");
      if( tracker ) {
      // The SBS trackers' origin really is the origin of the first plane
      Double_t z0 = tracker->GetOrigin().Z();
      cout << "z0 = " << z0 << endl;
      manager->SetZ0(z0);
      } else {
      cerr << "ERROR: cannot get tracker detector! z0 may be wrong" << endl;
      }
    */
    //manager->EmulateCalorimeter(false);
    
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
  //}

  //TFile* f =
  if( !fail )
    new TFile(rootfile,"r");
} 

int main()
{
  const char* filename = "digitized";
  const char* detsuffix = "ft";
  bool bkgd = false;
  Int_t nevent = -1;
  Int_t nseg = 0;
  bool do_cuts = true;
  
  ReplayMCDigitized(filename, detsuffix, bkgd, nevent, nseg, do_cuts);
  return 0;
};
