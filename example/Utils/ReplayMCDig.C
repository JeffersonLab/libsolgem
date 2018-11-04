#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TDatime.h"
#include "THaGlobals.h"
#include "THaInterface.h"
#include "THaTextvars.h"
#include "THaVarList.h"
#include "THaCutList.h"
#include "THaRunBase.h"
#include "THaAnalyzer.h"
#include "CodaDecoder.h"
#include "THaAnalysisObject.h"
#include "TSBSGeant4File.h"
#include "TSBSSimFile.h"
#include "TSBSDBManager.h"
#include "TSBSGEMChamber.h"
#include "TSBSGEMPlane.h"
#include "TSBSSpec.h"
#include "TSBSSimGEMDigitization.h"
#include "TSBSSimDecoder.h"
#include "THaInterface.h"
#include "SBSSpec.h"
#endif
/// #ifndef __CINT__

//#include "SBSSpec.h"
//#include "TSBSSimDecoder.h"
//#include "TSBSSimFile.h"

// #include "THaInterface.h"
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

void ReplayMCDigitized(string filename = "digitized", 
		       string detsuffix = "bbgem",//detector suffix: 
		       //"bbgem" for BigBite spectrometer (GMn, GEn, SIDIS);
		       //"FT" for Front Tracker spectrometer (GEp);
		       //"FPP" for Focal Plane Polarimeters (GEp).
		       bool bkgd = false,// flag to indicate if digitized file includes background or not.
		       UInt_t nevent = -1, // number of events to process
		       UInt_t nseg = 0, // number of segments
		       bool do_cuts = true )
{
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");
  //gSystem->Getenv("ANALYZER");
  
  if( nseg > 100 ){
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
  manager->LoadGeneralInfo(Form("db_generalinfo_%s.dat", detsuffix.c_str()));
  manager->LoadGeoInfo(Form("g4sbs_%s", detsuffix.c_str()));
  
  int Ns = manager->GetNSector();
  const int Nsect = Ns;

  cout << "detector " << detsuffix.c_str() << ": " << Nsect << " sectors." << endl;

  // dde = new TSBSSimDecoder();
  // dde->SetCrateMapName("db_sbssim_cratemap.dat");

  THaInterface::SetDecoder( TSBSSimDecoder::Class() );
  
  cout << "Reading " << detsuffix.c_str() << endl;
  THaApparatus* SBS_GEMdet = new SBS::SBSSpec( Form("sbs_%s",detsuffix.c_str()), Form("SBS / %s ", detsuffix.c_str()), Nsect );
  //THaApparatus* SBS_GEMdet = new SBSSpec( Form("sbs_%s",detsuffix.c_str()), "SBS / BB GEMs", Nsect );
  //TList* gHaApps = new TList();
  cout << "muh ? " << endl;
  gHaApps->Add( SBS_GEMdet );
  cout << "Just read " << detsuffix.c_str() << endl;

  SBS_GEMdet->Print("DET");

  TString db_prefix = SBS_GEMdet->GetName();
  db_prefix += ".tracker";
  gHaTextvars->Add( "DET", db_prefix.Data() );
  gHaTextvars->Add( "APP", SBS_GEMdet->GetName() );

  THaAnalyzer* analyzer = new THaAnalyzer;

  TString rootfile(Form("%s_%s_%s", filename.c_str(), detsuffix.c_str(), bg.c_str())), infile0(Form("%s_%s_%s", filename.c_str(), detsuffix.c_str(), bg.c_str()));
  TString odeffile("sbssim.odef"), cutfile(Form("sbs_%ssim.cuts",detsuffix.c_str()));
  rootfile.Append("_replayed_new.root");
  analyzer->EnableBenchmarks();
  analyzer->SetOutFile(rootfile);
  analyzer->SetOdefFile(odeffile);
  if( do_cuts ) analyzer->SetCutFile(cutfile);
  analyzer->SetSummaryFile(Form("%s_%s_new.sum", filename.c_str(), detsuffix.c_str()));
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
  //}

  //TFile* f =
  if( !fail )
    new TFile(rootfile,"r");
} 

int main(int argc, char *argv[])
{
  gHaVars    = new THaVarList;
  gHaCuts    = new THaCutList( gHaVars );
  gHaApps    = new TList;
  gHaPhysics = new TList;
  gHaEvtHandlers = new TList;
  // Use the standard CODA file decoder by default
  gHaDecoder = Decoder::CodaDecoder::Class();
  // File-based database by default
  // gHaDB      = new THaFileDB();
  gHaTextvars = new THaTextvars;
  
  string filename = "digitized";
  string detsuffix = "bbgem";
  bool bkgd = true;
  UInt_t nevent = -1;
  UInt_t nseg = 0;
  bool do_cuts = true;

  if(argc>7) {
    cout << "Usage: ReplayMCDig" << endl 
	 << " <filename>: file name prefix " << endl 
	 << " <detsuffix>: detector suffix (usually bbgem, ft, fpp)" << endl
	 << " <bkgd>: flag to analyze background (0: no; 1: yes)" << endl 
	 << " <nevent>: number of events to analyze " << endl
	 << " <nseg>: number of segments to analyze" << endl 
	 << " <do_cuts>: perform cuts" << endl;
    return -1;
  }
  
  cout << "run replay mc dig with following parameters:" << endl;
  if(argc>1)filename = argv[1];
  cout << " filename = " << filename.c_str() << endl;
  if(argc>2)detsuffix = argv[2];
  cout << " detsuffix = " << detsuffix.c_str() << endl;
  if(argc>3)bkgd = (bool)atoi(argv[3]);
  cout << " bkgd = " << bkgd << endl;
  if(argc>4)nevent = atoi(argv[4]);
  cout << " nevent = " << nevent << endl;
  if(argc>5)nseg = atoi(argv[5]);
  cout << " nseg = " << nseg << endl;
  if(argc>6)do_cuts = (bool)atoi(argv[6]);
  cout << " do_cuts  = " << do_cuts << endl;
  
  ReplayMCDigitized(filename, detsuffix, bkgd, nevent, nseg, do_cuts);
  return 0;  
};
