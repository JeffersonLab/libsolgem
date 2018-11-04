#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TDatime.h"
#include "TSBSGeant4File.h"
#include "TSBSDBManager.h"
#include "TSBSGEMChamber.h"
#include "TSBSSpec.h"
#include "TSBSSimGEMDigitization.h"
#include "TSBSSimDecoder.h"
#include "THaAnalysisObject.h"
#endif
// Example "replay" script
//#define DEBUG 1
using namespace std;

void DigitizationPass(UInt_t fspec = 1, // Spectrometer flag: 
		      // 1 for BB GEMs(GMn, GEn, SIDIS); 
		      // 3 for Front Tracker spectrometer (GEp);
		      // 4 for Focal Plane Polarimter (GEp).
		      UInt_t NsigFiles = 1,//Number of signal files to analize
		      string infile_sig_prefix = "/work/halla/sbs/efuchey/gmn13.5_elastic_sig_20180709_22",
		      UInt_t Nmin = 0, //first event to analyze
		      UInt_t Nmax = -1, //number of events to digitize
		      UInt_t nbacktoadd = 0, // number of background *files* to add to each event
		      string infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gmn13.5_beam_bkgd_20180718_22",
		      UInt_t firstbackfilenum = 0
		      ){
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");
  
  TDatime run_time = 991231;
  
  gSystem->AddDynamicPath("${LIBSBSGEM}");
  gSystem->Load("../libsolgem.so");
  
  ////////////////////////////////////////////////////////////////
  
  int Ngood = 0;
  
  TSBSGEMChamber *dGEM;
  TSBSSpec *dSpectro;
  TSBSSimGEMDigitization *dSimDigi;
  
  string outname;
  string bg = "bkgd";
  if(nbacktoadd==0)bg = "nobkgd";

  // string infile_sig_prefix;
  // string infile_bkgd_prefix;
  
  TSBSDBManager* manager = TSBSDBManager::GetInstance();
  switch(fspec){
  case(1):
    manager->LoadGeneralInfo("db_generalinfo_bbgem.dat");
    manager->LoadGeoInfo("g4sbs_bbgem"); //new plane--module Geo for SBS GEMs
    dSpectro = new TSBSSpec ("g4sbs_bbgem", "BB spectrometer");
    outname = Form("digitized_bbgem_%s.root", bg.c_str());
    dSpectro->Init(run_time);
    break;
  case(3):
    manager->LoadGeneralInfo("db_generalinfo_ft.dat");
    manager->LoadGeoInfo("g4sbs_ft");
    dSpectro = new TSBSSpec ("g4sbs_ft", "SBS spectrometer FT");
    outname = Form("digitized_ft_%s.root", bg.c_str());
    dSpectro->Init(run_time);
    break;
  case(4):
    manager->LoadGeneralInfo("db_generalinfo_fpp.dat");
    manager->LoadGeoInfo("g4sbs_fpp");
    dSpectro = new TSBSSpec ("g4sbs_fpp", "SBS spectrometer FPP");
    outname = Form("digitized_fpp_%s.root", bg.c_str());
    dSpectro->Init(run_time);
    break;
  default:
    cout << "No corresponding geometry; choose: " << endl 
	 << "1 (BBGEM)" << endl << "3 (FT)" << endl << "4 (FPP)" << endl;
    return;
    break;
  }
  
  cout << "outputfile:  " << outname << endl;
  
  //Init TSBSGEMChamber and add to TSBSSpec
  for(int i_plane=0;i_plane<manager->GetNGEMPlane();i_plane++)
    {
      for(int i_module=0;i_module<manager->GetNModule(i_plane);i_module++)
	{
	  dGEM = new TSBSGEMChamber (Form("plane%d.module%d",i_plane,i_module),Form("Test chamber on Plane: %d, Module: %d", i_plane, i_module));
	  dGEM->SetApparatus(dSpectro);
	  if( dGEM->Init() )//Loading GEM Geo information
	    return;
	  dSpectro->AddGEM (dGEM);
	}
    }
  
  dSimDigi = new TSBSSimGEMDigitization (*dSpectro,"ratedig");
 
  ////////////////////////////////////////////////////////////////
    
  UInt_t nevent = Nmin;
  
  TSBSGeant4File *f;

  TSBSGEMSimHitData *gd, *gb;
  g4sbsgendata *gen;
  
  cout << "creating file " << outname <<&outname<<endl;
  dSimDigi->InitTree (*dSpectro, outname);
    
  printf("Digitizing events\n");

  int hadback = 1;
  
  int N_bg_file_g = firstbackfilenum;
  
  cout << "start loop on sig files " << endl;
  
  //Add the loop on the signal files
  for(int i_sig = 0; i_sig<NsigFiles; i_sig++){
    f = new TSBSGeant4File(Form("%s/elastic_%d.root",infile_sig_prefix.c_str(), i_sig));
    printf("The filename returned is %s\n", f->GetFileName());
    f->SetSource(0);
    int res;
    
    res = f->Open();
    
    if( res != 1 ){
      printf("Opening g4sbs file returned %d\n", res);
      continue;
    }
    
    f->SetFirstEvNum(Nmin);
    cout << "about to go through events for file " << f->GetFileName() << endl;
    if(f->GetEvNum()!=Nmin-1){
      cout << "f->GetEvNum() = " << f->GetEvNum() << " != Nmin-1 =" << Nmin-1 << endl;
      exit(-1);
    }
    int d_flag_readevent = 0;
    while( f->ReadNextEvent(d_flag_readevent) && hadback && nevent<Nmax ){

      if(nevent%100==0 || (nbacktoadd>0 && nevent%10==0) ){
      cout << "Evt " << nevent << endl;
      }
    
      if(f->GetNData()==0){
	nevent++;
	continue;
      }    
      gd = f->GetGEMData();

      if(f->GetNGen()>0){
	gen = f->GetGenData(0);
	Ngood++;
      }else{
	cout << "No generated data for event " << nevent 
	     << ", skip it (Nhits = " << f->GetNData() << ")" << endl;
	nevent++;
	continue;
      }
      
      dSimDigi->SetTreeEvent((*gd), (*f), nevent);
    
      dSimDigi->Digitize(*gd, *dSpectro);
      
      // Access to generated vertex and momentum
      // gen->GetV();
      // gen->GetP();
      // gen->GetWeight();
      
      // Add some number of background files...
      int N_bg_file_g_post = N_bg_file_g+nbacktoadd;
      
      if(nbacktoadd){
	for(int Nfile = N_bg_file_g; Nfile < N_bg_file_g_post; Nfile++){
	  TSBSGeant4File *fback = new TSBSGeant4File(Form("%s/beam_bkgd_%d.root",infile_bkgd_prefix.c_str(), Nfile));
	  int open = fback->Open();
	  if(!open){
	    N_bg_file_g_post++;
	    N_bg_file_g++;
	    
	    if(N_bg_file_g>=10000){
	      int n_temp = Nfile;
	      Nfile = N_bg_file_g_post-n_temp;
	      N_bg_file_g_post = nbacktoadd;
	      N_bg_file_g = 0;
	    }
	    continue;
	  }
	  
	  fback->SetSource(1);
	  
	  int backidx = 0;
	  //while( hadback = fback->ReadNextEvent() && backidx < nbacktoadd ){
	  while( backidx < fback->GetEntries() ){
	    hadback = fback->ReadNextEvent();
	    gb = fback->GetGEMData();
	    
	    dSimDigi->AdditiveDigitize(*gb, *dSpectro);
	    
	    // //Randomize times based on gate width
	    // for( int bidx = 0; bidx < gb->GetNHit(); bidx++ ){
	    //   double timeshift = gRandom->Uniform(-dSimDigi->GetGateWidth(), 75.0 );//ns
	    //   gb->SetHitTime(bidx, gb->GetHitTime(bidx) + timeshift );
	    // }	
	    // //gd->AddGEMData(gb);
	    backidx++;
	  }
	  
	  // if( backidx != nbacktoadd ){
	  // printf("Warning:  Not enough background events to be added (%d)\n", backidx);
	  // }
	  
	  fback->Close();
	}
	N_bg_file_g = N_bg_file_g_post;
      }//end if nbacktoadd
    
      if(N_bg_file_g>=10000)N_bg_file_g = 0;
      
      dSimDigi->FillTree();
      
      delete gd;
      if(nbacktoadd)delete gb;//Also needs to delete gb...Or there will be huge memory leak, especially at high background
      nevent++;
    }
    f->Close();
  }
  printf("Completed %d events total: %d good events \n", nevent, Ngood);

  dSimDigi->WriteTree();
  dSimDigi->CloseTree();
  
  cout << "Tree closed" << endl;
}

int main(int argc, char *argv[])
{
  UInt_t fspec = 1; // Spectrometer flag: 
  // 1 for BB GEMs(GMn, GEn, SIDIS); 
  // 3 for Front Tracker spectrometer (GEp);
  // 4 for Focal Plane Polarimter (GEp).
  UInt_t NsigFiles = 1;//Number of signal files to analize
  string infile_sig_prefix = "/work/halla/sbs/efuchey/gmn13.5_elastic_sig_20180709_22";
  UInt_t Nmin = 0; //first event to analyze 
  UInt_t Nmax = -1; //number of events to digitize
  UInt_t nbacktoadd = 0; // number of background *files* to add to each event
  string infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gmn13.5_beam_bkgd_20180718_22";
  UInt_t firstbackfilenum = 0;
  
  if(argc>9) {
    cout << "Usage: DigiPass" << endl 
	 << " <fspec (1: BB, 3: FT, 4: FPP)>: spec flag " << endl 
	 << " <NsigFiles>: number of signal files to select" << endl
	 << " <infile_sig_prefix>: path to fetch signal files" << endl 
	 << " <Nmin>: first signal event to analyze " << endl
	 << " <Nmax>: maximum event number to analyze" << endl 
	 << " <nbacktoadd>: number of background files to add to signal" << endl
	 << " <infile_sig_prefix>: path to fetch background files" << endl;
    return -1;
  }
  cout << "run digi pass with following parameters:" << endl;
  if(argc>1)fspec = atoi(argv[1]);
  cout << " fspec = " << fspec << endl;
  if(argc>2)NsigFiles = atoi(argv[2]);
  cout << " NsigFiles = " << NsigFiles << endl;
  if(argc>3)infile_sig_prefix = argv[3];
  cout << " infile_sig_prefix = " << infile_sig_prefix.c_str() << endl;
  if(argc>4)Nmin = atoi(argv[4]);
  cout << " Nmin " << Nmin << endl;
  if(argc>5)Nmax = atoi(argv[5]);
  cout << " Nmax = " << Nmax << endl;
  if(argc>6)nbacktoadd = atoi(argv[6]);
  cout << " nbacktoadd = " << nbacktoadd << endl;
  if(argc>7)infile_bkgd_prefix = argv[7];
  cout << " infile_sig_prefix = " << infile_bkgd_prefix.c_str() << endl;
  if(argc>8)firstbackfilenum = atoi(argv[8]);
  cout << " firstbackfilenum = " << firstbackfilenum << endl;
 
  DigitizationPass(fspec, NsigFiles, infile_sig_prefix, Nmin, Nmax, nbacktoadd, infile_bkgd_prefix);
  return 0;
};
