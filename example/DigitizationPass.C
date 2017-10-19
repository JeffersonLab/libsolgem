// Example "replay" script
//#define DEBUG 1

void DigitizationPass(UInt_t fspec = 1, // Spectrometer flag: 
		      // 1 for BB GEMs(GMn, GEn, SIDIS); 
		      // 3 for Front Tracker spectrometer (GEp);
		      // 4 for Focal Plane Polarimter (GEp).
		      UInt_t NsigFiles = 1,//Number of signal files to analize
		      UInt_t Nmax = -1, //number of events to digitize
		      UInt_t nbacktoadd = 2, // number of background *files* to add to each event
		      bool mips = false,// boolean to digitize mips instead of regular signal.
		      bool print = false){
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");
  
  TDatime run_time = 991231;
  
  gSystem->Load("../libsolgem.so");
  
  ////////////////////////////////////////////////////////////////
  
  int Ngood = 0;
      
  TSBSGEMChamber *ddy;
  TSBSSpec *dds;
  TSBSSimGEMDigitization *ddd;
  TSBSSimDecoder *dde;

  char* outname;
  string bg = "bkgd";
  if(nbacktoadd==0)bg = "nobkgd";

  string infile_sig_prefix;
  string infile_bkgd_prefix;
  
  TSolDBManager* manager = TSolDBManager::GetInstance();
  switch(fspec){
  case(1):
    manager->LoadGeneralInfo("db_generalinfo_bbgem.dat");
    manager->LoadGeoInfo("g4sbs_bbgem");
    dds = new TSBSSpec ("g4sbs_bbgem", "BB spectrometer");
    outname = Form("digitized_bbgem_%s.root", bg.c_str());
    infile_sig_prefix = "/volatile/halla/sbs/efuchey/gmn_elastic/gmn13.5_elastic_sig_20171018_14";
    if(mips)infile_sig_prefix = "/volatile/halla/sbs/efuchey/misc/gmn13.5_BBgemMIP_20171018_14";
    infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gmn_beam_bkgd/gmn13.5_beam_bkgd_20170630_14";
    dds->Init(run_time);
    break;
  case(3):
    manager->LoadGeneralInfo("db_generalinfo_ft.dat");
    manager->LoadGeoInfo("g4sbs_ft");
    dds = new TSBSSpec ("g4sbs_ft", "SBS spectrometer FT");
    outname = Form("digitized_ft_%s.root", bg.c_str());
    infile_sig_prefix = "/volatile/halla/sbs/efuchey/gep_elastic/gep12_elastic_sig_20170727_14";
    if(mips)infile_sig_prefix = "/volatile/halla/sbs/efuchey/gep12_FTFPP_MIP_20170828_10";
    infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gep_beam_bkgd/gep12_beam_bkgd_20170114_11";
    dds->Init(run_time);
    break;
  case(4):
    manager->LoadGeneralInfo("db_generalinfo_fpp.dat");
    manager->LoadGeoInfo("g4sbs_fpp");
    dds = new TSBSSpec ("g4sbs_fpp", "SBS spectrometer FPP");
    outname = Form("digitized_fpp_%s.root", bg.c_str());
    infile_sig_prefix = "/volatile/halla/sbs/efuchey/gep_elastic/gep12_elastic_sig_20170727_14";
    if(mips)infile_sig_prefix = "/volatile/halla/sbs/efuchey/gep12_FTFPP_MIP_20170828_10";
    infile_bkgd_prefix = "/volatile/halla/sbs/efuchey/gep_beam_bkgd/gep12_beam_bkgd_20170114_11";
    dds->Init(run_time);
    break;
  default:
    cout << "No corresponding geometry; choose: " << endl 
	 << "1 (BBGEM)" << endl << "3 (FT)" << endl << "4 (FPP)" << endl;
    return;
    break;
  }
  
  cout << "1  " << outname << " " << &outname << endl;
  
  for(int i_ch = 0; i_ch<manager->GetNChamber()*manager->GetNSector(); i_ch++){
    ddy = new TSBSGEMChamber (Form("gem%d",i_ch),Form("Test chamber %d", i_ch));
    ddy->SetApparatus(dds);
    if( ddy->Init() )
      return;
    dds->AddGEM (ddy);
  }
  printf("\n");

  //cout << outname << " " << &outname << endl;
  
  if(print)dds->Print();
    
  ddd = new TSBSSimGEMDigitization (*dds,"ratedig");
  //ddd = new TSBSSimGEMDigitization (*dds,"testdigitization");
  ddd->SetMapSector(false);
    
  ////////////////////////////////////////////////////////////////
  
  int nevent = 0;

  TSBSGeant4File *f;
  
  int  ndata, i;
  TSolGEMData *gd, *gb;
  g4sbsgendata *gen;
  
  cout << "creating file " << outname << endl;
  ddd->InitTree (*dds, outname);
    
  printf("Digitizing events\n");
  ndata = 0;
  
  //Nmax = TMath::Min((Long64_t)Nmax, f->GetEntries());
    
  int hadback = 1;
  
  int N_bg_file_g = 0;
  
  //Add the loop on the signal files
  for(int i_sig = 0; i_sig<NsigFiles; i_sig++){
    if(mips){
      f = new TSBSGeant4File(Form("%s/gun_%d.root",infile_sig_prefix.c_str(), i_sig));
    }else{
      f = new TSBSGeant4File(Form("%s/elastic_%d.root",infile_sig_prefix.c_str(), i_sig));
    }
    printf("The filename returned is %s\n", f->GetFileName());
    f->SetSource(0);
    
    int res;
    
    res = f->Open();
    
    if( res != 1 ){
      printf("Opening g4sbs file returned %d\n", res);
      continue;
    }
    
    ////////////////////////////////////////////////////////////////
    int d_flag_readevent = 0;
    while( f->ReadNextEvent(d_flag_readevent) && hadback && nevent<Nmax ){
      
      if(nevent%100==0){
	cout << "Evt " << nevent << endl;
      }
      
      if(f->GetNData()==0){
	if(print)cout << "No hits, skip evt " << nevent << endl;
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
      
      ddd->SetTreeEvent((*gd), (*f), nevent);
      
      /*
	if(print){// || nevent==40
	cout << "number of hits in GEM data " << gd->GetNHit() << endl;
	while(ndata<gd->GetNHit()){
	
	//if(gd->GetParticleID(ndata)>1)continue;
	gd->Print();
	cout << "hit number " << ndata << endl;
	gd->PrintHit(ndata);
	ndata++;
	}
	}
      */
      ddd->Digitize(*gd, *dds);
      
      // Access to generated vertex and momentum
      // gen->GetV();
      // gen->GetP();
      // gen->GetWeight();
      
      // Add some number of background files...
      int N_bg_file_g_post = N_bg_file_g+nbacktoadd;
      
      if(nbacktoadd){
	for(int Nfile = N_bg_file_g; Nfile < N_bg_file_g_post; Nfile++){
	  //if(print)cout << N_bg_file_g << " <= " << Nfile << " < " << N_bg_file_g+nbacktoadd << endl;
	  TSBSGeant4File *fback = new TSBSGeant4File(Form("%s/beam_bkgd_%d.root",infile_bkgd_prefix.c_str(), Nfile));
	  int open = fback->Open();
	  if(!open){
	    N_bg_file_g_post++;
	    N_bg_file_g++;
	    
	    if(N_bg_file_g>=2000){
	      int n_temp = Nfile;
	      Nfile = N_bg_file_g_post-n_temp;
	      N_bg_file_g_post = nbacktoadd;
	      N_bg_file_g = 0;
	    }
	    //if(print)cout << Form("/group/exjpsi/eric/31722/beam_bkgd_%d.root does not exist", Nfile) << endl;
	  continue;
	  }
	  
	  fback->SetSource(1);
	  
	  int backidx = 0;
	  //while( hadback = fback->ReadNextEvent() && backidx < nbacktoadd ){
	  while( backidx < fback->GetEntries() ){
	    hadback = fback->ReadNextEvent();
	    gb = fback->GetGEMData();
	    
	    if(print && gb->GetNHit()>0){
	      cout << "Bkgd evt: " << gb->GetEvent() << ", number of hits " 
		   << gb->GetNHit() << endl;
	      int nback = 0;
	      while(nback<gb->GetNHit()){
		//if(gd->GetParticleID(ndata)>1)continue;
		gb->Print();
		cout << "hit number " << nback << endl;
		gb->PrintHit(nback);
		nback++;
	      }
	      
	    }
	    
	    ddd->AdditiveDigitize(*gb, *dds);
	    
	    // //Randomize times based on gate width
	    // for( int bidx = 0; bidx < gb->GetNHit(); bidx++ ){
	    //   double timeshift = gRandom->Uniform(-ddd->GetGateWidth(), 75.0 );//ns
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
	//if(print)cout << "new number of hits in GEM data " << gd->GetNHit() << endl;
	N_bg_file_g = N_bg_file_g_post;
      }//end if nbacktoadd
      
      if(N_bg_file_g>=2000)N_bg_file_g = 0;
      
      ddd->FillTree();
      
      //if(nevent==7)ddd->GetEvent()->Print("all");
    if(print)ddd->GetEvent()->Print("all");
    //ddd->GetEvent()->Print("clust");
    
    delete gd;
    nevent++;
    }// end read flag
    cout << "closing file " << endl;
    f->Close();
  }
  printf("Completed %d events total: %d good events \n", nevent, Ngood);
  ddd->WriteTree();
  ddd->CloseTree();
  
  cout << "Tree closed" << endl;
}
