// Example "replay" script

void DigDemo3(int fspec = 1, int Nmax = 1000, int nbacktoadd = 2, bool print = false){
  printf("\n** This gets called with 'analyzer' and not 'root' **\n");
  printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

  gSystem->Load("../libsolgem.so");
    
  ////////////////////////////////////////////////////////////////
    
  int Ngood = 0;
    
  TSBSGEMChamber *ddy;
  TSBSSpec *dds;
  TSBSSimGEMDigitization *ddd;
  TSBSSimDecoder *dde;
    
  char* outname;
    
  TSolDBManager* manager = TSolDBManager::GetInstance();
  switch(fspec){
  case(1):
    manager->LoadGeneralInfo("db_generalinfo_bbgem.dat");
    manager->LoadGeoInfo("g4sbs_bbgem");
    dds = new TSBSSpec ("g4sbs_bbgem", "BB spectrometer");
    outname = "digdemo3_bbgem.root";
    dds->Init();
    break;
  case(3):
    manager->LoadGeneralInfo("db_generalinfo_ft.dat");
    manager->LoadGeoInfo("g4sbs_ft");
    dds = new TSBSSpec ("g4sbs_ft", "SBS spectrometer FT");
    outname = "digdemo3_ft.root";
    dds->Init();
    break;
  case(4):
    manager->LoadGeneralInfo("db_generalinfo_fpp.dat");
    manager->LoadGeoInfo("g4sbs_fpp");
    dds = new TSBSSpec ("g4sbs_fpp", "SBS spectrometer FPP");
    outname = "digdemo3_fpp.root";
    dds->Init();
    break;
  default:
    cout << "No corresponding geometry; choose: " << endl 
	 << "1 (BBGEM)" << endl << "3 (FT)" << endl << "4 (FPP)" << endl;
    return;
    break;
  }
    
  for(int i_ch = 0; i_ch<manager->GetNChamber()*manager->GetNSector(); i_ch++){
    ddy = new TSBSGEMChamber (Form("gem%d",i_ch),Form("Test chamber %d", i_ch));
    ddy->SetApparatus(dds);
    if( ddy->Init() )
      return;
    dds->AddGEM (ddy);
  }
  printf("\n");
    
  if(print)dds->Print();
    
  ddd = new TSBSSimGEMDigitization (*dds,"ratedig");
  //ddd = new TSBSSimGEMDigitization (*dds,"testdigitization");
  ddd->SetMapSector(false);
    
  ////////////////////////////////////////////////////////////////
    
  TSBSGeant4File *f = new TSBSGeant4File("/work/halla/sbs/efuchey/gmn_elastic/gmn13.5_elastic_sig_20170504_10/elastic_0.root");
  printf("The filename returned is %s\n", f->GetFileName());
  f->SetSource(0);
    
  int res;
    
  res = f->Open();

  if( res != 1 ){
    printf("Opening g4sbs file returned %s\n", res);
    return;
  }

  ////////////////////////////////////////////////////////////////
    
  int nevent = 0;
    
  int  ndata, i;
  TSolGEMData *gd, *gb;
  g4sbsgendata *gen;
    
  ddd->InitTree (*dds, outname);
    
  printf("Digitizing events\n");
  ndata = 0;
  Nmax = TMath::Min((Long64_t)Nmax, f->GetEntries());
    
  int hadback = 1;
  
  int N_bg_file_g = 0;
  
  while( f->ReadNextEvent() && hadback && nevent<Nmax ){

    if(nevent%100==0){
      cout << "Evt " << nevent << endl;
    }
    
    if(f->GetNData()==0){
      //if(print)cout << "No hits, skip evt " << nevent << endl;
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
    
    ddd->Digitize(*gd, *dds);
      
    // Access to generated vertex and momentum
    // gen->GetV();
    // gen->GetP();
    // gen->GetWeight();
      
    // Add some number of background files...
    int N_bg_file_g_post = N_bg_file_g+nbacktoadd;
 
    if(nbacktoadd){
      for(int Nfile = N_bg_file_g; Nfile < N_bg_file_g_post; Nfile++){
	TSBSGeant4File *fback = new TSBSGeant4File(Form("/work/halla/sbs/efuchey/gmn_beam_bkgd/gmn13.5_beam_bkgd_20170512_18/beam_bkgd_%d.root", Nfile));
	int open = fback->Open();
	if(!open){
	  N_bg_file_g_post++;
	  N_bg_file_g++;
	  
	  if(N_bg_file_g>=500){
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
    
    if(N_bg_file_g>=500)N_bg_file_g = 0;
    
    ddd->FillTree();
      
    //if(nevent==25)ddd->GetEvent()->Print("all");
    //if(print)ddd->GetEvent()->Print("all");
      
    delete gd;
    nevent++;
  }
  printf("Completed %d events total: %d good events \n", nevent, Ngood);

  ddd->WriteTree();
  ddd->CloseTree();
        
}
