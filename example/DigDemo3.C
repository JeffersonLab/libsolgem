// Example "replay" script

void DigDemo3(int Nmax = 10000, bool print = false){
    printf("\n** This gets called with 'analyzer' and not 'root' **\n");
    printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

    gSystem->Load("../libsolgem.so");
    
    ////////////////////////////////////////////////////////////////

    TSBSGEMChamber *ddy;
    TSBSSpec *dds;
    TSBSSimGEMDigitization *ddd;
    //TSBSSimDecoder *dde;
    TSolSimDecoder *dde;
    
    dds = new TSBSSpec ("g4sbs", "SBS spectrometer");
    dds->Init();
    
    for(int i_ch = 0; i_ch<10; i_ch++){
      ddy = new TSBSGEMChamber (Form("gem%d",i_ch),Form("Test chamber %d", i_ch));
      ddy->SetApparatus(dds);
      if( ddy->Init() )
	return;
      dds->AddGEM (ddy);
    }
    printf("\n");

    for(int i_ch = 10; i_ch<16; i_ch++){
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
    
    //dde = new TSBSSimDecoder ();
    dde = new TSolSimDecoder ();
    dde->SetCrateMapName("db_solsim_cratemap.dat");
    
    ////////////////////////////////////////////////////////////////
    
    TSBSGeant4File *f = new TSBSGeant4File("/group/exjpsi/eric/14001/elastic_0.root");
    //TSBSGeant4File *f = new TSBSGeant4File("/group/exjpsi/eric/14001/elastic_0.root");
    printf("The filename returned is %s\n", f->GetFileName());
    f->SetSource(0);
    
    int res;
    
    res = f->Open();

    if( res != 1 ){
	printf("Opening g4sbs file returned %s\n", res);
	return;
    }

    // Hypothetical background run
    TSBSGeant4File *fback = new TSBSGeant4File("/group/exjpsi/eric/11001/beam_bkgd_0.root");
    //TSBSGeant4File *fback = new TSBSGeant4File("/group/exjpsi/eric/11001/beam_bkgd_0.root");
    fback->Open();
    fback->SetSource(1);
 
    ////////////////////////////////////////////////////////////////
    
    int nevent = 0;
    
    int  ndata, i;
    TSolGEMData *gd, *gb;
    g4sbsgendata *gen;

    ddd->InitTree (*dds, "digdemo3.root");
    
    printf("Digitizing events\n");
    ndata = 0;
        
    int hadback = 1;

    while( f->ReadNextEvent() && hadback && nevent<Nmax){
      
      if(nevent%100==0)cout << "Evt " << nevent << endl;
            
      gd = f->GetGEMData();
      gen = f->GetGenData(0);
      
      ddd->SetTreeEvent((*gd), (*f), nevent);
      
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
      
      ddd->Digitize(*gd, *dds);
      
      // Access to generated vertex and momentum
      // gen->GetV();
      // gen->GetP();
      // gen->GetWeight();
      
      // Add some number of background events
      int nbacktoadd = 0;
      int backidx = 0;
      
      //while( hadback = fback->ReadNextEvent() && backidx < nbacktoadd ){
      while( backidx < nbacktoadd ){
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
	
	// // Randomize times based on gate width
	//  for( int bidx = 0; bidx < gb->GetNHit(); bidx++ ){
	//    double timeshift = gRandom->Uniform(-ddd->GetGateWidth(), 75.0 );//ns
	//    gb->SetHitTime(bidx, gb->GetHitTime(bidx) + timeshift );
	//  }	
	//  //gd->AddGEMData(gb);
	backidx++;
      }
      
      if(print)cout << "new number of hits in GEM data " << gd->GetNHit() << endl;
      
      if( backidx != nbacktoadd ){
        printf("Warning:  Not enough background events to be added (%d)\n", backidx);
      }
      
      //if(nevent==25)
      //ddd->GetEvent()->Print("");
      
      ddd->FillTree();
      
      //if(nevent==25)ddd->GetEvent()->Print("all");
      
      delete gd;
      nevent++;
    }
    printf("Completed %d events\n", nevent);

    ddd->WriteTree();
    ddd->CloseTree();

}
