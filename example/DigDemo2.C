// Example "replay" script

void DigDemo2(){
    printf("\n** This gets called with 'analyzer' and not 'root' **\n");
    printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

    gSystem->Load("../libsolgem.so");

    ////////////////////////////////////////////////////////////////


    TSolGEMChamber *ddy;
    TSolSpec *dds;
    TSolSimGEMDigitization *ddd;

    dds = new TSolSpec ("spectrometer", "SOLiD spectrometer");
    dds->Init();

    ddy = new TSolGEMChamber ("gemccham1","Test chamber");
    ddy->Init();
    dds->AddGEM (*ddy);

    printf("\n");

    ddy = new TSolGEMChamber ("gemccham2","Test chamber");
    ddy->Init();
    dds->AddGEM (*ddy);

    printf("\n");

    ddy = new TSolGEMChamber ("gemccham3","Test chamber");
    ddy->Init();
    dds->AddGEM (*ddy);

    printf("\n");

    ddy = new TSolGEMChamber ("gemccham4","Test chamber");
    ddy->Init();
    dds->AddGEM (*ddy);

    dds->Print();

    ddd = new TSolSimGEMDigitization (*dds);

    ////////////////////////////////////////////////////////////////

    TSolEVIOFile *f = new TSolEVIOFile("testgem.ev");
    printf("The filename returned is %s\n", f->GetFileName());

    int res;
    
    res = f->Open();

    if( res != 1 ){
	printf("Opening EVIO returned %s\n", res);
	return;
    }

    // Hypothetical background run
    TSolEVIOFile *fback = new TSolEVIOFile("testback.ev");
    fback->Open();
    
    ////////////////////////////////////////////////////////////////

    int  ndata, i;
    TSolGEMData *gd, *gb;
    gendata *gen;

    ddd->InitTree (*dds, "digdemo2.root");
    
    printf("Digitizing events\n");
    ndata = 0;

    int hadback = 1;

    while( f->ReadNextEvent() && hadback ){
	gd = f->GetGEMData();
	gen = f->GetGenData(0);

	// Access to generated vertex and momentum
	gen->GetV();
	gen->GetP();
	gen->GetWeight();

	int nbacktoadd = 1;
	int backidx = 0;
	// Add some number of background events
	while( hadback = fback->ReadNextEvent() && backidx < nbacktoadd ){
	    gb = fback->GetGEMData();

	    // Randomize times based on gate width
	    for( int bidx = 0; bidx < gb->GetNHit(); bidx++ ){
		double timeshift = gRandom->Uniform(-ddd->GetGateWidth(), 75.0 );//ns
		gb->SetHitTime(bidx, gb->GetHitTime(bidx) + timeshift );
	    }

	    gd->AddGEMData(gb);
	    backidx++;
	}

	if( backidx != nbacktoadd ){
	    printf("Warning:  Not enough background events to be added (%d)\n", backidx);
	}

	ddd->Digitize(*gd, *dds);
	ndata++;
	delete gd;
    }
    printf("Completed %d events\n", ndata);

    ddd->WriteTree();
    ddd->CloseTree();

}
