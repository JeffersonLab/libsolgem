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
    
    ////////////////////////////////////////////////////////////////

    int  ndata, i;
    TSolGEMData *gd;
    gendata *gen;

    ddd->InitTree (*dds, "digdemo2.root");
    
    printf("Digitizing events\n");
    ndata = 0;
    while( f->ReadNextEvent() ){
	gd = f->GetGEMData();
	gen = f->GetGenData(0);

	// Access to generated vertex and momentum
	gen->GetV();
	gen->GetP();

	ddd->Digitize(*gd, *dds);
	ndata++;
	delete gd;
    }
    printf("Completed %d events\n", ndata);

    ddd->WriteTree();
    ddd->CloseTree();

}
