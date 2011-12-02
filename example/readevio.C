
//  Simple instantiation of classes
//  Gonna flesh this out later to show
//  how to run over data.  This will be
//  our example "replay" script

void readevio(){
    printf("\n** This gets called with 'analyzer' and not 'root' **\n");
    printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

    gSystem->Load("../libsolgem.so");

    TSolEVIOFile *f = new TSolEVIOFile("testgem.ev");
    printf("The filename returned is %s\n", f->GetFileName());

    int res;
    
    res = f->Open();

    if( res != 1 ){
	printf("Opening EVIO returned %s\n", res);
	return;
    }

    int  ndata, i;

    TSolGEMData *gd;
    
    while( f->ReadNextEvent() ){
	printf("Event %d\n", f->GetEvNum());

	gd = f->GetGEMData();

	printf("ndata = %d\n", (int) gd->GetNHit() );
	ndata = gd->GetNHit();
	for( i = 0; i < ndata; i++ ){
	    gd->PrintHit(i);
	}
	printf("\n");
    }
}
