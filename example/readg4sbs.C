#define DEBUG true
//  Simple instantiation of classes
//  Gonna flesh this out later to show
//  how to run over data.  This will be
//  our example "replay" script

void readg4sbs(){
    printf("\n** This gets called with 'analyzer' and not 'root' **\n");
    printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

    gSystem->Load("../libsolgem.so");
    gSystem->Load("../../libg4sbsroot.so");

    TSolSimG4SBSFile *f = new TSolSimG4SBSFile("/group/exjpsi/eric/14001/elastic_0.root");
    printf("The filename returned is %s\n", f->GetFileName());
    
    int res;
    
    //cout << DEBUG << endl;
    
    res = f->Open();

    if( res != 1 ){
	printf("Opening G4SBS file returned %s\n", res);
	return;
    }

    int  ndata, i;

    TSolGEMData *gd;
    
    while( f->ReadNextEvent() && f->GetEvNum()<1){
      //if(f->GetEvNum()%1000==0)
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
