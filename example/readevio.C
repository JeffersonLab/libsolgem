
//  Simple instantiation of classes
//  Gonna flush this out later to show
//  how to run over data.  This will be
//  our example "replay" script

void readevio(){
    printf("\n** This gets called with 'analyzer' and not 'root' **\n");
    printf("** If you're getting missing symbol errors, this is likely the cause **\n\n");

    gSystem->Load("../libsolgem.so");


    TSolSpec *s = new TSolSpec("solspec", "SoLID Spectrometer with GEMs");
    TSolGEMPlane *p = new TSolGEMPlane("gemplane", "A SolID GEM plane");

    // These will be called elsewhere, this is just demonstratory
    TSolGEMCluster *c = new TSolGEMCluster();

}
