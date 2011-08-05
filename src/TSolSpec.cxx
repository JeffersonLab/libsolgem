#include "TSolSpec.h"
#include "TSolGEMCluster.h"
#include "TSolGEMPlane.h"
#include <stdio.h>

TSolSpec::TSolSpec(const char* name, const char* desc )
    :THaSpectrometer(name,desc) {

	printf("Hello, I'm a spectrometer named %s\n", name);
	



	return;
}

Int_t TSolSpec::ReadDatabase( const TDatime& date ){
	// Make a bunch of planes based on a database specification
	// (simple interface is written in the analyzer, just need
	// to implement it)

    fNplanes = 0;
    fPlanes = new  TSolGEMPlane *[fNplanes];

    // ...

    return 0;
}

Int_t TSolSpec::CoarseTrack(){
    int i,j;

    // Assume decoding is done.  You can get the clustered hits with
    // calles like

    TSolGEMCluster *c;
    Double_t z;

    // loop over X planes
    if( GetPlane(i)->GetDirection() == kGEMX ){
	c = (TSolGEMCluster *) GetPlane(i)->GetClusters()->AddrAt(j);
	z = GetPlane(i)->GetOrigin().Z();
	
	c->GetPos();
	c->GetE();
    }

    return 0;
}

Int_t TSolSpec::CoarseReconstruct(){
    return 0;
}

Int_t TSolSpec::Track(){
    return 0;
}


Int_t TSolSpec::Reconstruct(){
    return 0;
}

Int_t TSolSpec::FindVertices(TClonesArray &a){
    return 0;
}
