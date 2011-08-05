#include "TSolSpec.h"
#include <stdio.h>

TSolSpec::TSolSpec(const char* name, const char* desc )
    :THaSpectrometer(name,desc) {

	printf("Hello, I'm a spectrometer named %s\n", name);

	return;
}

Int_t TSolSpec::CoarseTrack(){
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
