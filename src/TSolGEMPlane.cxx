#include "TSolGEMPlane.h"
#include "THaEvData.h"

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc ){
    printf("I'm a GEM plane named %s\n", name );

    return;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    return 0;
}
