#include "TClonesArray.h"
#include "TSolGEMPlane.h"
#include "TSolGEMCluster.h"
#include "THaEvData.h"

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc ){
    printf("I'm a GEM plane named %s\n", name );

    fClusters = new TClonesArray("TSolGEMCluster", 100);

    return;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

    int i;

    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}
