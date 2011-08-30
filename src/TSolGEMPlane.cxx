#include "TClonesArray.h"
#include "TSolGEMPlane.h"
#include "TSolGEMCluster.h"
#include "THaEvData.h"

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc )
  : fNStrips (100), // some dummy values for now, will be filled from database
    fSBeg (-1.0),
    fSPitch (0.0001),
    fSAngle (0.0)
{
    printf("I'm a GEM plane named %s\n", name );

    fClusters = new TClonesArray("TSolGEMCluster", 100);

    return;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

    int i = 0;

    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}
