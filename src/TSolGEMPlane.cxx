#include "TSolGEMPlane.h"

#include <iostream>

#include "TClonesArray.h"
#include "TSolGEMCluster.h"
#include "THaEvData.h"

using namespace std;

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc )
  : fDir (kGEMX),
    fNStrips (100), // some dummy values for now, will be filled from database
    fSBeg (-1.0),
    fSPitch (0.0001),
    fPixelFactor (2)
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

Double_t 
TSolGEMPlane::GetSAngle()   const
{
  if (fDir == kGEMX)
    return 3.14159/2;
  else if (fDir == kGEMY)
    return 0.0;
  else
    {
      cerr << __FUNCTION__ << " Strip angle undefined" << endl;
      return 9999.0;
    }
}


