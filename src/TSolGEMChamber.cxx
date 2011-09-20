#include "TSolGEMChamber.h"

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
    printf("I'm a GEM chamber named %s\n", name );

    fNPlanes = 2;
    fPlanes = new TSolGEMPlane*[fNPlanes];
    for (UInt_t i = 0; i < fNPlanes; ++i)
      fPlanes[i] = NULL;

    return;
}

TSolGEMChamber::~TSolGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
}
  

Int_t 
TSolGEMChamber::Decode (const THaEvData& ed)
{
  for (UInt_t i = 0; i < GetNPlanes(); ++i)
    {
      GetPlane (i).Decode (ed);
    }
  return 0;
}

void 
TSolGEMChamber::InitPlane (UInt_t i, char* name, char* desc)
{
  fPlanes[i] = new TSolGEMPlane (name, desc, this);
}

