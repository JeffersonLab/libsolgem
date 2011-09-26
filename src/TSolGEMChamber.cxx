#include <iostream>

#include "TSolGEMChamber.h"

using namespace std;

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  // For now at least we just hard wire two chambers
    fNPlanes = 2;
    fPlanes = new TSolGEMPlane*[fNPlanes];

    return;
}

TSolGEMChamber::~TSolGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
}
  
Int_t 
TSolGEMChamber::ReadDatabase (const TDatime& date)
{
  // This should be a database read but actually we're just hard wiring
  // this for now

  InitPlane (0, TString (GetName()) + "x", TString (GetTitle()) +" x");
  InitPlane (1, TString (GetName()) + "y", TString (GetTitle()) +" y");

  return kOK;
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
TSolGEMChamber::InitPlane (const UInt_t i, const char* name, const char* desc)
{
  fPlanes[i] = new TSolGEMPlane (name, desc, this);
  fPlanes[i]->SetName (name);
  fPlanes[i]->Init();
}

void
TSolGEMChamber::Print()
{
  cout << "I'm a GEM chamber named " << GetName() << endl;
  fPlanes[0]->Print();
  fPlanes[1]->Print();
}

  
