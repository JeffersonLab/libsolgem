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
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  Print(kFALSE);
  cout << "^^^^^" << endl;

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
TSolGEMChamber::Print (const Bool_t printplanes)
{
  cout << "I'm a GEM chamber named " << GetName() << endl;
  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2) << endl;

  const Float_t* s = GetSize();
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;

  if (printplanes)
    for (UInt_t i = 0; i < GetNPlanes(); ++i)
      {
	fPlanes[i]->Print();
      }
}
