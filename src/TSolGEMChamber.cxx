#include <iostream>

#include "TSolGEMChamber.h"
#include "TSolWedge.h"

using namespace std;

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  // For now at least we just hard wire two chambers
  fNPlanes = 2;
  fPlanes = new TSolGEMPlane*[fNPlanes];
  fWedge = new TSolWedge;
  
  return;
}

TSolGEMChamber::~TSolGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
  delete fWedge;
}


const char* TSolGEMChamber::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
	return app->GetName();
    else
	return fPrefix;
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

  InitPlane (0, TString (GetName()) + "x", TString (GetTitle()) +" x");
  InitPlane (1, TString (GetName()) + "y", TString (GetTitle()) +" y");

  return kOK;
}

Int_t 
TSolGEMChamber::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{

  Int_t err = THaDetector::ReadGeometry (file, date, required);
  if (err)
    return err;

  Double_t r0 = -999.0;
  Double_t r1 = -999.0;
  Double_t phi0 = -999.0;
  Double_t dphi = -999.0;
  Double_t z0 = -999.0;
  Double_t depth = -999.0;
  const DBRequest request[] = 
    {
      {"r0",          &r0,           kDouble, 0, 1},
      {"r1",          &r1,           kDouble, 0, 1},
      {"phi0",        &phi0,         kDouble, 0, 1},
      {"dphi",        &dphi,         kDouble, 0, 1},
      {"z0",          &z0,           kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB (file, date, request, fPrefix);
  
  if (err)
    return err;
  
  // Database specifies angles in degrees, convert to radians
  Double_t torad = atan(1) / 45.0;
  phi0 *= torad;
  dphi *= torad;

  fWedge->SetGeometry (r0, r1, phi0, dphi);

  fOrigin[0] = (fWedge->GetOrigin())[0];
  fOrigin[1] = (fWedge->GetOrigin())[1];
  fOrigin[2] = z0;
  fSize[0] = (fWedge->GetSize())[0];
  fSize[1] = (fWedge->GetSize())[1];
  fSize[2] = depth;

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

  cout << "  Wedge geometry: r0: " << fWedge->GetR0()
       << " r1: " << fWedge->GetR1()
       << " phi0: " << fWedge->GetPhi0()
       << " dphi: " << fWedge->GetDPhi() << endl;

  if (printplanes)
    for (UInt_t i = 0; i < GetNPlanes(); ++i)
      {
	fPlanes[i]->Print();
      }
}
