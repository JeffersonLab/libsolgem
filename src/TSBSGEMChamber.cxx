#include <iostream>

#include "TSBSGEMChamber.h"
#include "TSBSGEMPlane.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TMath.h"
#include "ha_compiledata.h"

using namespace std;

//Recommanded constructor
TSBSGEMChamber::TSBSGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  // For now at least we just hard wire two chambers
  fNPlanes = 2;
  fPlanes = new TSBSGEMPlane*[fNPlanes];
  fBox = new TSBSBox;

  return;
}

TSBSGEMChamber::~TSBSGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
  delete fBox;
}


const char* TSBSGEMChamber::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
      return Form ("%s.", app->GetName());
    else
      return fPrefix;
}
  
Int_t
TSBSGEMChamber::ReadDatabase (const TDatime& date)
{
  //Read the geometry for the TSBSBox
  //Calls read geometry which, as it name indicates, actually reads the parameters
 
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  err = InitPlane (0, "x", TString (GetTitle()) +" x");
  if( err != kOK ) return err;
  err = InitPlane (1, "y", TString (GetTitle()) +" y");
  if( err != kOK ) return err;

  return kOK;
}

Int_t
TSBSGEMChamber::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{

  Int_t err = THaDetector::ReadGeometry (file, date, required);
  if (err)
    return err;
  
  Double_t dmag = -999.0;
  Double_t d0 = -999.0;
  Double_t xoffset = -999.0;
  Double_t dx = -999.0;
  Double_t dy = -999.0;
  //Double_t thetaH = -999.0;
  Double_t thetaV = -999.0;
  Double_t depth = -999.0;
  const DBRequest request[] =
    {
      {"dmag",        &dmag,         kDouble, 0, 1},
      {"d0",          &d0,           kDouble, 0, 1},
      {"xoffset",     &xoffset,      kDouble, 0, 1},
      {"dx",          &dx,           kDouble, 0, 1},
      {"dy",          &dy,           kDouble, 0, 1},
      //{"thetaH",      &thetaH,       kDouble, 0, 1},
      {"thetaV",      &thetaV,       kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB (file, date, request, fPrefix);
  if (err)
    return err;
  
  // Database specifies angles in degrees, convert to radians
  Double_t torad = atan(1) / 45.0;
  //thetaH *= torad;
  thetaV *= torad;
  
  // Set the geometry for the dependent TSBSBox
  fBox->SetGeometry (dmag, d0, xoffset, dx, dy, //thetaH, 
		     thetaV);

  fOrigin[0] = (fBox->GetOrigin())[0];
  fOrigin[1] = (fBox->GetOrigin())[1];
  fOrigin[2] = (fBox->GetOrigin())[2];
  fSize[0] = fBox->GetDX();
  fSize[1] = fBox->GetDY();
  fSize[2] = depth;

  return kOK;
}

//
Int_t
TSBSGEMChamber::Decode (const THaEvData& ed )
{
  for (UInt_t i = 0; i < GetNPlanes(); ++i)
    {
      GetPlane (i).Decode (ed);//"Neutralized" function: does nothing and returns 0.
    }
  return 0;
}

Int_t
TSBSGEMChamber::InitPlane (const UInt_t i, const char* name, const char* desc)
{
  //Initialize TSBSGEMPlane number i
  fPlanes[i] = new TSBSGEMPlane (name, desc, this);
  fPlanes[i]->SetName (name);
  return fPlanes[i]->Init();
}

void
TSBSGEMChamber::Print (const Bool_t printplanes)
{
  //Print TSBSGEMChamber (and dependent TSBSGEMPlanes) info
  cout << "I'm a GEM chamber named " << GetName() << endl;
  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2)
       << " (rho,theta,phi)=(" << o.Mag() << "," << o.Theta()*TMath::RadToDeg()
       << "," << o.Phi()*TMath::RadToDeg() << ")"
       << endl;

#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
  const Double_t* s = GetSize();
#else
  const Float_t* s = GetSize();
#endif
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;

  cout << "  Box geometry: D0: " << fBox->GetD0()
       << " XOffset: " << fBox->GetXOffset()
       << " DX: " << fBox->GetDX()
       << " DY: " << fBox->GetDY()
    //<< " thetaH: " << fBox->GetThetaH()*TMath::RadToDeg()
       << " thetaV: " << fBox->GetThetaV()*TMath::RadToDeg()
       << endl;
  
  double x = o(0);
  double y = o(1);
  double z = o(2);
  fBox->HallCenterToSpec(x, y, z);
  cout << "Origin: recomputed: " << x << " " << y << " " << z << endl;
  
  if (printplanes){
    cout << "I contain the " << GetNPlanes() << " following planes: " << endl;
    for (UInt_t i = 0; i < GetNPlanes(); ++i)
      {
	fPlanes[i]->Print(!printplanes);
      }
  }
}

//Transformations: frame conversions
void  
TSBSGEMChamber::HallCenterToPlane (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  HallCenterToPlane (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in mm

void  
TSBSGEMChamber::HallCenterToSpec (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  HallCenterToSpec (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in mm

void 
TSBSGEMChamber::LabToPlane (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  LabToPlane (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void  
TSBSGEMChamber::HallCenterToLab (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  HallCenterToLab (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void 
TSBSGEMChamber::LabToSpec (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  LabToSpec (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void 
TSBSGEMChamber::SpecToPlane (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  SpecToPlane (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void 
TSBSGEMChamber::PlaneToSpec (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  PlaneToSpec (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void 
TSBSGEMChamber::SpecToLab (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  SpecToLab (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void  
TSBSGEMChamber::LabToHallCenter (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  LabToHallCenter (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in mm

void 
TSBSGEMChamber::PlaneToLab (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  PlaneToLab (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in meters

void  
TSBSGEMChamber::SpecToHallCenter (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  SpecToHallCenter (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in mm

void  
TSBSGEMChamber::PlaneToHallCenter (TVector3& X_) const {
  double x = X_[0];
  double y = X_[1];
  double z = X_[2];
  
  PlaneToHallCenter (x, y, z);
  
  X_[0] = x;
  X_[1] = y;
  X_[2] = z;
};  // input and output in mm
