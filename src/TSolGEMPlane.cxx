#include "TSolGEMPlane.h"

#include <iostream>

#include "TClonesArray.h"

#include "TSolGEMChamber.h"
#include "TSolGEMCluster.h"
#include "THaEvData.h"

using namespace std;

TSolGEMPlane::TSolGEMPlane()
  : THaSubDetector()
{
//   fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
//   fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

Int_t 
TSolGEMPlane::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  return kOK;
}

Int_t 
TSolGEMPlane::ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required)
{
  // Get x/y position, size, and frame angle from database if and only
  // if parent is null otherwise copy from parent

  // Note that origin is in lab frame, size is in chamber frame

  Int_t err = kOK;
  TSolGEMChamber* parent = (TSolGEMChamber*) GetParent();
  if (parent == NULL)
    {
      THaSubDetector::ReadGeometry (file, date, false);
      const DBRequest request[] = 
	{
	  {"angle",       &fAngle,       kDouble,    0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      if (err)
	return err;
    }
  else
    {
      fOrigin = parent->GetOrigin();
      fSize[0] = (parent->GetSize())[0];
      fSize[1] = (parent->GetSize())[1];
      fSize[2] = (parent->GetSize())[2];
      fAngle = (parent->GetAngle());
      Double_t z0 = fOrigin[2];
      const DBRequest request[] = 
	{
	  {"z0",          &z0,           kDouble, 0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      if (err)
	return err;
      fOrigin[2] = z0;
    }

  const DBRequest request[] = 
    {
      {"direction",   &fDir,         kInt,    0, 1},
      {"pitch",       &fSPitch,      kDouble, 0, 1},
      {0}
    };
  err = LoadDB( file, date, request, fPrefix );

  if (err)
    return err;

  SetRotations();

  fNStrips = 2 * (GetSize())[fDir] / fSPitch;
  if (2 * (GetSize())[fDir] - fNStrips * fSPitch > 1E-9) 
    fNStrips++;
  fSBeg = -(GetSize())[fDir];

  return kOK;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

    int i = 0;

//     new ((*fClusters)[i]) TSolGEMCluster();

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

void
TSolGEMPlane::LabToChamber (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLC * x - fSLC * y;
  y = fSLC * temp + fCLC * y;
  return;
}

void
TSolGEMPlane::ChamberToStrip (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCCS * x - fSCS * y;
  y = fSCS * temp + fCCS * y;
  return;
}

void
TSolGEMPlane::LabToStrip (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLS * x - fSLS * y;
  y = fSLS * temp + fCLS * y;
  return;
}

void
TSolGEMPlane::StripToChamber (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCCS * x + fSCS * y;
  y = -fSCS * temp + fCCS * y;
  return;
}

void
TSolGEMPlane::ChamberToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLC * x + fSLC * y;
  y = -fSLC * temp + fCLC * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

void
TSolGEMPlane::StripToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLS * x + fSLS * y;
  y = -fSLS * temp + fCLS * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

Int_t
TSolGEMPlane::GetStrip (Double_t x, Double_t y) const
{
  // Strip number corresponding to coordinates x, y in 
  // strip frame, or -1 if outside (2-d) bounds

  // For now strips are either horizontal or vertical, so this is easy

  Double_t xc, yc; // chamber frame
  if (fDir == kGEMX)
    {
      xc = x;
      yc = y;
    }
  else if (fDir == kGEMY)
    {
      xc = y;
      yc = -x;
    }
  else
    {
      cerr << __FUNCTION__ << " Strip angle undefined" << endl;
      return -1;
    }

  // Check if in bounds (with a grace margin)

  if (xc < GetLowerEdgeX()-1e-3 || xc > GetUpperEdgeX()+1e-3 ||
      yc < GetLowerEdgeY()-1e-3 || yc > GetUpperEdgeY()+1e-3)
    return -1;

  Int_t s = (Int_t) ((x - fSBeg) / GetSPitch());
  return s >= 0 ? s : 0;
}

void 
TSolGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2) << endl;

  const Float_t* s = GetSize();
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;
  cout << "  " << GetNStrips() << " strips, pitch " << GetSPitch() << endl;
}

void
TSolGEMPlane::SetRotations()
{
  // Set rotation angle trig functions

  fCLC = -0.157;
  fCLC = cos (-GetAngle());
  fCCS = cos (-GetSAngle());
  fCLS = cos (-GetAngle()-GetSAngle());
  fSLC = sin (-GetAngle());
  fSCS = sin (-GetSAngle());
  fSLS = sin (-GetAngle()-GetSAngle());
}
