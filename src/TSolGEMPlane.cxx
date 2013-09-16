#include "TSolGEMPlane.h"

#include "TClonesArray.h"

#include "TSolGEMChamber.h"
#include "TSolGEMCluster.h"
#include "TSolWedge.h"
#include "THaEvData.h"
#include "TMath.h"

#include <iostream>
#include <cassert>

using namespace std;

TSolGEMPlane::TSolGEMPlane()
  : THaSubDetector()
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fWedge = new TSolWedge;
  return;
}

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fWedge = new TSolWedge;
  return;
}

TSolGEMPlane::~TSolGEMPlane()
{
  //  delete fClusters;
  delete fWedge;
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
  // Get x/y position, size, and angles from database if and only
  // if parent is null otherwise copy from parent

  // Note that origin is in lab frame, size is in wedge frame.

  Int_t err;
  Double_t torad = atan(1) / 45.0;

  TSolGEMChamber* parent = (TSolGEMChamber*) GetParent();
  Double_t z0;
  Double_t depth;
  if (parent != NULL)
    {
      fOrigin = parent->GetOrigin();
      fSize[0] = (parent->GetSize())[0];
      fSize[1] = (parent->GetSize())[1];
      fSize[2] = (parent->GetSize())[2];
      fWedge->SetGeometry (parent->GetWedge().GetR0(),
			   parent->GetWedge().GetR1(),
			   parent->GetWedge().GetPhi0(),
			   parent->GetWedge().GetDPhi());

      z0 = fOrigin[2];
      depth = fSize[2];
    }
  else
    {
      Double_t r0 = -999.0;
      Double_t r1 = -999.0;
      Double_t phi0 = -999.0;
      Double_t dphi = -999.0;
      const DBRequest request[] = 
	{
	  {"r0",          &r0,           kDouble, 0, 1},
	  {"r1",          &r1,           kDouble, 0, 1},
	  {"phi0",        &phi0,         kDouble, 0, 1},
	  {"dphi",        &dphi,         kDouble, 0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      
      if (err)
	return err;

      // Database specifies angles in degrees, convert to radians
      phi0 *= torad;
      dphi *= torad;
      
      fWedge->SetGeometry (r0, r1, phi0, dphi);

      fOrigin[0] = (fWedge->GetOrigin())[0];
      fOrigin[1] = (fWedge->GetOrigin())[1];
      fSize[0] = (fWedge->GetSize())[0];
      fSize[1] = (fWedge->GetSize())[1];
      z0 = -999.0;
      depth = -999.0;
    }

  const DBRequest request[] = 
    {
      {"stripangle",  &fSAngle,      kDouble, 0, 1},
      {"pitch",       &fSPitch,      kDouble, 0, 1},
      {"z0",          &z0,           kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB( file, date, request, fPrefix );

  if (err)
    return err;

  fSAngle *= torad;

  SetRotations();
  fOrigin[2] = z0;
  fSize[2] = depth;
  
  // Get numbers of strips

  Double_t xs0 = (GetSize())[0];
  Double_t ys0 = (GetSize())[1];
  Double_t xs[4] = {xs0, xs0, -xs0, -xs0};
  Double_t ys[4] = {ys0, -ys0, ys0, -ys0};

  Int_t smin = 1e9;
  Int_t smax = 1e-9;
  for (UInt_t i = 0; i < 4; ++i)
    {
      PlaneToStrip (xs[i], ys[i]);
      Int_t s = (Int_t) (xs[i] / GetSPitch());
      if (s < smin) smin = s;
      if (s > smax) smax = s;
    }
  fNStrips = smax - smin + 1;
  fSBeg = -fNStrips * fSPitch * 0.5;

  return kOK;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

  //    int i = 0;

    //    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}

Double_t 
TSolGEMPlane::GetSAngle()   const
{
  return fSAngle;
}

void
TSolGEMPlane::PlaneToStrip (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCWS * x - fSWS * y;
  y = fSWS * temp + fCWS * y;
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
TSolGEMPlane::StripToPlane (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCWS * x + fSWS * y;
  y = -fSWS * temp + fCWS * y;
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

Double_t TSolGEMPlane::StripNumtoStrip( Int_t strip )
{
    // Gives x coordinate in strip frame of a wire
    return (strip - GetStrip(0.,0.))*GetSPitch();
}


Double_t TSolGEMPlane::StriptoProj( Double_t s )
{
    // Gives coordinate in projection frame from strip frame x
    Double_t r = (GetWedge().GetR1()-GetWedge().GetR0())/2.0;
    return s + r*fCWS;
}


Double_t TSolGEMPlane::StripNumtoProj( Int_t s ){
    // Gives coordinate in projection frame from strip number
    return StriptoProj( StripNumtoStrip(s) );
}

Double_t 
TSolGEMPlane::GetStripLowerEdge (UInt_t is) const {return (fSBeg + is * GetSPitch());}

Double_t 
TSolGEMPlane::GetStripUpperEdge (UInt_t is) const {return GetStripLowerEdge (is) + GetSPitch();}


Int_t
TSolGEMPlane::GetStripUnchecked( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame,
  // no questions asked. Caller must check return value

  return (Int_t) ((x - fSBeg) / GetSPitch());
}

Int_t
TSolGEMPlane::GetStripInRange( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame
  // and, if out of range, limit it to allowable values.

  Int_t s = GetStripUnchecked(x);
  if( s < 0 )              s = 0;
  if( s >= GetNStrips() )  s = GetNStrips()-1;
  return s;
}
    
Int_t
TSolGEMPlane::GetStrip (Double_t x, Double_t yc) const
{
  // Strip number corresponding to coordinates x, y in 
  // strip frame, or -1 if outside (2-d) bounds

  Double_t xc = x;
  StripToLab (xc, yc);

  if (!fWedge->Contains (xc, yc))
    return -1;

  Int_t s = GetStripUnchecked(x);
  assert( s >= 0 && s < GetNStrips() ); // by construction in ReadGeometry()
  return s;
}

void 
TSolGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2)
       << " (rho,theta,phi)=(" << o.Mag() << "," << o.Theta()*TMath::RadToDeg()
       << "," << o.Phi()*TMath::RadToDeg() << ")"
       << endl;

  const Float_t* s = GetSize();
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;

  cout << "  Wedge geometry: r0: " << fWedge->GetR0()
       << " r1: " << fWedge->GetR1()
       << " phi0: " << fWedge->GetPhi0()*TMath::RadToDeg()
       << " dphi: " << fWedge->GetDPhi()*TMath::RadToDeg()
       << endl;

  cout << "  " << GetNStrips() << " strips"
       << ", angle " << GetSAngle()*TMath::RadToDeg()
       << ", start " << fSBeg << " " << 0.5*(GetStripLowerEdge(0)+GetStripUpperEdge(0))
       << ", pitch " << GetSPitch()
       << endl;
}

void
TSolGEMPlane::SetRotations()
{
  // Set rotation angle trig functions

  fCWS = cos (-GetSAngle());
  fCLS = cos (-GetAngle()-GetSAngle());
  fSWS = sin (-GetSAngle());
  fSLS = sin (-GetAngle()-GetSAngle());
}
