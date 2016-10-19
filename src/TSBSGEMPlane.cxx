#include "TSBSGEMPlane.h"

#include "TClonesArray.h"

#include "TSBSGEMChamber.h"
#include "TSolGEMCluster.h"
#include "TSBSBox.h"
#include "THaEvData.h"
#include "TMath.h"
#include "ha_compiledata.h"

#include <iostream>
#include <cassert>

using namespace std;

TSBSGEMPlane::TSBSGEMPlane()
  : THaSubDetector()
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fBox = new TSBSBox;
  return;
}

TSBSGEMPlane::TSBSGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fBox = new TSBSBox;
  return;
}

TSBSGEMPlane::~TSBSGEMPlane()
{
  //  delete fClusters;
  delete fBox;
}

Int_t 
TSBSGEMPlane::ReadDatabase (const TDatime& date)
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
TSBSGEMPlane::ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required)
{
  // Get x/y position, size, and angles from database if and only
  // if parent is null otherwise copy from parent

  // Note that origin is in lab frame, size is in wedge frame.

  Int_t err;
  Double_t torad = atan(1) / 45.0;

  TSBSGEMChamber* parent = (TSBSGEMChamber*) GetParent();
  // Double_t z0;
  Double_t depth;
  if (parent != NULL)
    {
      // fOrigin = parent->GetOrigin();
      // fSize[0] = (parent->GetSize())[0];
      // fSize[1] = (parent->GetSize())[1];
      // fSize[2] = (parent->GetSize())[2];
      fBox->SetGeometry (parent->GetBox().GetD0(),
      			 parent->GetBox().GetDX(),
      			 parent->GetBox().GetDY(),
      			 parent->GetBox().GetThetaH(),
      			 parent->GetBox().GetThetaV());
      
      // z0 = fOrigin[2];
      depth = fSize[2];
    }
  else
    {
      Double_t d0 = -999.0;
      Double_t dx = -999.0;
      Double_t dy = -999.0;
      Double_t thetaH = -999.0;
      Double_t thetaV = -999.0;
      const DBRequest request[] = 
	{
	  {"d0",          &d0,           kDouble, 0, 1},
	  {"dx",          &dx,           kDouble, 0, 1},
	  {"dy",          &dy,           kDouble, 0, 1},
	  {"thetaH",      &thetaH,       kDouble, 0, 1},
	  {"thetaV",      &thetaV,       kDouble, 0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      
      if (err)
	return err;

      // Database specifies angles in degrees, convert to radians
      thetaH *= torad;
      thetaV *= torad;
      
      fBox->SetGeometry (d0, dx, dy, thetaH, thetaV);

      fOrigin[0] = (fBox->GetOrigin())[0];
      fOrigin[1] = (fBox->GetOrigin())[1];
      fOrigin[2] = (fBox->GetOrigin())[2];
      fSize[0] = (fBox->GetSize())[0];
      fSize[1] = (fBox->GetSize())[1];
      depth = -999.0;
    }

  const DBRequest request[] = 
    {
      {"stripangle",  &fSAngle,      kDouble, 0, 1},
      {"pitch",       &fSPitch,      kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB( file, date, request, fPrefix );

  if (err)
    return err;

  fSAngle *= torad;

  SetRotations();
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

Int_t TSBSGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

  //    int i = 0;

    //    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}

Double_t 
TSBSGEMPlane::GetSAngle()   const
{
  return fSAngle;
}

void
TSBSGEMPlane::LabToStrip (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x, y, z};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_LS), m_temp);
  
  x = m_res(0, 0) - fOrigin.X();
  y = m_res(1, 0) - fOrigin.Y();
  z = m_res(2, 0) - fOrigin.Z();
  
  return;
}

void
TSBSGEMPlane::StripToPlane (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCBS * x + fSBS * y;
  y = -fSBS * temp + fCBS * y;
  return;
}

void
TSBSGEMPlane::StripToLab (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x+fOrigin.X(), y+fOrigin.Y(), z+fOrigin.Z()};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_LS), m_temp);
  
  x = m_res(0, 0);
  y = m_res(1, 0);
  z = m_res(2, 0);
  
  return;
}

Double_t TSBSGEMPlane::StripNumtoStrip( Int_t strip )
{
    // Gives x coordinate in strip frame of a wire
  return (strip - GetStrip(0.,0., 0.))*GetSPitch();// attention, c'est peut-etre debile...
  //ou peut-etre pas...
}


Double_t TSBSGEMPlane::StriptoProj( Double_t s )
{
    // Gives coordinate in projection frame from strip frame x
  Double_t x = (GetBox().GetDX())/2.0;
    return s + x*fCBS;
}


Double_t TSBSGEMPlane::StripNumtoProj( Int_t s ){
    // Gives coordinate in projection frame from strip number
    return StriptoProj( StripNumtoStrip(s) );
}

Double_t 
TSBSGEMPlane::GetStripLowerEdge (UInt_t is) const {return (fSBeg + is * GetSPitch());}

Double_t 
TSBSGEMPlane::GetStripUpperEdge (UInt_t is) const {return GetStripLowerEdge (is) + GetSPitch();}


Int_t
TSBSGEMPlane::GetStripUnchecked( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame,
  // no questions asked. Caller must check return value

  return (Int_t) ((x - fSBeg) / GetSPitch());
}

Int_t
TSBSGEMPlane::GetStripInRange( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame
  // and, if out of range, limit it to allowable values.

  Int_t s = GetStripUnchecked(x);
  if( s < 0 )              s = 0;
  if( s >= GetNStrips() )  s = GetNStrips()-1;
  return s;
}
    
Int_t
TSBSGEMPlane::GetStrip (Double_t x, Double_t yc, Double_t zc) const
{
  // Strip number corresponding to coordinates x, y in 
  // strip frame, or -1 if outside (2-d) bounds

  Double_t xc = x;
  StripToLab (xc, yc, zc);

  if (!fBox->Contains (xc, yc, zc))
    return -1;

  Int_t s = GetStripUnchecked(x);
  assert( s >= 0 && s < GetNStrips() ); // by construction in ReadGeometry()
  return s;
}

void 
TSBSGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

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

  cout << "  Box geometry: D0 " << fBox->GetD0()
       << " theta_H: " << fBox->GetThetaH()*TMath::RadToDeg()
       << " theta_V: " << fBox->GetThetaV()*TMath::RadToDeg() 
       << " DX: " << fBox->GetDX()
       << " DY: " << fBox->GetDY()
       << endl;

  cout << "  " << GetNStrips() << " strips"
       << ", angle " << GetSAngle()*TMath::RadToDeg()
       << ", start " << fSBeg << " " << 0.5*(GetStripLowerEdge(0)+GetStripUpperEdge(0))
       << ", pitch " << GetSPitch()
       << endl;
}

void
TSBSGEMPlane::SetRotations()
{
  // Set rotation angle trig functions
  fCBS = cos (-GetSAngle());
  fSBS = sin (-GetSAngle());

  Double_t thetaH = fBox->GetThetaH();
  Double_t thetaV = fBox->GetThetaV();

  Double_t arr_roty[9] = {cos(thetaH),  0, sin(thetaH),
			  0,             1,            0,
			  -sin(thetaH), 0, cos(thetaH)};//rotation around hall pivot
  Double_t arr_rotxp[9] = {1,            0,            0,
			   0, cos(thetaV), -sin(thetaV),
			   0, sin(thetaV),  cos(thetaV)};//spectrometer "bending"
  Double_t arr_rotsp[9] = {fCBS, fSBS, 0,
			   fSBS, -fCBS, 0,
			   0,        0, 1};
  //rotation to strip; assuming that by convention, x is orthogonal to the strips' extension 
    
  TMatrixD Roty(3,3,arr_roty);// rotation along y: spectrometer theta
  TMatrixD Rotxp(3,3,arr_rotxp);// ratotion along x': spectrometer bending
  TMatrixD Rotsp(3,3,arr_rotsp);// ratotion along x': spectrometer bending
  TMatrixD Rotlp(3,3, arr_roty);
  fRotMat_SL = new TMatrixD(3,3, arr_roty);
  
  // Set rotation angle trig functions
  Rotlp.Mult(Roty, Rotxp);
  fRotMat_SL->Mult(Rotlp, Rotsp);
  
  fRotMat_LS = fRotMat_SL;
  fRotMat_LS->Invert();
}
