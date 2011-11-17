#include "TSolWedge.h"

#include <iostream>
using namespace std;

#include <TDatime.h>

using namespace NSolWedge;

TSolWedge::TSolWedge (Double_t r0, Double_t r1, Double_t phi0, Double_t dphi)
  : fR0R0(0),
    fR1R1(1),
    fPhi0(0),
    fTPhi0(0),
    fPhi1(gPI),
    fTPhi1(0)
{
  SetGeometry (r0, r1, phi0, dphi);
}

std::vector < Double_t > 
TSolWedge::Bounds() const
{
  // Return bounds: xmin, ymin, xmax, ymax
  // Each is either +-r1 or a coordinate of a corner.

  std::vector < Double_t > v(4); 

  Double_t r0 = sqrt (fR0R0);
  Double_t r1 = sqrt (fR1R1);

  Double_t x[4] = {r0 * cos (fPhi0),
		   r0 * cos (fPhi1),
		   r1 * cos (fPhi0),
		   r1 * cos (fPhi1)};
  Double_t y[4] = {r0 * sin (fPhi0),
		   r0 * sin (fPhi1),
		   r1 * sin (fPhi0),
		   r1 * sin (fPhi1)};
  v[0] = x[0];
  v[1] = y[0];
  v[2] = x[0];
  v[3] = y[0];
  for (UInt_t i = 1; i < 4; ++i)
    {
      if (x[i] < v[0]) v[0] = x[i];
      if (y[i] < v[1]) v[1] = y[i];
      if (x[i] > v[2]) v[2] = x[i];
      if (y[i] > v[3]) v[3] = y[i];
    }

  if (fPhi0 < gPI && fPhi1 > gPI) v[0] = -r1;
  if (fPhi0 < 3*gPI/2 && fPhi1 > 3*gPI/2) v[1] = -r1;
  if (fPhi0 < 0 && fPhi1 > 0) v[2] = r1;
  if (fPhi0 < gPI/2 && fPhi1 > gPI/2) v[3] = r1;

  return v;
}

Bool_t
TSolWedge::Contains (Double_t x, Double_t y) const
{
  // Is (x, y) within the wedge?

  Double_t r2 = x * x + y * y;
  if (r2 < fR0R0 || r2 > fR1R1)
    return false;

  if (fPhi0 < 0)
    // Wedge straddles +x axis, need OR
    return (CmpAng (x, y, fPhi0 + 2 * gPI, fTPhi0) >= 0 ||
	    CmpAng (x, y, fPhi1, fTPhi1) <= 0);
  else
    // Doesn't, need AND
    return (CmpAng (x, y, fPhi0, fTPhi0) >= 0 &&
	    CmpAng (x, y, fPhi1, fTPhi1) <= 0);
}


Int_t
NSolWedge::CmpAng (Double_t x, Double_t y, Double_t phi, Double_t tphi)
{
  // x and y are coordinates; phi is an angle in 0 <= phi < 2pi; tphi
  // is its tangent.  Return -1, 0, 1 if angle between (x, y) and +x
  // axis is less than, equal to, greater than phi.

  // First check quadrants and return if not same -- or if point is
  // on y axis.

  if (x > 0 && y >= 0)
    {
      if (phi >= gPI/2) return -1;
    }
  else if (x <= 0 && y > 0)
    {
      if (phi < gPI/2) return 1;
      if (phi >= gPI) return -1;
      if (x == 0) return phi == gPI/2 ? 0 : -1;
    }
  else if (x < 0 && y <= 0)
    {
      if (phi < gPI) return 1;
      if (phi >= 3*gPI/2) return -1;
    }
  else if (x >= 0 && y < 0)
    {
      if (phi < 3*gPI/2) return 1;
      if (x == 0) return phi == 3*gPI/2 ? 0 : -1;
    }
  // In the same quadrant, compare tangent

  Double_t r = y / x;
  return r < tphi ? -1 : (r == tphi ? 0 : 1);
}

void
TSolWedge::SetGeometry (const Double_t r0,
			const Double_t r1,
			const Double_t phi0,
			const Double_t dphi)
{
  // Note that origin is in lab frame, size is in wedge frame.

  fR0R0 = r0 * r0;
  fR1R1 = r1 * r1;

  fPhi0 = -dphi/2; // for now, for bounds calculation in wedge frame
  fPhi1 = fPhi0 + dphi;
  vector <Double_t> bounds = Bounds();

  fPhi0 = phi0;  // now set the actual phi0 in lab frame
  fTPhi0 = cos (fPhi0) == 0 ? 9999. : sin (fPhi0) / cos (fPhi0);
  fPhi1 = phi0 + dphi;
  fTPhi1 = cos (fPhi1) == 0 ? 9999. : sin (fPhi1) / cos (fPhi1);
  SetRotations();

  // Origin is center of bounding box, converted to lab frame
  fOrigin.resize (2);
  fOrigin[0] = 0; // ... but without origin shift of course!
  fOrigin[1] = 0;
  Double_t xc = (bounds[2]+bounds[0])/2;
  Double_t yc = (bounds[3]+bounds[1])/2;
  WedgeToLab (xc, yc);
  fOrigin[0] = xc;
  fOrigin[1] = yc;
  
  // Size is half size of bounding box, in wedge frame
  fSize.resize (2);
  fSize[0] = (bounds[2] - bounds[0]) / 2;
  fSize[1] = (bounds[3] - bounds[1]) / 2;
}

void
TSolWedge::LabToWedge (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLW * x - fSLW * y;
  y = fSLW * temp + fCLW * y;
  return;
}

void
TSolWedge::WedgeToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLW * x + fSLW * y;
  y = -fSLW * temp + fCLW * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

void
TSolWedge::SetRotations()
{
  // Set rotation angle trig functions

  fCLW = cos (-GetAngle());
  fSLW = sin (-GetAngle());
}
