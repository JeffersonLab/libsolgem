
#ifndef __TSOLWEDGE_H
#define __TSOLWEDGE_H

#include <cmath>
#include <vector>

#include <Rtypes.h>

class TDatime;

// A "wedge" is a section of an annulus. It is characterized by
// the minimum and maximum radius (r0 and r1), the minimum angle (phi0),
// and the angular width (dphi).

// The inner and outer arcs are always centered on (x, y) = (0, 0) in the
// lab frame. 

// Derived from these quantities is a rectangular prism bounding box,
// one of whose sides is parallel to the symmetry axis of the wedge,
// described by a "size" which is half widths, and an "origin" which 
// is the center of the bounding box.

// The "wedge frame" is the frame whose x/y origin is the center of
// the bounding box and whose x axis lies along the symmetry axis of
// the wedge.

// The origin and phi0 are specified in the lab frame. The size is in the
// wedge frame.

class TSolWedge
{
 public:
  TSolWedge (Double_t r0 = 0, Double_t r1 = 999., Double_t phi0 = 0, Double_t dphi = 6.28);
  virtual ~TSolWedge() {};

  Double_t GetR0() const {return sqrt (fR0R0);};
  Double_t GetR1() const {return sqrt (fR1R1);};
  Double_t GetPhi0() const {return fPhi0;};
  Double_t GetDPhi() const {return fPhi1 - fPhi0;};
  std::vector < Double_t > GetOrigin() const {return fOrigin;};
  std::vector < Double_t > GetSize() const {return fSize;};

  std::vector < Double_t > Bounds() const;
  Bool_t Contains (Double_t x, Double_t y) const;

  void SetGeometry (const Double_t r0,
		    const Double_t r1,
		    const Double_t phi0,
		    const Double_t dphi);
    
  Double_t GetAngle() const {return GetPhi0() + GetDPhi()/2;}; // rotation angle between lab and wedge frame

  // Frame conversions
  void LabToWedge (Double_t& x, Double_t& y) const;  // input and output in meters
  void WedgeToLab (Double_t& x, Double_t& y) const;  // input and output in meters

 private:
  void SetRotations();

  Double_t fR0R0;
  Double_t fR1R1;
  Double_t fPhi0;
  Double_t fTPhi0;
  Double_t fPhi1;
  Double_t fTPhi1;

  std::vector < Double_t > fOrigin;   // x, y
  std::vector < Double_t > fSize;     // x, y

  // Trig functions for rotations
  Double_t fCLW; // cos (lab to wedge angle)
  Double_t fSLW; // sin...
};

namespace NSolWedge
{
  const Double_t gPI = 4 * atan (1.0);;
  Int_t CmpAng (Double_t x, Double_t y, Double_t phi, Double_t tphi);
};

#endif
