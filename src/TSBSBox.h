
#ifndef __TSBSBOX_H
#define __TSBSBOX_H

#include <cmath>
//#include <vector>

#include <Rtypes.h>
#include <TVector3.h>
#include <TMatrixD.h>

class TDatime;

// This class is derivated from the "TSolWedge" class, 
// to implement the SBS geometry, instead of the SoLID one. 
// Its main utility is to provide the transformation matrices from the GEM frames to the lab.

// It is basically just a box, characterized by an extension in x (dx) and y (dy)
// and by two angles: 
// one of horizontal rotation (thetaH), which translates the SBS angle,
// one of vertical rotation (thetaV), translating the "bending" of the SBS wrt the xOz plane.

// the box is centered on the "central ray" which is then rotated by thetaH and thetaV.
// To illustrate this, let's take a dumb example: 
// if thetaH = thetaV = 0; then the center of the box would be located on the Z axis, 
// and the box coverage on the x (resp y) direction would be from -dx/2 to +dx/2 
// (resp -dy/2 to +dy/2)

// The rotation is made the following way: 
// first the box is rotated by thetaH wrt y direction (i.e. x, z, modified, y conserved)
// then, the box is rotated by thetaV wrt x' direction obtained at the previous step 
// (i.e. y, z' modified, x' conserved)

// The box z location is to be understood as the box location on the z" axis 
// obtained with the rotation defined before, and is made wrt x, y = 0, 0.


class TSBSBox
{
 public:
  TSBSBox (Double_t d0 = 1.0, Double_t dx = 1.0, Double_t dy = 1.0, 
	   Double_t thetaH = 0.0, Double_t thetaV = 0.0);
  virtual ~TSBSBox() {};
  
  Double_t GetD0() const {return fD0;};
  Double_t GetDX() const {return fDX;};
  Double_t GetDY() const {return fDY;};
  Double_t GetThetaH() const {return fThetaH;};
  Double_t GetThetaV() const {return fThetaV;};
  
  TVector3 GetOrigin() const {return fOrigin;};
  TVector3 GetSize() const {return fSize;};
  
  void SetGeometry (const Double_t d0,
		    const Double_t dx,
		    const Double_t dy,
		    const Double_t thetaH,
		    const Double_t thetaV);
    
  // Frame conversions
  void LabToBox (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in meters
  void BoxToLab (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in meters

  Bool_t Contains (Double_t x, Double_t y, Double_t z) const;
  
 private:
  void SetRotations();

  Double_t fD0;
  Double_t fDX;
  Double_t fDY;
  Double_t fThetaH;
  Double_t fThetaV;

  TVector3 fOrigin;
  TVector3 fSize;
  
  // matrices for rotations
  TMatrixD* fRotMat_BL; 
  TMatrixD* fRotMat_LB; 
};

#endif
