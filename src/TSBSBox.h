
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
// Its main utility is to provide the transformation matrices 
// from the GEM plane frames to the spectrometer frames or to the lab.

// It is basically just a box, characterized by an extension in x (dx) and y (dy)
// and by five parameters for its location: 
// dmag, the distance of the magnet to the target, 
// which is taken as a refence for the spectrometer coordinates 
// d0, the distance of the box center to the spectrometer reference, 
// xOffset, the translation in X (in transport coordinates) of the chamber, and two angles: 
// one of horizontal rotation (thetaH), which translates the SBS angle,
// one of vertical rotation (thetaV), translating the "bending" of the SBS wrt the xOz plane.
// The "Zero" of the spectrometer is defined as: ( -dmag*sin(thetaH), 0, dmag*cos(thetaH) ); 
//
// Warning: the x direction of the box is actually in the -y direction of the lab, 
//          the y direction of the box being in the + x direction of the lab 

// the box is centered on the "central ray" which is then rotated by thetaH and thetaV.
// To illustrate this, let's take a dumb example: 
// if thetaH = thetaV = 0; then the center of the box would be located on the Z axis, 
// and the box coverage on the x (resp y) direction would be from -dx/2 to +dx/2 
// (resp -dy/2 to +dy/2)

// The rotation to place the box in the lab is made the following way.
// first the box is rotated by thetaH wrt y direction (i.e. x, z, modified, y conserved)
// then, the box is rotated by thetaV wrt x' direction obtained at the previous step 
// (i.e. y, z' modified, x' conserved)
// lastly, the frame is rotated by 90 degrees wrt z" direction (i.e. x', y', modified, z conserved)

// The box z location is to be understood as the box location on the z" axis 
// obtained with the rotation defined before, and is made wrt x, y = 0, 0.


class TSBSBox
{
 public:
  // Default constructors and destructors.
  TSBSBox (Double_t dmag = 1.0, Double_t d0 = 1.0, Double_t xoffset = 0.0, 
	   Double_t dx = 1.0, Double_t dy = 1.0, 
	   Double_t thetaH = 0.0, Double_t thetaV = 0.0);
  virtual ~TSBSBox() {};
  
  // Members getters
  Double_t GetDMag() const {return fDMag;};
  Double_t GetD0() const {return fD0;};
  Double_t GetXOffset() const {return fXOffset;};
  Double_t GetDX() const {return fDX;};
  Double_t GetDY() const {return fDY;};
  Double_t GetThetaH() const {return fThetaH;};
  Double_t GetThetaV() const {return fThetaV;};
  
  TVector3 GetOrigin() const {return fOrigin;};
  //TVector3 GetSize() const {return fSize;};
  
  // One global member setter. 
  // The box geometry is not supposed to be modified along the way
  void SetGeometry (const Double_t dmag,
		    const Double_t d0,
		    const Double_t xoffset,
		    const Double_t dx,
		    const Double_t dy,
		    const Double_t thetaH,
		    const Double_t thetaV);
    
  // Frame conversions methods:  Lab -> Box, Lab -> Spec -> Box, and reverse transformations
  // NB: "Spec" coordinates are equivalent to "transport" coordinates.
  void HallCenterToBox (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void HallCenterToSpec (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void LabToBox (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void HallCenterToLab (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void LabToSpec (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void SpecToBox (Double_t& x, Double_t& y) const; // const {return;};  // input and output in mm
  void BoxToSpec (Double_t& x, Double_t& y) const; // const {return;};  // input and output in mm
  void SpecToLab (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void LabToHallCenter (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void BoxToLab (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void SpecToHallCenter (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  void BoxToHallCenter (Double_t& x, Double_t& y, Double_t& z) const;  // input and output in mm
  
  // Checker: does the plane contain the hit position ?
  Bool_t Contains (Double_t x, Double_t y) const;//hit position in the box
  Bool_t Contains (Double_t x, Double_t y, Double_t z) const;//hit position in the lab.
  
 private:
  // evaluate the rotation matrices
  void SetRotations();
  
  // Members
  Double_t fDMag;
  Double_t fD0;
  Double_t fXOffset;
  Double_t fDX;
  Double_t fDY;
  Double_t fThetaH;
  Double_t fThetaV;

  TVector3 fOrigin;
  //TVector3 fSize;
  
  // Matrices for rotations
  TMatrixD* fRotMat_BL; 
  TMatrixD* fRotMat_LB; 
};

#endif
