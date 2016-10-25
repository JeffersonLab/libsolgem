#ifndef __TSBSGEMCHAMBER_H
#define __TSBSGEMCHAMBER_H

#include "THaDetector.h"
#include "TSBSBox.h"

class THaEvData;
class TSBSGEMPlane;

// A chamber is a "wedge" (section of an annulus). It is characterized by
// the minimum and maximum radius (r0 and r1), the minimum angle (phi0),
// the angular width (dphi), the z coordinate of the front face (z0),
// and the thickness in z (depth)

// The inner and outer arcs are always centered on (x, y) = (0, 0) in the
// lab frame. 

// Derived from these quantities is a rectangular prism bounding box,
// one of whose sides is parallel to the symmetry axis of the wedge,
// described by a "size" which is half the transverse sizes and the
// full z size, and an "origin" which is the center of the front face of the 
// bounding box.

// The "wedge frame" or "plane frame" is the frame whose x/y origin is
// the center of the bounding box and whose x axis lies along the
// symmetry axis of the wedge.

// The origin and phi0 are specified in the lab frame. The size is in the
// wedge frame.

class TSBSGEMChamber : public THaDetector {
 public:
  TSBSGEMChamber(const char *name, const char *desc);
  TSBSGEMChamber() : fPlanes(0), fNPlanes(0), fBox(0) {} // for ROOT RTTI

  virtual ~TSBSGEMChamber();

  Int_t ReadDatabase (const TDatime& date);
  Int_t ReadGeometry (FILE* file, const TDatime& date, Bool_t required);
  const char* GetDBFileName() const;

  Int_t Decode( const THaEvData & );

  // Return positions of chamber edges, in its own reference frame, in meters
  Double_t GetLowerEdgeX() const {return -(GetSize())[0]/2;}
  Double_t GetLowerEdgeY() const {return -(GetSize())[1]/2;}
  Double_t GetUpperEdgeX() const {return +(GetSize())[0]/2;}
  Double_t GetUpperEdgeY() const {return +(GetSize())[1]/2;}

  TSBSBox& GetBox() const {return *fBox;};
  
  // Frame conversions
  void LabToPlane (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->LabToBox (x, y, z);
  };  // input and output in meters
  void PlaneToLab (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->BoxToLab (x, y, z);
  };  // input and output in meters

  void LabToSpec (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->LabToSpec (x, y, z);
  };  // input and output in meters
  void SpecToPlane (Double_t& x, Double_t& y, Double_t& z) const {
    z = z-fBox->GetD0();
    fBox->SpecToBox (x, y);
  };  // input and output in meters
  void PlaneToSpec (Double_t& x, Double_t& y, Double_t& z) const {
    z = z+fBox->GetD0();
    fBox->BoxToSpec (x, y);
  };  // input and output in meters
  void SpecToLab (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->SpecToLab (x, y, z);
  };  // input and output in meters

  void LabToPlane (TVector3& X_) const;
  void PlaneToLab (TVector3& X_) const;
  
  void LabToSpec (TVector3& X_) const;
  void SpecToPlane (TVector3& X_) const;
  void PlaneToSpec (TVector3& X_) const;
  void SpecToLab (TVector3& X_) const;
  
  UInt_t GetNPlanes() const {return fNPlanes;};

  TSBSGEMPlane& GetPlane (UInt_t i) const {return *(fPlanes[i]);};
  Int_t InitPlane (const UInt_t i, const char* name, const char* desc);
  void Print (const Bool_t printplanes = kTRUE);

 private:
  TSBSGEMPlane** fPlanes; // List of chambers
  UInt_t fNPlanes;
  TSBSBox* fBox;  // Box geometry

  ClassDef(TSBSGEMChamber,0)

};

#endif//__TSBSGEMCHAMBER_H
