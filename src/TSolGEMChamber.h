#ifndef __TSOLGEMCHAMBER_H
#define __TSOLGEMCHAMBER_H

#include "THaDetector.h"
#include "TSolWedge.h"

#include "THaEvData.h"
#include "TSolGEMPlane.h"

class TSolWedge;

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

class TSolGEMChamber : public THaDetector {
 public:
  TSolGEMChamber(const char *name, const char *desc);
  virtual ~TSolGEMChamber();

  Int_t ReadDatabase (const TDatime& date);
  Int_t ReadGeometry (FILE* file, const TDatime& date, Bool_t required);
  const char* GetDBFileName() const;

  Int_t Decode( const THaEvData & );

  // Return positions of chamber edges, in its own reference frame, in meters
  Double_t GetLowerEdgeX() const {return -(GetSize())[0];}
  Double_t GetLowerEdgeY() const {return -(GetSize())[1];}
  Double_t GetUpperEdgeX() const {return +(GetSize())[0];}
  Double_t GetUpperEdgeY() const {return +(GetSize())[1];}

  TSolWedge& GetWedge() const {return *fWedge;};
  Double_t GetAngle() const {return fWedge->GetAngle();}; // rotation angle between lab and wedge frame

  // Frame conversions
  void LabToPlane (Double_t& x, Double_t& y) const {fWedge->LabToWedge (x, y);};  // input and output in meters
  void PlaneToLab (Double_t& x, Double_t& y) const {fWedge->WedgeToLab (x, y);};  // input and output in meters

  UInt_t GetNPlanes() const {return fNPlanes;};

  TSolGEMPlane& GetPlane (UInt_t i) const {return *(fPlanes[i]);};
  void InitPlane (const UInt_t i, const char* name, const char* desc);
  void Print (const Bool_t printplanes = kTRUE);

 private:
  TSolGEMPlane** fPlanes; // List of chambers
  UInt_t fNPlanes;
  TSolWedge* fWedge;  // Wedge geometry

 public:
  ClassDef(TSolGEMChamber,0)

};

#endif//__TSOLGEMCHAMBER_H
