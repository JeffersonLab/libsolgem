#ifndef __TSOLGEMCHAMBER_H
#define __TSOLGEMCHAMBER_H

#include "THaDetector.h"
#include "THaEvData.h"

#include "TSolGEMPlane.h"

// A chamber contains one or more (currently hardwired at two) GEM planes.
// It is rectangular, and oriented within a geometric plane of constant Z
// at an arbitrary angle from the x axis. (This is a separate rotation
// angle from the angle of the strips in each plane, which is with respect
// to the chamber's x axis.)

// The geometry is characterized by: Origin (in lab frame), rotation angle
// (with respect to lab x axis), and size (in rotated frame).

class TSolGEMChamber : public THaDetector {
 public:
  TSolGEMChamber(const char *name, const char *desc);
  virtual ~TSolGEMChamber();

  Int_t ReadDatabase (const TDatime& date);

  Int_t Decode( const THaEvData & );

  // Return positions of chamber edges, in its own reference frame
  Double_t GetLowerEdgeX() const {return (GetOrigin())[0] - (GetSize())[0];}
  Double_t GetLowerEdgeY() const {return (GetOrigin())[1] - (GetSize())[1];}
  Double_t GetUpperEdgeX() const {return (GetOrigin())[0] + (GetSize())[0];}
  Double_t GetUpperEdgeY() const {return (GetOrigin())[1] + (GetSize())[1];}

  UInt_t GetNPlanes() const {return fNPlanes;};
  Double_t GetAngle() const {return fAngle;};

  TSolGEMPlane& GetPlane (UInt_t i) const {return *(fPlanes[i]);};
  void InitPlane (const UInt_t i, const char* name, const char* desc);
  void Print (const Bool_t printplanes = kTRUE);

 private:
  TSolGEMPlane** fPlanes; // List of chambers
  UInt_t fNPlanes;
  Double_t fAngle;        // Angle of orientation

 public:
  ClassDef(TSolGEMChamber,0)

    };

#endif//__TSOLGEMCHAMBER_H
