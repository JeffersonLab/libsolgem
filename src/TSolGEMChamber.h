#ifndef __TSOLGEMCHAMBER_H
#define __TSOLGEMCHAMBER_H

#include "THaDetector.h"
#include "TSolWedge.h"

#include "THaEvData.h"
#include "TSolGEMPlane.h"

class TSolWedge;

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
