#ifndef __TSOLGEMCHAMBER_H
#define __TSOLGEMCHAMBER_H

#include "THaDetector.h"
#include "THaEvData.h"

#include "TSolGEMPlane.h"

class TSolGEMChamber : public THaDetector {
 public:
  TSolGEMChamber(const char *name, const char *desc);
  virtual ~TSolGEMChamber();

  Int_t ReadDatabase (const TDatime& date);

  Int_t Decode( const THaEvData & );

  Double_t GetLowerEdgeX() const {return (GetOrigin())[0] - (GetSize())[0];}
  Double_t GetLowerEdgeY() const {return (GetOrigin())[1] - (GetSize())[1];}
  Double_t GetUpperEdgeX() const {return (GetOrigin())[0] + (GetSize())[0];}
  Double_t GetUpperEdgeY() const {return (GetOrigin())[1] + (GetSize())[1];}
  UInt_t GetNPlanes() const {return fNPlanes;};
  TSolGEMPlane& GetPlane (UInt_t i) const {return *(fPlanes[i]);};
  void InitPlane (const UInt_t i, const char* name, const char* desc);
  void Print (const Bool_t printplanes = kTRUE);

 private:
  TSolGEMPlane** fPlanes; // List of chambers
  UInt_t fNPlanes;

 public:
  ClassDef(TSolGEMChamber,0)

    };

#endif//__TSOLGEMCHAMBER_H
