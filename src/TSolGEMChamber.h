#ifndef __TSOLGEMCHAMBER_H
#define __TSOLGEMCHAMBER_H

#include "THaDetector.h"

#include "types.h"

class TSolGEMPlane;
class TClonesArray;

class TSolGEMChamber : public THaDetector {
 public:
  TSolGEMChamber(const char *name, const char *desc);
  virtual ~TSolGEMChamber() {;}

  Int_t Decode( const THaEvData &);

  Double_t GetLowerEdgeX() const {return (GetOrigin())[0] - (GetSize())[0];}
  Double_t GetLowerEdgeY() const {return (GetOrigin())[1] - (GetSize())[1];}
  Double_t GetUpperEdgeX() const {return (GetOrigin())[0] + (GetSize())[0];}
  Double_t GetUpperEdgeY() const {return (GetOrigin())[1] + (GetSize())[1];}
  TSolGEMPlane& GetPlane (UInt_t i) const {return *((TSolGEMPlane*)((*fPlanes)[i]));};

 private:
  TClonesArray *fPlanes; // List of chambers

 public:
  ClassDef(TSolGEMChamber,1)

    };

#endif//__TSOLGEMCHAMBER_H
