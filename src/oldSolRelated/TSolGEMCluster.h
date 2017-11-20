#ifndef __TSOLGEMCLUSTER_
#define __TSOLGEMCLUSTER_

#include "TROOT.h"

class TSolGEMCluster : public TObject {
public:
  TSolGEMCluster() {printf("I'm a cluster!\n");}
  virtual ~TSolGEMCluster() {;}
  
  Double_t GetPos() const { return fPos; }
  Double_t GetE()   const { return fE; }

private:
  Double_t fPos;
  Double_t fE;
  
public:
  ClassDef(TSolGEMCluster,1)

};
#endif//__TSOLGEMCLUSTER_
