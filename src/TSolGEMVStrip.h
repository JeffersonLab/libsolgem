// Define the virtual strip class, used by digi
// N.B. The strip is virtual, not the class!
//

#ifndef __TSolGEMVStrip__
#define __TSolGEMVStrip__

#include <TArrayD.h>
#include <TArrayS.h>

class TSolGEMVStrip {

 private:
  // working variable of digi method (virtual strips)

  TArrayS *fIdx;  // index of the strip
  TArrayS *fADC; // value of the ADC sample
  TArrayD *fCharge; // total charge in strip (sampled)

  Float_t fTime; // time of first sampling
  Int_t fNsample; // number of samples

  Float_t fHitCharge;

  Short_t fSize; // effective size, maybe different from size of the array (some elements maybe not used)

 public:

  TSolGEMVStrip(Short_t n = 10, // number of strips in hit 
	       Short_t nsample = 10); 
  ~TSolGEMVStrip();

  void AddStripAt(Short_t idx, Short_t n) {fIdx->AddAt(idx,n);};
  void AddChargeAt(Double_t val, Short_t n) {fCharge->AddAt(val,n);};
  void AddSampleAt(Short_t val, Short_t sample, Short_t n) {fADC->AddAt(val,n*fNsample+sample);};
  void setTime(Double_t val) {fTime = val;};
  void setHitCharge(Float_t val) { fHitCharge = val; }; // total charge of the avalanche
  void setSize(Short_t n) { fSize = n; };

  Short_t getSize() { return fSize; };
  Short_t getIdx(Short_t n) { return fIdx->At(n); }; 
  Short_t getADC(Short_t n, Short_t sample) { return fADC->At(n*fNsample+sample); }; 
  Double_t getCharge(Short_t n) { return fCharge->At(n); };
  Float_t getTime() { return fTime; };
  Double_t getHitCharge() { return fHitCharge; };

  void Print();
};

#endif
