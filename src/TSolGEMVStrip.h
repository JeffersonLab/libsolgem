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

  // TODO: rewrite this using STL containers ...
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
  void SetTime(Double_t val) {fTime = val;};
  void SetHitCharge(Float_t val) { fHitCharge = val; }; // total charge of the avalanche
  void SetSize(Short_t n) { fSize = n; };

  Short_t GetSize() const { return fSize; };
  Short_t GetIdx(Short_t n) const { return fIdx->At(n); }; 
  Short_t GetADC(Short_t n, Short_t sample) const { return fADC->At(n*fNsample+sample); }; 
  Double_t GetCharge(Short_t n) const { return fCharge->At(n); };
  Float_t GetTime() const { return fTime; };
  Double_t GetHitCharge() const { return fHitCharge; };

  void Print();
};

#endif
