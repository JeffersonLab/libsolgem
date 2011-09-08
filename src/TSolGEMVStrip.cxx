#include <iostream>

#include "TSolGEMVStrip.h"

using namespace std;

TSolGEMVStrip::TSolGEMVStrip(Short_t n, // number of strips in hit 
			   Short_t nsample) {
  fSize = n;
  fNsample = nsample;
  fIdx = new TArrayS(n);
  fCharge = new TArrayD(n);
  fTime = -1.;
  fADC = new TArrayS(fNsample*n);
  fHitCharge = 0;
};

TSolGEMVStrip::~TSolGEMVStrip() {
  delete fIdx;
  delete fADC;
  delete fCharge;
};

void 
TSolGEMVStrip::Print() {
  Int_t i,k;
  
  cerr << "Virtual strips sampled starting at " << GetTime() << endl;
  
  for (i=0;i<GetSize();i++) {
    
    cerr << i << ") " << GetIdx(i) << " : ";
    cerr << GetCharge(i) << " = ";
    
    for (k=0;k<fNsample;k++) {
      cerr << GetADC(i,k) << " ";
    }
    
    cerr << endl;
  }
  
};
