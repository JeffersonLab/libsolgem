#ifndef __TSolSimAux__
#define __TSolSimAux__

// 
// Non member functions for simulation
//

#include <Rtypes.h>

namespace TSolSimAux
{
  Short_t ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits);
  Double_t PulseShape(Double_t t, 
		      Double_t C,  // normalization factor
		      Double_t Tp); // shaping time 
  Double_t PulseShape(Double_t t,
		      Double_t A,  // normalization
		      Double_t tau0,  // first time constant
		      Double_t tau1);  // second time constant
  Double_t Gaus2D(Double_t *x, Double_t *par);
  Double_t MultiGaus2D(Double_t *x, Double_t *par);
  Double_t SimpleCircle(Double_t *x, Double_t *par);
  Double_t MultiCircle(Double_t *x, Double_t *par);
}

#endif
