#include <iostream>

#include "TSolGEMData.h"

using namespace std;

TSolGEMData::TSolGEMData (Int_t h)
{
  if (h > 0) InitHit (h);
}

TSolGEMData::~TSolGEMData()
{
}

void 
TSolGEMData::ClearHit() 
{
  if (fGem)
  {
    delete[] fIdxV;
    delete[] fGem;
    delete[] fEdep;
    delete[] fPID;
    
    for (Int_t k = 0; k < fNHit; k++) 
    {
      delete fXi[k];
      delete fXo[k];
      delete fXr[k];
      delete fMom[k];
    }
    delete[] fXi;
    delete[] fXo;
    delete[] fXr;
    delete[] fMom;
    
    cerr << __FUNCTION__ << "fGEM deleted " << endl;
    
    fGem = 0;
    fNHit = 0;    
  }
};

void
TSolGEMData::InitHit (Int_t h)
{
  if (h <= 0)
    return;

  fNHit = h; // number of hits in event

  // Hits

  fGem = new Int_t[h]; // chamber with hit
  fEdep = new Double_t[h]; // energy lost in drift
  fIdxV = new Int_t[h]; // entry number
  fPID = new Int_t[h]; // particle ID
  fXi = new TVector3*[h];
  fXo = new TVector3*[h];
  fXr = new TVector3*[h];
  fMom = new TVector3*[h];
  for (Int_t k = 0; k < fNHit; k++) 
  {
  	fXi[k] = new TVector3;
  	fXo[k] = new TVector3;
  	fXr[k] = new TVector3;
  	fMom[k] = new TVector3;
  }	  		

  cerr << __FUNCTION__ << "fGEM initialized with length " << h << endl;
}
