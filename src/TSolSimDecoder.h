#ifndef __TSolSimDecoder_h
#define __TSolSimDecoder_

/////////////////////////////////////////////////////////////////////
//
//   TSolSimDecoder
//
/////////////////////////////////////////////////////////////////////

#include "THaEvData.h"
#include "THaAnalysisObject.h"
#include "TList.h"

class THaCrateMap;

class TSolSimDecoder : public THaEvData {
 public:
  TSolSimDecoder();
  virtual ~TSolSimDecoder();

  Int_t  LoadEvent( const int*evbuffer, THaCrateMap* usermap );

  void   Clear( Option_t* opt="" );
  Int_t  GetNTracks() const { return fTracks.GetSize(); }
  Int_t  DefineVariables( THaAnalysisObject::EMode mode = 
			  THaAnalysisObject::kDefine );

  static const Int_t CHAN_PER_SLOT = 128;

 protected:

  TList   fTracks;       // Monte Carlo tracks

  bool    fIsSetup;

  ClassDef(TSolSimDecoder,0) // Decoder for simulated SoLID spectrometer data
};

#endif
