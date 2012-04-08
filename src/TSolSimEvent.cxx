//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSolSimEvent
//
//   Common class definitions for SoLID simulation decoder
//
/////////////////////////////////////////////////////////////////////

#include "TSolSimEvent.h"
#include "TClonesArray.h"
#include "TString.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
void TSolSimTrack::Print( const Option_t* opt ) const
{
  //TObject::Print(opt);
  cout << "track: type = " << fType << ", hit1 = " << fHit1 
       << ", PID = " << fPID << ", t0 = " << fTimeOffset << endl;
  cout << "  Origin    = (" << fOrigin.X() << ", " << fOrigin.Y() << ", " 
       << fOrigin.Z() << ")" << endl;
  cout << "  Momentum  = (" << fMomentum.Px() << ", " 
       << fMomentum.Py() << ", " << fMomentum.Pz() << ")" << endl;
}

//-----------------------------------------------------------------------------
TSolSimEvent::TSolSimEvent()
  : fRunID(0), fEvtID(0), fRefFile(0), fMCTracks(0), 
    fNClust(0), fNSignal(0), fNStrips(0)
{

}			     

//-----------------------------------------------------------------------------
TSolSimEvent::~TSolSimEvent()
{
  delete fMCTracks;
}

//-----------------------------------------------------------------------------
void TSolSimEvent::Clear( const Option_t* /*opt*/ )
{
  // Clear the event in preparation for reading next tree entry

  fRunID = fEvtID = fRefFile = 0;
  fNClust = fNSignal = fNStrips = 0;
  if( fMCTracks ) fMCTracks->Clear();
}

//-----------------------------------------------------------------------------
void TSolSimEvent::Print( const Option_t* opt ) const
{
  // Print current event info

  cout << ">>>>> =====================================" << endl;
  cout << "Event number:               " << fEvtID << endl;
  cout << "Number of fired GEM strips: " << fNStrips  << endl;

  TString sopt(opt);
  bool do_hit    = sopt.Contains("hit", TString::kIgnoreCase);
  // bool do_signal = sopt.Contains("signal", TString::kIgnoreCase);
  // bool do_backgr = sopt.Contains("backgr", TString::kIgnoreCase);
  // bool do_pileup = sopt.Contains("pileup", TString::kIgnoreCase);
  if( do_hit ) {
    for( Int_t i=0; i<fNStrips; ++i ) {
      cout << "hit = " << i
	   << ", GEM = "    << fStpGEM[i]
	   << ", plane = "  << fStpPlane[i]
	   << ", strip = "  << fStpNum[i]
	   << ", type = "   << fStpSigType[i]
	   << ", adc = ";
      for( int k=0; k<MC_MAXSAMP; k++ ) {
	cout << fStpADC[k][i];
	if( k+1 != MC_MAXSAMP ) cout << ", ";
      }
      cout << endl;
    }
  }

  if( sopt.Contains("track", TString::kIgnoreCase) && fMCTracks ) {
    for( Int_t i=0; i<fMCTracks->GetLast()+1; ++i ) {
      fMCTracks->UncheckedAt(i)->Print(opt);
    }
  }
}

//-----------------------------------------------------------------------------
ClassImp(TSolSimEvent)
ClassImp(TSolSimTrack)
