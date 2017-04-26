//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSBSSimEvent
//
//   Common class definitions for SoLID simulation decoder
//
/////////////////////////////////////////////////////////////////////

#include "TSBSSimEvent.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
TSBSSimTrack::TSBSSimTrack( Int_t number, Int_t pid,
			    const TVector3& vertex, const TVector3& momentum )
  : Podd::MCTrack( number, pid, vertex, momentum )
{
}

// Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
// refer to comment in TSBSSimEvent.h l. 30-32
/*
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::MCFitR() const
{
  return TMath::Sqrt(fMCFitPar[0]*fMCFitPar[0] + fMCFitPar[2]*fMCFitPar[2] );
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::MCFitPhi() const
{
  return TVector3( fMCFitPar[0], fMCFitPar[2], 0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::MCFitThetaDir() const
{
  return TVector3( fMCFitPar[1], fMCFitPar[3], 1.0 ).Theta();
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::MCFitPhiDir() const
{
  return TVector3( fMCFitPar[1], fMCFitPar[3], 1.0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::RcFitR() const
{
  return TMath::Sqrt(fRcFitPar[0]*fRcFitPar[0] + fRcFitPar[2]*fRcFitPar[2] );
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::RcFitPhi() const
{
  return TVector3( fRcFitPar[0], fRcFitPar[2], 0 ).Phi();
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::RcFitThetaDir() const
{
  return TVector3( fRcFitPar[1], fRcFitPar[3], 1.0 ).Theta();
}
//-----------------------------------------------------------------------------
Double_t TSBSSimTrack::RcFitPhiDir() const
{
  return TVector3( fRcFitPar[1], fRcFitPar[3], 1.0 ).Phi();
}
*/

//-----------------------------------------------------------------------------
TSBSSimEvent::TSBSSimEvent()
  : fRunID(0), fEvtID(0), fMCTracks(0), fNSignal(0), fSectorsMapped(false),
    fSignalSector(0)
{
}

//-----------------------------------------------------------------------------
TSBSSimEvent::TSBSSimEvent( UInt_t ntracks )
  : fRunID(0), fEvtID(0), fNSignal(0), fSectorsMapped(false),
    fSignalSector(0)
{
  if( ntracks == 0 ) ntracks = 1;
  fMCTracks = new TClonesArray( "TSBSSimTrack", ntracks );
}

//-----------------------------------------------------------------------------
TSBSSimEvent::~TSBSSimEvent()
{
  delete fMCTracks;
}

//-----------------------------------------------------------------------------
TSBSSimTrack* TSBSSimEvent::AddTrack( Int_t number, Int_t pid,
				      const TVector3& vertex,
				      const TVector3& momentum )
{
  // Add a physics track with the given parameters

  return
    new( (*fMCTracks)[GetNtracks()] ) TSBSSimTrack( number, pid,
						    vertex, momentum );
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Clear( const Option_t* opt )
{
  // Clear the event in preparation for reading next tree entry

  TString sopt(opt);

  fNSignal = 0;
  fGEMClust.clear();
  fGEMStrips.clear();
  if( sopt.Contains("all",TString::kIgnoreCase) ) {
    if( fMCTracks ) fMCTracks->Clear(opt);
  }
}

//-----------------------------------------------------------------------------
Int_t TSBSSimEvent::GetNtracks() const
{
  // Get number of physics tracks

  return fMCTracks->GetLast()+1;
}

//-----------------------------------------------------------------------------
void TSBSSimEvent::Print( const Option_t* opt ) const
{
  // Print current event info

  cout << ">>>>> =====================================" << endl;
  cout << "Event number:               " << fEvtID << endl;
  cout << "Number of hits:             " << fGEMClust.size()   << endl;
  cout << "Number of fired GEM strips: " << fGEMStrips.size()  << endl;

  TString sopt(opt);
  bool do_all    = sopt.Contains("all",   TString::kIgnoreCase);
  bool do_hit    = sopt.Contains("hit",   TString::kIgnoreCase) || do_all;
  bool do_clust  = sopt.Contains("clust", TString::kIgnoreCase) || do_all;
  bool do_track  = sopt.Contains("track", TString::kIgnoreCase) || do_all;
  // bool do_signal = sopt.Contains("signal", TString::kIgnoreCase);
  // bool do_backgr = sopt.Contains("backgr", TString::kIgnoreCase);
  // bool do_pileup = sopt.Contains("pileup", TString::kIgnoreCase);

  if( do_track && fMCTracks ) {
    for( Int_t i=0; i<GetNtracks(); ++i ) {
      fMCTracks->UncheckedAt(i)->Print(opt);
    }
  }
  if( do_clust ) {
    for( vector<GEMCluster>::const_iterator ic = fGEMClust.begin();
	 ic != fGEMClust.end(); ++ic ) {
      const GEMCluster& c = *ic;
      cout << "hit = " << c.fID
	   << ", src = "    << c.fSource
	   << ", sect = "   << c.fRealSector;
      if( fSectorsMapped ) cout << "->" << c.fSector;
      cout << ", plane = "  << c.fPlane
	   << ", GID = "    << c.fType
	   << ", TRID = "   << c.fTRID
	   << ", PID = "    << c.fPID
	   << ", chrg = "   << c.fCharge
	   << ", time = "   << c.fTime
	   << ", size = ("  << c.fSize[0]  << ", " << c.fSize[1]  << ")"
	   << ", start = (" << c.fStart[0] << ", " << c.fStart[1] << ")"
	   << ", coord = (" << c.fXProj[0] << ", " << c.fXProj[1] << ")"
	   << endl;
      cout << " mom lab  = "; c.fP.Print();
      cout << " mom spec = "; c.fPspec.Print();
      cout << " mcpos  = "; c.fMCpos.Print();
      cout << " hitpos = "; c.fHitpos.Print();
    }
    if( do_hit && !fGEMClust.empty() )
      cout << "-------------" << endl;
  }

  if( do_hit ) {
    UInt_t i = 0;
    for( vector<DigiGEMStrip>::const_iterator is = fGEMStrips.begin();
	 is != fGEMStrips.end(); ++is ) {
      const DigiGEMStrip& s = *is;
      cout << "strip = " << i++
	   << ", sect = "   << s.fSector
	   << ", plane = "  << s.fPlane
	   << ", proj = "   << s.fProj
	   << ", strip = "  << s.fChan
	   << ", type = "   << s.fSigType
	   << ", chrg = "   << s.fCharge
	   << ", adc = ";
      for( int k=0; k<s.fNsamp; k++ ) {
	cout << s.fADC[k];
	if( k+1 != s.fNsamp ) cout << ", ";
      }
      cout << ", hits = ";
      for( Int_t isc = 0; isc < s.fClusters.GetSize(); isc++ ) {
	cout << s.fClusters[isc];
	if( isc+1 != s.fClusters.GetSize() ) cout << ", ";
      }
      cout << endl;
    }
  }
}

//-----------------------------------------------------------------------------
ClassImp(TSBSSimEvent)
ClassImp(TSBSSimTrack)
