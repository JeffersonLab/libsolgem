//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSBSSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as TSBSSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "TSBSSimDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "TSolDBManager.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

//EFuchey: 2016/12/10: it is necessary to declare the TSolDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
static TSolDBManager* fManager = TSolDBManager::GetInstance();
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane };

Double_t TSBSSimDecoder::fgCaloZ  = 6.8;
Double_t TSBSSimDecoder::fgCaloRes  = 0.01;
Bool_t   TSBSSimDecoder::fgDoCalo = false;
Double_t TSBSSimDecoder::fgZ0 = 3.435510;

typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
TSBSSimDecoder::TSBSSimDecoder()
{
  // Constructor

  fMCHits     = new TClonesArray( "TSBSSimGEMHit",    200 );
  fMCTracks   = new TClonesArray( "TSBSSimTrack",       1 );
  fBackTracks = new TClonesArray( "TSBSSimBackTrack",   5 );
  
  fgZ0 = fManager->GetZ0();
  fgDoCalo = fManager->GetDoCalo();
  fgCaloZ = fManager->GetCaloZ();
  fgCaloRes = fManager->GetCaloRes();
  
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
}

//-----------------------------------------------------------------------------
TSBSSimDecoder::~TSBSSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );

  delete fBackTracks;
  // fMCHits and fMCTracks are deleted by SimDecoder destructor
}

//-----------------------------------------------------------------------------
Int_t TSBSSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char* const here = "TSBSSimDecoder::DefineVariables";

  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;

  SimDecoder::DefineVariables( mode );

  RVarDef vars[] = {
    // Generated track info
    //    { "tr.n",      "Number of tracks",      "GetNMCTracks()" },
    { "tr.x",      "Track origin x (m)",    "fMCTracks.TSBSSimTrack.VX()" },
    { "tr.y",      "Track origin y (m)",    "fMCTracks.TSBSSimTrack.VY()" },
    { "tr.z",      "Track origin z (m)",    "fMCTracks.TSBSSimTrack.VZ()" },
    { "tr.p",      "Track momentum (GeV)",  "fMCTracks.TSBSSimTrack.P() "},
    { "tr.theta",  "Track theta_p (rad)",   "fMCTracks.TSBSSimTrack.PTheta()" },
    { "tr.phi",    "Track phi_p (rad)",     "fMCTracks.TSBSSimTrack.PPhi()" },
    { "tr.pid",    "Track PID (PDG)",       "fMCTracks.TSBSSimTrack.fPID" },
    { "tr.num",    "GEANT track number",    "fMCTracks.TSBSSimTrack.fNumber" },
    { "tr.planes", "Bitpattern of planes hit", "fMCTracks.TSBSSimTrack.fHitBits" },
    { "tr.nhits",  "Number of tracker hits","fMCTracks.TSBSSimTrack.fNHits" },
    { "tr.nfound", "Number of hits found",  "fMCTracks.TSBSSimTrack.fNHitsFound" },
    { "tr.flags",  "Reconstruction status", "fMCTracks.TSBSSimTrack.fReconFlags" },

    // Results of fit to MC points - measures multiple scattering
    { "tr.mcfit.r",     "Track x from MC fit [m]",
                                               "fMCTracks.TSBSSimTrack.MCFitR()" },
    { "tr.mcfit.phi",   "Track phi from MC fit [rad]",
                                             "fMCTracks.TSBSSimTrack.MCFitPhi()" },
    { "tr.mcfit.thdir", "Track dir theta from MC fit [rad]",
                                        "fMCTracks.TSBSSimTrack.MCFitThetaDir()" },
    { "tr.mcfit.phdir", "Track x from MC fit [rad]",
                                          "fMCTracks.TSBSSimTrack.MCFitPhiDir()" },
    { "tr.mcfit.chi2",  "Chi2 of MC fit",
                                           "fMCTracks.TSBSSimTrack.fMCFitPar[4]" },
    { "tr.mcfit.ndof",  "NDoF of MC fit",
                                           "fMCTracks.TSBSSimTrack.fMCFitPar[5]" },
    { "tr.mcfit.vx",    "Vertex x from MC fit [m]",
                                           "fMCTracks.TSBSSimTrack.fMCFitPar[6]" },
    { "tr.mcfit.vy",    "Vertex y from MC fit [m]",
                                           "fMCTracks.TSBSSimTrack.fMCFitPar[7]" },
    { "tr.mcfit.vz",    "Vertex z from MC fit [m]",
                                           "fMCTracks.TSBSSimTrack.fMCFitPar[8]" },

    // Results of fit to reconstructed MC hits - checks hit resolution effects
    // independent of track finding
    { "tr.fit.r",     "Track x from rec hit fit [m]",
                                               "fMCTracks.TSBSSimTrack.RcFitR()" },
    { "tr.fit.phi",   "Track phi from rec hit fit [rad]",
                                             "fMCTracks.TSBSSimTrack.RcFitPhi()" },
    { "tr.fit.thdir", "Track dir theta from rec hit fit [rad]",
                                        "fMCTracks.TSBSSimTrack.RcFitThetaDir()" },
    { "tr.fit.phdir", "Track x from rec hit fit [rad]",
                                          "fMCTracks.TSBSSimTrack.RcFitPhiDir()" },
    { "tr.fit.chi2",  "Chi2 of rec hit fit",
                                           "fMCTracks.TSBSSimTrack.fRcFitPar[4]" },
    { "tr.fit.ndof",  "NDoF of rec hit fit",
                                           "fMCTracks.TSBSSimTrack.fRcFitPar[5]" },
    { "tr.fit.vx",    "Vertex x from rec hit fit [m]",
                                           "fMCTracks.TSBSSimTrack.fRcFitPar[6]" },
    { "tr.fit.vy",    "Vertex y from rec hit fit [m]",
                                           "fMCTracks.TSBSSimTrack.fRcFitPar[7]" },
    { "tr.fit.vz",    "Vertex z from rec hit fit [m]",
                                           "fMCTracks.TSBSSimTrack.fRcFitPar[8]" },

    // "Back tracks": hits of the primary particle in the first tracker plane
    { "btr.n",     "Number of back tracks",     "GetNBackTracks()" },
    { "btr.pid",   "Track PID (PDG)",           "fBackTracks.TSBSSimBackTrack.fPID" },
    { "btr.num",   "GEANT particle number",     "fBackTracks.TSBSSimBackTrack.fType" },
    { "btr.planes","Bitpattern of planes hit",  "fBackTracks.TSBSSimBackTrack.fHitBits" },
    { "btr.ufail", "Undigitized u planes",      "fBackTracks.TSBSSimBackTrack.fUfailBits" },
    { "btr.vfail", "Undigitized v planes",      "fBackTracks.TSBSSimBackTrack.fVfailBits" },
    { "btr.sect",  "Sector number",             "fBackTracks.TSBSSimBackTrack.fSector" },
    { "btr.p",     "Track momentum (GeV)",      "fBackTracks.TSBSSimBackTrack.P() "},
    // Track position in Cartesian/TRANSPORT coordinates, not optimal for SoLID
    { "btr.x",     "Track pos lab x [m]",       "fBackTracks.TSBSSimBackTrack.X()" },
    { "btr.y",     "Track pos lab y [m]",       "fBackTracks.TSBSSimBackTrack.Y()" },
    { "btr.th",    "Track dir tan(theta)",      "fBackTracks.TSBSSimBackTrack.ThetaT()" },
    { "btr.ph",    "Track dir tan(phi)",        "fBackTracks.TSBSSimBackTrack.PhiT()" },
    // Track position and direction in cylindrical coordinates, good for SoLID
    { "btr.r",     "Track pos lab r_trans (m)", "fBackTracks.TSBSSimBackTrack.R()" },
    { "btr.theta", "Track pos lab theta [rad]", "fBackTracks.TSBSSimBackTrack.Theta()" },
    { "btr.phi",   "Track pos lab phi [rad]",   "fBackTracks.TSBSSimBackTrack.Phi()" },
    { "btr.thdir", "Track dir theta [rad]",     "fBackTracks.TSBSSimBackTrack.ThetaDir()" },
    { "btr.phdir", "Track dir phi [rad]",       "fBackTracks.TSBSSimBackTrack.PhiDir()" },
    // Hit coordinates in first tracker plane, relative to plane origin
    { "btr.hx",    "Track pos plane x [m]",     "fBackTracks.TSBSSimBackTrack.HX()" },
    { "btr.hy",    "Track pos plane y [m]",     "fBackTracks.TSBSSimBackTrack.HY()" },

    // Digitized hits registered in the GEMs
    //    { "hit.n",     "Number of MC hits",          "GetNMCHits()" },
    { "hit.id",    "MC hit number",              "fMCHits.TSBSSimGEMHit.fID" },
    { "hit.sect",  "MC hit sector",              "fMCHits.TSBSSimGEMHit.fSector" },
    { "hit.rsect", "MC hit non-mapped sector",   "fMCHits.TSBSSimGEMHit.fRealSector" },
    { "hit.plane", "MC hit plane",               "fMCHits.TSBSSimGEMHit.fPlane" },
    { "hit.src",   "MC data set source",         "fMCHits.TSBSSimGEMHit.fSource" },
    { "hit.type",  "MC hit GEANT counter",       "fMCHits.TSBSSimGEMHit.fType" },
    { "hit.pid",   "MC hit PID (PDG)",           "fMCHits.TSBSSimGEMHit.fPID" },
    { "hit.p",     "MC hit particle mom [GeV]",  "fMCHits.TSBSSimGEMHit.P()" },
    { "hit.x",     "MC hit lab x position [m]",  "fMCHits.TSBSSimGEMHit.X()" },
    { "hit.y",     "MC hit lab y position [m]",  "fMCHits.TSBSSimGEMHit.Y()" },
    { "hit.z",     "MC hit lab z position [m]",  "fMCHits.TSBSSimGEMHit.Z()" },
    // Hit position in cylindrical/spherical coordinates, good for SoLID
    { "hit.r",     "MC hit lab r [m]",           "fMCHits.TSBSSimGEMHit.R()" },
    { "hit.theta", "MC hit lab theta [rad]",     "fMCHits.TSBSSimGEMHit.Theta()" },
    { "hit.phi",   "MC hit lab phi [rad]",       "fMCHits.TSBSSimGEMHit.Phi()" },
    { "hit.charge","MC hit cluster charge",      "fMCHits.TSBSSimGEMHit.fCharge" },
    { "hit.time",  "MC hit time offset [s]",     "fMCHits.TSBSSimGEMHit.fTime" },
    { "hit.usz",   "MC hit u cluster size",      "fMCHits.TSBSSimGEMHit.fUSize" },
    { "hit.ustart","MC hit u cluster 1st strip", "fMCHits.TSBSSimGEMHit.fUStart" },
    { "hit.upos",  "MC hit u cluster center [m]","fMCHits.TSBSSimGEMHit.fUPos" },
    { "hit.vsz",   "MC hit v cluster size",      "fMCHits.TSBSSimGEMHit.fVSize" },
    { "hit.vstart","MC hit v cluster 1st strip", "fMCHits.TSBSSimGEMHit.fVStart" },
    { "hit.vpos",  "MC hit v cluster center [m]","fMCHits.TSBSSimGEMHit.fVPos" },

    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void TSBSSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCHits, fMCTracks and fMCPoints

  fBackTracks->Clear(opt);
  fStripMap.clear();
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int TSBSSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int TSBSSimDecoder::LoadEvent(const Int_t* evbuffer )
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

//-----------------------------------------------------------------------------
static inline
void StripToROC( Int_t s_plane, Int_t s_sector, Int_t s_proj,
		 Int_t s_chan,
		 Int_t& crate, Int_t& slot, Int_t& chan )
{
  // Convert location parameters (plane,sector,proj,chan) of the given strip
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx

  div_t d = div( s_chan, fManager->GetChanPerSlot() );
  Int_t module = d.quot;
  chan = d.rem;
  Int_t ix = module +
    fManager->GetModulesPerReadOut()*( s_proj + fManager->GetNReadOut()*( s_plane + fManager->GetNTracker()*s_sector ));
  d = div( ix, fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan +
    fManager->GetChanPerSlot()*( slot + fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber()*crate );
}

//-----------------------------------------------------------------------------
Int_t TSBSSimDecoder::StripFromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fStripMap.empty() )
    return -1;

  StripMap_t::const_iterator found = fStripMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fStripMap.end() )
    return -1;

  return found->second;
}

//-----------------------------------------------------------------------------
MCHitInfo TSBSSimDecoder::GetMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel

  const char* const here = "TSBSSimDecoder::GetMCHitInfo";

  Int_t istrip = StripFromROC( crate, slot, chan );
  assert( istrip >= 0 );  // else logic error in caller or bad fStripMap

  assert( buffer );       // Must still have the event buffer
  const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);

  assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  const TSBSSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  assert( strip.fProj >= 0 && strip.fProj < fManager->GetNReadOut() );

  MCHitInfo mc;
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
	Warning( Here(here), "Event %d: Multiple hits of primary particle "
		 "in plane %d\nShould never happen. Call expert.",
		 simEvent->fEvtID, strip.fPlane );
	continue;
      }
      // Strip contains a contribution from a primary particle hit :)
      mc.fMCTrack = 1;    // currently only ever one primary particle per event
      mc.fMCPos   = c.fXProj[strip.fProj];
      mc.fMCTime  = c.fTime;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj];
      }
    }
  }
  assert( strip.fClusters.GetSize() == 0 || mc.fMCTrack > 0 || mc.fContam > 0 );

  if( mc.fMCTrack == 0 ) {
    if( mc.fContam > 1 ) {
      // If only background hits, report the mean position of all those hits
      mc.fMCPos /= static_cast<Double_t>(mc.fContam);
    }
    mc.fMCTime = strip.fTime1;
  }
  return mc;
}

//-----------------------------------------------------------------------------
static inline Int_t NumberOfSetBits( UInt_t v )
{
  // Count number of bits set in 32-bit integer. From
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
Int_t TSBSSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t TSBSSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  static const char* const here = "TSBSSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSBSSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);

  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    if( (ret = init_cmap()) != HED_OK )
      return ret;
    if( (ret = init_slotdata(fMap)) != HED_OK)
      return ret;
    first_decode = false;
  }
  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");

  // FIXED? needed?
  // evscaler = 0;
  // event_length = 0;
  
  event_type = 1;
  event_num = simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  fWeight = simEvent->fWeight;

  if( fDoBench ) fBench->Begin("physics_decode");

  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<TSBSSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    const TSBSSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    Int_t crate, slot, chan;
    StripToROC( s.fPlane, s.fSector, s.fProj, s.fChan, crate, slot, chan );
    for( Int_t k = 0; k < s.fNsamp; k++ ) {
      Int_t raw = s.fADC[k];
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
    // Build map from ROC address to strip index. This is needed to extract
    // the MC truth info later in the tracking detector decoder via GetMCChanInfo.
#ifndef NDEBUG
    pair<StripMap_t::const_iterator,bool> ins =
#endif
      fStripMap.insert( make_pair( MakeROCKey(crate,slot,chan), i ) );
    assert( ins.second );
  }

  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks. We need to copy them here so we can export them as global
  // variables.
  TClonesArray* tracks = simEvent->fMCTracks;
  assert( tracks );
  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    TSBSSimTrack* trk = static_cast<TSBSSimTrack*>(tracks->UncheckedAt(i));
    new( (*fMCTracks)[i] ) TSBSSimTrack(*trk);
  }
  assert( GetNMCTracks() > 0 );

  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fManager->GetNTracker(), primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<TSBSSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fManager->GetNTracker() ) {
      Error( here, "Illegal plane number = %d in cluster. "
	     "Should never happen. Call expert.", c.fPlane );
      simEvent->Print("clust");
      return HED_FATAL;
    }

    // Save hits in the GEMs
    new( (*fMCHits)[GetNMCHits()] ) TSBSSimGEMHit(c);

    // Extra bookkeeping for primary tracks, used for making back tracks below
    if( c.fType == kPrimaryType && c.fSource == kPrimarySource ) {
      // Record the primary track's points for access via the SimDecoder interface.
      // Record one point per projection so that we can study residuals.
      Int_t itrack = 1;
      primary_sector = c.fSector;
      MCTrackPoint* upt =
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kUPlane,
							  c.fMCpos, c.fP );
      upt->fMCTime = c.fTime;
      MCTrackPoint* vpt =
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kVPlane,
							  c.fMCpos, c.fP );
      vpt->fMCTime = c.fTime;

      // Keep bitpattern of planes crossed by this primary
      SETBIT(primary_hitbits,c.fPlane);
      // Save index of the primary particle hit closest to plane 0
      if( c.fPlane < best_primary_plane ) {
	best_primary = i;
	best_primary_plane = c.fPlane;
      }
      // Determine digitization hit inefficiency: Check if this MC hit
      // activated GEM strips in both readout planes
      if( c.fSize[0] == 0 ) {
	SETBIT(ufail, c.fPlane);
	CLRBIT(upt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(upt->fStatus, MCTrackPoint::kDigitized);
      }
      if( c.fSize[1] == 0 ) {
	SETBIT(vfail, c.fPlane);
	CLRBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      } else {
	SETBIT(vpt->fStatus, MCTrackPoint::kDigitized);
      }
    }
  }

  // Sort fMCPoints by type (u,v) and plane number, then calculate plane-to-plane
  // differences. The following assumes that all points are from the same track
  // (ensured above). If that is no longer so one day, fMCPoints will need to
  // be sorted by track number as well, and the algo below needs to be changed.
  fMCPoints->Sort();
  Double_t mass = 0;
  TSBSSimTrack* trk = static_cast<TSBSSimTrack*>(fMCTracks->UncheckedAt(0));
  assert(trk);
  if( TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(trk->fPID) )
    mass = particle->Mass();
  else
    Warning( "LoadEvent", "No enrty in PDG database for PID = %d", trk->fPID );

  MCTrackPoint* prev_pt = 0;
  for( Int_t i = 0; i < GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(i) );
    assert(pt);
    if( prev_pt && prev_pt->fType == pt->fType ) {
      assert( pt->fMCTrack == prev_pt->fMCTrack );
      if( prev_pt->fPlane+1 == pt->fPlane ) {
	pt->fDeltaE = TMath::Sqrt(prev_pt->fMCP.Mag2() + mass*mass) -
	  TMath::Sqrt(pt->fMCP.Mag2() + mass*mass);
	pt->fDeflect = prev_pt->fMCP.Angle(pt->fMCP);
	pt->fToF = pt->fMCTime - prev_pt->fMCTime;
      }
    }
    prev_pt = pt;
  }

  // Keep statistics in the MC track
  trk->fNHits = 2*NumberOfSetBits(primary_hitbits);
  trk->fHitBits = primary_hitbits;

  // "Back tracks"
  // Record the apparent track from the primary particle
  // of the signal data here, i.e. type == 1 and source == 0.
  // There is only ever one primary particle per event.
  if( best_primary >= 0 ) {

    Int_t nback = GetNBackTracks();
    assert( nback == 0 );

    TSBSSimBackTrack* btr = new( (*fBackTracks)[nback] )
      TSBSSimBackTrack(simEvent->fGEMClust[best_primary]);
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fManager->GetNTracker() GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    if( fgDoCalo && trk->fNHits == 2*fManager->GetNTracker() ) {
      // Retrieve last MC track point
      assert( GetNMCPoints() == 2*fManager->GetNTracker() );
      MCTrackPoint* pt =
	static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fManager->GetNTracker()-1) );
      assert( pt );
      const TVector3& pos = pt->fMCPoint;
      TVector3 dir = pt->fMCP.Unit();
      if( fgCaloZ <= pos.Z() ) {
	Error( here, "Calorimeter z = %lf less than z of last GEM plane = %lf. "
	       "Set correct value with SetCaloZ() or turn off calo emulation.",
	       fgCaloZ, pos.Z() );
	return HED_FATAL;
      }
      if( TMath::Abs(dir.Z()) < 1e-6 ) {
	Error( here, "Illegal primary track direction (%lf,%lf,%lf). "
	       "Should never happen. Call expert.", dir.X(), dir.Y(), dir.Z() );
	return HED_ERR;
      }
      dir *= 1.0/dir.Z();  // Make dir a transport vector
      TVector3 hitpos = pos + (fgCaloZ-pos.Z()) * dir;

      // Smear the position with the given resolution
      // Assumes z-axis normal to calorimeter plane. Otherwise we would have to
      // get the plane's fXax and fYax
      TVector3 res( gRandom->Gaus(0.0, fgCaloRes),
		    gRandom->Gaus(0.0, fgCaloRes), 0.0 );
      hitpos += res;

      // Encode the raw hit data for the dummy GEM planes.
      // The actual coordinate transformation to u or v takes place in each
      // plane's Decode() where all the required geometry information is at hand.
      // This bypasses any type of digitization. That should be fine for dummy
      // planes where we want to inject known lab hit positions into the tracking.
      //
      // Because of the way the detector map is layed out at the moment,
      // we place the calorimeter in fake sector 31, so the data are in two slots
      // (for u and v coordinates, respectively) in the ROC immediately
      // following the GEM trackers for sector 30. In each slot, channels
      // 0-29 correspond to the sector of the MC track sector (should always be 0
      // if mapping sectors. Each "hit" corresponds to one measured position.
      // Currently, there is only ever one hit per channel since there is only
      // one MC track. The hit's raw data are hitpos.X(), the data, hitpos.Y(),
      // each a Float_t value interpreted as Int_t.
      assert( primary_sector == 0 );

      union FloatIntUnion {
	Float_t f;
	Int_t   i;
      } datx, daty;
      datx.f = static_cast<Float_t>(hitpos.X());
      daty.f = static_cast<Float_t>(hitpos.Y());

      Int_t crate, slot, chan;
      StripToROC( 0, fManager->GetNSector(), kUPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
      StripToROC( 0, fManager->GetNSector(), kVPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
    }
  }

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNMCTracks() << endl;
  //fMCTracks.Print();

  return HED_OK;
}

//-----------------------------------------------------------------------------
TSBSSimGEMHit::TSBSSimGEMHit( const TSBSSimEvent::GEMCluster& c )
  : fID(c.fID), fSector(c.fSector), fPlane(c.fPlane),
    fRealSector(c.fRealSector), fSource(c.fSource), fType(c.fType),
    fPID(c.fPID), fP(c.fP), fXEntry(c.fXEntry), fMCpos(c.fMCpos),
    fHitpos(c.fHitpos), fCharge(c.fCharge), fTime(c.fTime),
    fUSize(c.fSize[0]), fUStart(c.fStart[0]), fUPos(c.fXProj[0]),
    fVSize(c.fSize[1]), fVStart(c.fStart[1]), fVPos(c.fXProj[1])
{
  // Construct hit from cluster
}

//-----------------------------------------------------------------------------
void TSBSSimGEMHit::Print( const Option_t* ) const
{
  // Print TSBSSimGEMHit info

}

//-----------------------------------------------------------------------------
TSBSSimBackTrack::TSBSSimBackTrack( const TSBSSimEvent::GEMCluster& c )
  : fType(c.fType), fPID(c.fPID), fSector(c.fSector), fSource(c.fSource),
    fHitBits(0), fUfailBits(0), fVfailBits(0)
{
  // Construct track from cluster info

  Update( c );
  SetHitBit( c.fPlane );
}

//-----------------------------------------------------------------------------
Int_t TSBSSimBackTrack::Update( const TSBSSimEvent::GEMCluster& c )
{
  // Project track coordinates to first tracker plane

  static const char* const here = "TSBSSimBackTrack::Update";

  // Currently not needed since Update only called from constructor
  // if( fType != c.fType || fPID != c.fPID || fSector != c.fSector ) {
  //   Error( here, "Updating with inconsistent GEMCluster data: "
  // 	   "type = %d/%d, pid = %d/%d, sector = %d/%d.\n"
  // 	   "Should never happen. Call expert.",
  // 	   fType, c.fType, fPID, c.fPID, fSector, c.fSector );
  //   return -1;
  // }

  if( c.fPlane > 0 ) {
    Double_t dz = c.fMCpos.Z() - TSBSSimDecoder::GetZ0();
    if( dz <= 0 ) {
      Error( here, "Illegal fMCpos z-coordinate in plane = %d. "
	     "Should never happen. Call expert.", c.fPlane );
      c.fMCpos.Print();
      return -2;
    }
    fOrigin = c.fMCpos - dz * c.fP.Unit();
  } else {
    fOrigin = c.fMCpos;
  }
  fHitpos = c.fHitpos; // FIXME: project this, too?
  fMomentum = c.fP;

  return 0;
}

//-----------------------------------------------------------------------------
void TSBSSimBackTrack::Print( const Option_t* ) const
{
  // Print TSBSSimBackTrack info

  cout << "track: type = " << fType
       << ", PID = "       << fPID
       << ", sector = "    << fSector
       << endl;
  cout << "  Origin    = ";  fOrigin.Print();
  cout << "  Momentum  = ";  fMomentum.Print();
  cout << "  Hitpos    = ";  fHitpos.Print();
  cout << "  HitBits   = " << fHitBits << endl;
}

//-----------------------------------------------------------------------------
ClassImp(TSBSSimDecoder)
ClassImp(TSBSSimGEMHit)
ClassImp(TSBSSimBackTrack)
