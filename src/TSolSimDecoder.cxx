//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSolSimDecoder
//
//   Decoder for SoLID simulation data
//
//   Interprets event buffer from input as TSolSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "TSolSimDecoder.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"

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

// Constants for conversion of strip (plane,sector,proj,chan) to (crate,slot,chan)
static Int_t fNPLANES;
static Int_t fNSECTORS;
static Int_t fNPROJ;
static Int_t fCHAN_PER_SLOT;
static Int_t fmodules_per_readout;
static Int_t fmodules_per_chamber;
static Int_t fchambers_per_crate;

static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane };

Double_t TSolSimDecoder::fgZ0 = 1.571913;

Bool_t   TSolSimDecoder::fgDoCalo = false;
Double_t TSolSimDecoder::fgCaloZ  = 0.32;
Double_t TSolSimDecoder::fgCaloRes  = 0.01;

typedef vector<int>::size_type vsiz_t;

//-----------------------------------------------------------------------------
TSolSimDecoder::TSolSimDecoder(const char* filedbpath)
{
  // Constructor

  fMCHits     = new TClonesArray( "TSolSimGEMHit",    200 );
  fMCTracks   = new TClonesArray( "TSolSimTrack",       1 );
  fBackTracks = new TClonesArray( "TSolSimBackTrack",   5 );
  
  InitGeomParam(filedbpath);
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
}

//-----------------------------------------------------------------------------
TSolSimDecoder::~TSolSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );

  delete fBackTracks;
  // fMCHits and fMCTracks are deleted by SimDecoder destructor
}

//-----------------------------------------------------------------------------
// Reading database for miscellaneous parameters.
// This is done without proper DB request, but there is a set of default parameters.
// If those have to be used, user will be warned by a warning message.
// Data should be sorted as in file db/db_g4sbsmiscdata.dat
// This is the user's responsibility to make sure his input file is read correctly.
void TSolSimDecoder::InitGeomParam(const char* dbpath) {
  ifstream in(dbpath);
  if(!in.is_open()){
    printf("may not read geometry database at %s\n", dbpath);
    printf("using solid default params");
    
    fNPLANES = 5;
    fNSECTORS = 30;
    fNPROJ = 2;
    fCHAN_PER_SLOT = 1500;
    fmodules_per_readout = 1;
    fmodules_per_chamber = fNPROJ*fmodules_per_readout;
    fchambers_per_crate = (GetMAXSLOT()/fmodules_per_chamber/fNPLANES)*fNPLANES;
    fgZ0 = 1.571913;
    fgDoCalo = false;
    fgCaloZ  = 0.32;
    fgCaloRes  = 0.01;
  }else{
    cout << "Info: reading database at location " << dbpath << endl;
    cout <<" This file should be written the same way db/db_g4sbsmiscdata.dat "<< endl;
    cout << "(same structure, same order of parameters)" << endl;
    Float_t dummy;
    in.ignore(100,'=');
    in >> fNSECTORS;
    in.ignore(100,'=');
    in >> fNPLANES;
    in.ignore(100,'=');
    in >> fNPROJ;
    in.ignore(100,'=');
    in >> fCHAN_PER_SLOT;
    in.ignore(100,':');
    in >> fmodules_per_readout;
    fmodules_per_chamber = fNPROJ*fmodules_per_readout;
    fchambers_per_crate = (GetMAXSLOT()/fmodules_per_chamber/fNPLANES)*fNPLANES;

    in.ignore(100,'=');
    in >> fgZ0;
    in.ignore(100,'=');
    in >> fgDoCalo;
    in.ignore(100,'=');
    in >> fgCaloZ;
    in.ignore(100,'=');
    in >> fgCaloRes;
    
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
  }
  
  // cout << fNSECTORS << " " << fNPLANES << " " << fNPROJ << " " << fCHAN_PER_SLOT << " "
  // 	 << fmodules_per_readout << "; " << endl;
  //   << fgZ0 << " " << fgDoCalo << " " << fgCaloZ << " " << fgCaloRes << endl;
  
  in.close();
}

//-----------------------------------------------------------------------------
Int_t TSolSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char* const here = "TSolSimDecoder::DefineVariables";

  if( mode == THaAnalysisObject::kDefine && fIsSetup )
    return THaAnalysisObject::kOK;

  SimDecoder::DefineVariables( mode );

  RVarDef vars[] = {
    // Generated track info
    //    { "tr.n",      "Number of tracks",      "GetNMCTracks()" },
    { "tr.x",      "Track origin x (m)",    "fMCTracks.TSolSimTrack.VX()" },
    { "tr.y",      "Track origin y (m)",    "fMCTracks.TSolSimTrack.VY()" },
    { "tr.z",      "Track origin z (m)",    "fMCTracks.TSolSimTrack.VZ()" },
    { "tr.p",      "Track momentum (GeV)",  "fMCTracks.TSolSimTrack.P() "},
    { "tr.theta",  "Track theta_p (rad)",   "fMCTracks.TSolSimTrack.PTheta()" },
    { "tr.phi",    "Track phi_p (rad)",     "fMCTracks.TSolSimTrack.PPhi()" },
    { "tr.pid",    "Track PID (PDG)",       "fMCTracks.TSolSimTrack.fPID" },
    { "tr.num",    "GEANT track number",    "fMCTracks.TSolSimTrack.fNumber" },
    { "tr.planes", "Bitpattern of planes hit", "fMCTracks.TSolSimTrack.fHitBits" },
    { "tr.nhits",  "Number of tracker hits","fMCTracks.TSolSimTrack.fNHits" },
    { "tr.nfound", "Number of hits found",  "fMCTracks.TSolSimTrack.fNHitsFound" },
    { "tr.flags",  "Reconstruction status", "fMCTracks.TSolSimTrack.fReconFlags" },

    // Results of fit to MC points - measures multiple scattering
    { "tr.mcfit.r",     "Track x from MC fit [m]",
                                               "fMCTracks.TSolSimTrack.MCFitR()" },
    { "tr.mcfit.phi",   "Track phi from MC fit [rad]",
                                             "fMCTracks.TSolSimTrack.MCFitPhi()" },
    { "tr.mcfit.thdir", "Track dir theta from MC fit [rad]",
                                        "fMCTracks.TSolSimTrack.MCFitThetaDir()" },
    { "tr.mcfit.phdir", "Track x from MC fit [rad]",
                                          "fMCTracks.TSolSimTrack.MCFitPhiDir()" },
    { "tr.mcfit.chi2",  "Chi2 of MC fit",
                                           "fMCTracks.TSolSimTrack.fMCFitPar[4]" },
    { "tr.mcfit.ndof",  "NDoF of MC fit",
                                           "fMCTracks.TSolSimTrack.fMCFitPar[5]" },
    { "tr.mcfit.vx",    "Vertex x from MC fit [m]",
                                           "fMCTracks.TSolSimTrack.fMCFitPar[6]" },
    { "tr.mcfit.vy",    "Vertex y from MC fit [m]",
                                           "fMCTracks.TSolSimTrack.fMCFitPar[7]" },
    { "tr.mcfit.vz",    "Vertex z from MC fit [m]",
                                           "fMCTracks.TSolSimTrack.fMCFitPar[8]" },

    // Results of fit to reconstructed MC hits - checks hit resolution effects
    // independent of track finding
    { "tr.fit.r",     "Track x from rec hit fit [m]",
                                               "fMCTracks.TSolSimTrack.RcFitR()" },
    { "tr.fit.phi",   "Track phi from rec hit fit [rad]",
                                             "fMCTracks.TSolSimTrack.RcFitPhi()" },
    { "tr.fit.thdir", "Track dir theta from rec hit fit [rad]",
                                        "fMCTracks.TSolSimTrack.RcFitThetaDir()" },
    { "tr.fit.phdir", "Track x from rec hit fit [rad]",
                                          "fMCTracks.TSolSimTrack.RcFitPhiDir()" },
    { "tr.fit.chi2",  "Chi2 of rec hit fit",
                                           "fMCTracks.TSolSimTrack.fRcFitPar[4]" },
    { "tr.fit.ndof",  "NDoF of rec hit fit",
                                           "fMCTracks.TSolSimTrack.fRcFitPar[5]" },
    { "tr.fit.vx",    "Vertex x from rec hit fit [m]",
                                           "fMCTracks.TSolSimTrack.fRcFitPar[6]" },
    { "tr.fit.vy",    "Vertex y from rec hit fit [m]",
                                           "fMCTracks.TSolSimTrack.fRcFitPar[7]" },
    { "tr.fit.vz",    "Vertex z from rec hit fit [m]",
                                           "fMCTracks.TSolSimTrack.fRcFitPar[8]" },

    // "Back tracks": hits of the primary particle in the first tracker plane
    { "btr.n",     "Number of back tracks",     "GetNBackTracks()" },
    { "btr.pid",   "Track PID (PDG)",           "fBackTracks.TSolSimBackTrack.fPID" },
    { "btr.num",   "GEANT particle number",     "fBackTracks.TSolSimBackTrack.fType" },
    { "btr.planes","Bitpattern of planes hit",  "fBackTracks.TSolSimBackTrack.fHitBits" },
    { "btr.ufail", "Undigitized u planes",      "fBackTracks.TSolSimBackTrack.fUfailBits" },
    { "btr.vfail", "Undigitized v planes",      "fBackTracks.TSolSimBackTrack.fVfailBits" },
    { "btr.sect",  "Sector number",             "fBackTracks.TSolSimBackTrack.fSector" },
    { "btr.p",     "Track momentum (GeV)",      "fBackTracks.TSolSimBackTrack.P() "},
    // Track position in Cartesian/TRANSPORT coordinates, not optimal for SoLID
    { "btr.x",     "Track pos lab x [m]",       "fBackTracks.TSolSimBackTrack.X()" },
    { "btr.y",     "Track pos lab y [m]",       "fBackTracks.TSolSimBackTrack.Y()" },
    { "btr.th",    "Track dir tan(theta)",      "fBackTracks.TSolSimBackTrack.ThetaT()" },
    { "btr.ph",    "Track dir tan(phi)",        "fBackTracks.TSolSimBackTrack.PhiT()" },
    // Track position and direction in cylindrical coordinates, good for SoLID
    { "btr.r",     "Track pos lab r_trans (m)", "fBackTracks.TSolSimBackTrack.R()" },
    { "btr.theta", "Track pos lab theta [rad]", "fBackTracks.TSolSimBackTrack.Theta()" },
    { "btr.phi",   "Track pos lab phi [rad]",   "fBackTracks.TSolSimBackTrack.Phi()" },
    { "btr.thdir", "Track dir theta [rad]",     "fBackTracks.TSolSimBackTrack.ThetaDir()" },
    { "btr.phdir", "Track dir phi [rad]",       "fBackTracks.TSolSimBackTrack.PhiDir()" },
    // Hit coordinates in first tracker plane, relative to plane origin
    { "btr.hx",    "Track pos plane x [m]",     "fBackTracks.TSolSimBackTrack.HX()" },
    { "btr.hy",    "Track pos plane y [m]",     "fBackTracks.TSolSimBackTrack.HY()" },

    // Digitized hits registered in the GEMs
    //    { "hit.n",     "Number of MC hits",          "GetNMCHits()" },
    { "hit.id",    "MC hit number",              "fMCHits.TSolSimGEMHit.fID" },
    { "hit.sect",  "MC hit sector",              "fMCHits.TSolSimGEMHit.fSector" },
    { "hit.rsect", "MC hit non-mapped sector",   "fMCHits.TSolSimGEMHit.fRealSector" },
    { "hit.plane", "MC hit plane",               "fMCHits.TSolSimGEMHit.fPlane" },
    { "hit.src",   "MC data set source",         "fMCHits.TSolSimGEMHit.fSource" },
    { "hit.type",  "MC hit GEANT counter",       "fMCHits.TSolSimGEMHit.fType" },
    { "hit.pid",   "MC hit PID (PDG)",           "fMCHits.TSolSimGEMHit.fPID" },
    { "hit.p",     "MC hit particle mom [GeV]",  "fMCHits.TSolSimGEMHit.P()" },
    { "hit.x",     "MC hit lab x position [m]",  "fMCHits.TSolSimGEMHit.X()" },
    { "hit.y",     "MC hit lab y position [m]",  "fMCHits.TSolSimGEMHit.Y()" },
    { "hit.z",     "MC hit lab z position [m]",  "fMCHits.TSolSimGEMHit.Z()" },
    // Hit position in cylindrical/spherical coordinates, good for SoLID
    { "hit.r",     "MC hit lab r [m]",           "fMCHits.TSolSimGEMHit.R()" },
    { "hit.theta", "MC hit lab theta [rad]",     "fMCHits.TSolSimGEMHit.Theta()" },
    { "hit.phi",   "MC hit lab phi [rad]",       "fMCHits.TSolSimGEMHit.Phi()" },
    { "hit.charge","MC hit cluster charge",      "fMCHits.TSolSimGEMHit.fCharge" },
    { "hit.time",  "MC hit time offset [s]",     "fMCHits.TSolSimGEMHit.fTime" },
    { "hit.usz",   "MC hit u cluster size",      "fMCHits.TSolSimGEMHit.fUSize" },
    { "hit.ustart","MC hit u cluster 1st strip", "fMCHits.TSolSimGEMHit.fUStart" },
    { "hit.upos",  "MC hit u cluster center [m]","fMCHits.TSolSimGEMHit.fUPos" },
    { "hit.vsz",   "MC hit v cluster size",      "fMCHits.TSolSimGEMHit.fVSize" },
    { "hit.vstart","MC hit v cluster 1st strip", "fMCHits.TSolSimGEMHit.fVStart" },
    { "hit.vpos",  "MC hit v cluster center [m]","fMCHits.TSolSimGEMHit.fVPos" },

    { 0 }
  };

  return THaAnalysisObject::
    DefineVarsFromList( vars, THaAnalysisObject::kRVarDef,
			mode, "", this, MC_PREFIX, here );
}

//-----------------------------------------------------------------------------
void TSolSimDecoder::Clear( Option_t* opt )
{
  // Clear track and plane data

  SimDecoder::Clear(opt);   // clears fMCHits, fMCTracks and fMCPoints

  fBackTracks->Clear(opt);
  fStripMap.clear();
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
int TSolSimDecoder::LoadEvent(const UInt_t* evbuffer )
#else
int TSolSimDecoder::LoadEvent(const Int_t* evbuffer )
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

  div_t d = div( s_chan, fCHAN_PER_SLOT );
  Int_t module = d.quot;
  chan = d.rem;
  Int_t ix = module +
    fmodules_per_readout*( s_proj + fNPROJ*( s_plane + fNPLANES*s_sector ));
  d = div( ix, fchambers_per_crate*fmodules_per_chamber );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan +
    fCHAN_PER_SLOT*( slot + fchambers_per_crate*fmodules_per_chamber*crate );
}

//-----------------------------------------------------------------------------
Int_t TSolSimDecoder::StripFromROC( Int_t crate, Int_t slot, Int_t chan ) const
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
MCHitInfo TSolSimDecoder::GetMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel

  const char* const here = "TSolSimDecoder::GetMCHitInfo";

  Int_t istrip = StripFromROC( crate, slot, chan );
  assert( istrip >= 0 );  // else logic error in caller or bad fStripMap

  assert( buffer );       // Must still have the event buffer
  const TSolSimEvent* simEvent = reinterpret_cast<const TSolSimEvent*>(buffer);

  assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  const TSolSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  assert( strip.fProj >= 0 && strip.fProj < fNPROJ );

  MCHitInfo mc;
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1
    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TSolSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
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
Int_t TSolSimDecoder::DoLoadEvent(const UInt_t* evbuffer )
#else
Int_t TSolSimDecoder::DoLoadEvent(const Int_t* evbuffer )
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  static const char* const here = "TSolSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1,6,0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert( fMap || fNeedInit );

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSolSimFile). The pointer-to-unsigned integer is
  // needed compatibility with the standard decoder.
  const TSolSimEvent* simEvent = reinterpret_cast<const TSolSimEvent*>(buffer);

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

  // FIXED?: needed?
  //evscaler = 0;
  //event_length = 0;

  event_type = 1;
  event_num = simEvent->fEvtID;
  recent_event = event_num;

  // Event weight
  fWeight = simEvent->fWeight;

  if( fDoBench ) fBench->Begin("physics_decode");

  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<TSolSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    const TSolSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
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
    TSolSimTrack* trk = static_cast<TSolSimTrack*>(tracks->UncheckedAt(i));
    new( (*fMCTracks)[i] ) TSolSimTrack(*trk);
  }
  assert( GetNMCTracks() > 0 );

  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fNPLANES, primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<TSolSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TSolSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fNPLANES ) {
      Error( here, "Illegal plane number = %d in cluster. "
	     "Should never happen. Call expert.", c.fPlane );
      simEvent->Print("clust");
      return HED_FATAL;
    }

    // Save hits in the GEMs
    new( (*fMCHits)[GetNMCHits()] ) TSolSimGEMHit(c);

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
  TSolSimTrack* trk = static_cast<TSolSimTrack*>(fMCTracks->UncheckedAt(0));
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

    TSolSimBackTrack* btr = new( (*fBackTracks)[nback] )
      TSolSimBackTrack(simEvent->fGEMClust[best_primary]);
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fNPLANES GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    if( fgDoCalo && trk->fNHits == 2*fNPLANES ) {
      // Retrieve last MC track point
      assert( GetNMCPoints() == 2*fNPLANES );
      MCTrackPoint* pt =
	static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fNPLANES-1) );
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
      StripToROC( 0, fNSECTORS, kUPlane, primary_sector, crate, slot, chan );
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i)
	  == SD_ERR )
	return HED_ERR;
      StripToROC( 0, fNSECTORS, kVPlane, primary_sector, crate, slot, chan );
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
TSolSimGEMHit::TSolSimGEMHit( const TSolSimEvent::GEMCluster& c )
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
void TSolSimGEMHit::Print( const Option_t* ) const
{
  // Print TSolSimGEMHit info

}

//-----------------------------------------------------------------------------
TSolSimBackTrack::TSolSimBackTrack( const TSolSimEvent::GEMCluster& c )
  : fType(c.fType), fPID(c.fPID), fSector(c.fSector), fSource(c.fSource),
    fHitBits(0), fUfailBits(0), fVfailBits(0)
{
  // Construct track from cluster info

  Update( c );
  SetHitBit( c.fPlane );
}

//-----------------------------------------------------------------------------
Int_t TSolSimBackTrack::Update( const TSolSimEvent::GEMCluster& c )
{
  // Project track coordinates to first tracker plane

  static const char* const here = "TSolSimBackTrack::Update";

  // Currently not needed since Update only called from constructor
  // if( fType != c.fType || fPID != c.fPID || fSector != c.fSector ) {
  //   Error( here, "Updating with inconsistent GEMCluster data: "
  // 	   "type = %d/%d, pid = %d/%d, sector = %d/%d.\n"
  // 	   "Should never happen. Call expert.",
  // 	   fType, c.fType, fPID, c.fPID, fSector, c.fSector );
  //   return -1;
  // }

  if( c.fPlane > 0 ) {
    Double_t dz = c.fMCpos.Z() - TSolSimDecoder::GetZ0();
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
void TSolSimBackTrack::Print( const Option_t* ) const
{
  // Print TSolSimBackTrack info

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
ClassImp(TSolSimDecoder)
ClassImp(TSolSimGEMHit)
ClassImp(TSolSimBackTrack)
