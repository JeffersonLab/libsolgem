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
#include "TClonesArray.h"
#include "TError.h"
#include <cstdlib>

using namespace std;

// Prefix of our own global variables (MC truth data)
static const char* const MC_PREFIX = "MC.";

//-----------------------------------------------------------------------------
TSolSimDecoder::TSolSimDecoder() : fIsSetup(false)
{
  // Constructor

  fHits = new TClonesArray( "TSolSimGEMHit", 200 );
  fBackTracks = new TClonesArray( "TSolSimBackTrack", 5 );

  DefineVariables();
}

//-----------------------------------------------------------------------------
TSolSimDecoder::~TSolSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );

  delete fBackTracks;
  delete fHits;
}

//-----------------------------------------------------------------------------
Int_t TSolSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  const char* const here = "TSolSimDecoder::DefineVariables";

  if( mode == THaAnalysisObject::kDefine && fIsSetup ) 
    return THaAnalysisObject::kOK;
  fIsSetup = ( mode == THaAnalysisObject::kDefine );

  RVarDef vars[] = {
    // Generated track info
    { "tr.n",      "Number of tracks",      "GetNTracks()" },
    { "tr.x",      "Track origin x (m)",    "fTracks.TSolSimTrack.VX()" },
    { "tr.y",      "Track origin y (m)",    "fTracks.TSolSimTrack.VY()" },
    { "tr.z",      "Track origin z (m)",    "fTracks.TSolSimTrack.VZ()" },
    { "tr.p",      "Track momentum (GeV)",  "fTracks.TSolSimTrack.P() "},
    { "tr.theta",  "Track theta_p (rad)",   "fTracks.TSolSimTrack.PTheta()" },
    { "tr.phi",    "Track phi_p (rad)",     "fTracks.TSolSimTrack.PPhi()" },
    { "tr.pid",    "Track PID (PDG)",       "fTracks.TSolSimTrack.fPID" },
    { "tr.num",    "GEANT track number",    "fTracks.TSolSimTrack.fNumber" },

    // "Back tracks": hits of the primary particle in the first tracker plane
    { "btr.n",     "Number of back tracks",     "GetNBackTracks()" },
    { "btr.pid",   "Track PID (PDG)",           "fBackTracks.TSolSimBackTrack.fPID" },
    { "btr.num",   "GEANT particle number",     "fBackTracks.TSolSimBackTrack.fType" },
    { "btr.planes","Bitpattern of planes hit",  "fBackTracks.TSolSimBackTrack.fHitBits" },
    { "btr.sect",  "Sector number",             "fBackTracks.TSolSimBackTrack.fSector" },
    { "btr.p",     "Track momentum (GeV)",      "fBackTracks.TSolSimBackTrack.P() "},
    // Track position in Cartesian/TRANSPORT coordinates, not optimal for SoLID
    { "btr.x",     "Track pos lab x (m)",       "fBackTracks.TSolSimBackTrack.X()" },
    { "btr.y",     "Track pos lab y (m)",       "fBackTracks.TSolSimBackTrack.Y()" },
    { "btr.th",    "Track dir tan(theta)",      "fBackTracks.TSolSimBackTrack.ThetaT()" },
    { "btr.ph",    "Track dir tan(phi)",        "fBackTracks.TSolSimBackTrack.PhiT()" },
    // Track position and direction in cylindrical coordinates, good for SoLID
    { "btr.r",     "Track pos lab r_trans (m)", "fBackTracks.TSolSimBackTrack.R()" },
    { "btr.theta", "Track pos lab theta (rad)", "fBackTracks.TSolSimBackTrack.Theta()" },
    { "btr.phi",   "Track pos lab phi (rad)",   "fBackTracks.TSolSimBackTrack.Phi()" },
    { "btr.thdir", "Track dir theta (rad)",     "fBackTracks.TSolSimBackTrack.ThetaDir()" },
    { "btr.phdir", "Track dir phi (rad)",       "fBackTracks.TSolSimBackTrack.PhiDir()" },
    // Hit coordinates in first tracker plane, relative to plane origin
    { "btr.hx",    "Track pos plane x (m)",     "fBackTracks.TSolSimBackTrack.HX()" },
    { "btr.hy",    "Track pos plane y (m)",     "fBackTracks.TSolSimBackTrack.HY()" },

    // All hits registered in the GEMs
    { "hit.n",     "Number of MC hits",          "GetNHits()" },
    { "hit.id",    "MC hit number",              "fHits.TSolSimGEMHit.fID" },
    { "hit.sect",  "MC hit sector",              "fHits.TSolSimGEMHit.fSector" },
    { "hit.plane", "MC hit plane",               "fHits.TSolSimGEMHit.fPlane" },
    { "hit.type",  "MC hit GEANT counter",       "fHits.TSolSimGEMHit.fType" },
    { "hit.pid",   "MC hit PID (PDG)",           "fHits.TSolSimGEMHit.fPID" },
    { "hit.p",     "MC hit particle p [GeV]",    "fHits.TSolSimGEMHit.P()" },
    { "hit.x",     "MC hit lab x position [m]",  "fHits.TSolSimGEMHit.X()" },
    { "hit.y",     "MC hit lab y position [m]",  "fHits.TSolSimGEMHit.Y()" },
    { "hit.z",     "MC hit lab z position [m]",  "fHits.TSolSimGEMHit.Z()" },
    { "hit.charge","MC hit cluster charge",      "fHits.TSolSimGEMHit.fCharge" },
    { "hit.time",  "MC hit time offset [s]",     "fHits.TSolSimGEMHit.fTime" },
    { "hit.usz",   "MC hit u cluster size",      "fHits.TSolSimGEMHit.fUSize" },
    { "hit.uwire", "MC hit u cluster 1st wire",  "fHits.TSolSimGEMHit.fUStart" },
    { "hit.upos",  "MC hit u cluster center [m]","fHits.TSolSimGEMHit.fUPos" },
    { "hit.vsz",   "MC hit v cluster size",      "fHits.TSolSimGEMHit.fVSize" },
    { "hit.vwire", "MC hit v cluster 1st wire",  "fHits.TSolSimGEMHit.fVStart" },
    { "hit.vpos",  "MC hit v cluster center [m]","fHits.TSolSimGEMHit.fVPos" },

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

  THaEvData::Clear();

  fHits->Clear();
  fBackTracks->Clear();
  // Never delete the physics tracks, only clear the list. The tracks are part
  // of a TClonesArray which is Clear()ed  in TSolSimFile::ReadEvent() by the
  // call to fEvent->Clear().
  fTracks.Clear("nodelete");
}

//-----------------------------------------------------------------------------
int TSolSimDecoder::LoadEvent(const int* evbuffer, THaCrateMap* map)
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  int ret = DoLoadEvent( evbuffer, map );

  if( fDoBench ) fBench->Stop("physics_decode");

  return ret;
}

//-----------------------------------------------------------------------------
int TSolSimDecoder::DoLoadEvent(const int* evbuffer, THaCrateMap* map)
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  static const char* const here = "TSolSimDecoder::LoadEvent";

  fMap = map;

  // Local copy of evbuffer pointer - any good use for it?
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSolSimFile). The pointer-to-integer is there
  // just for compatibility with the standard decoder.
  const TSolSimEvent* simEvent = reinterpret_cast<const TSolSimEvent*>(evbuffer);

  if (first_decode) {
    init_cmap();     
    if (init_slotdata(map) == HED_ERR) return HED_ERR;
    first_decode = false;
  }
  if( fDoBench ) fBench->Begin("clearEvent");
  Clear();
  for( int i=0; i<fNSlotClear; i++ )
    crateslot[fSlotClear[i]]->clearEvent();
  if( fDoBench ) fBench->Stop("clearEvent");
  
  // FIXME: needed?
  evscaler = 0;
  event_length = 0;

  event_type = 1;
  event_num = simEvent->fEvtID;
  recent_event = event_num;

  if( fDoBench ) fBench->Begin("physics_decode");

  // Decode the digitized strip data.  Populate crateslot array.
  const Int_t NPLANES = 4; // FIXME: get from database
  for( vector<TSolSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {

    const TSolSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    // The (roc,slot,chan) assignment must match the detmap definition.
    // See TreeSearch/dbconvert.cxx
    // FIXME: make parameters configurable
    const int NPROJ = 2, CHAN_PER_SLOT = 1280;
    const int modules_per_readout = 1;
    const int modules_per_chamber = 2*modules_per_readout;
    const int chambers_per_crate = (MAXSLOT/modules_per_chamber/NPLANES)*NPLANES;
    div_t d = div( s.fChan, CHAN_PER_SLOT );
    Int_t module = d.quot;
    Int_t chan   = d.rem;
    Int_t ix = module +
      modules_per_readout*( s.fProj + NPROJ*( s.fPlane + NPLANES*s.fSector ));
    d = div( ix, chambers_per_crate*modules_per_chamber );
    Int_t roc  = d.quot;
    Int_t slot = d.rem;

    for( Int_t k = 0; k < s.fNsamp; k++ ) {
      Int_t raw = s.fADC[k];
      if( crateslot[idx(roc,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
  }

  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks
  TClonesArray* tracks = simEvent->fMCTracks;
  if( tracks ) {
    for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
      fTracks.Add( tracks->UncheckedAt(i) );
    }
  }

  // MC hit data ("clusters") and "back tracks"
  for( vector<TSolSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TSolSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    new( (*fHits)[GetNHits()] ) TSolSimGEMHit(c);

    // Keep all hits of the primary track in the GEM planes
    if( c.fType == 1 ) {
      if( c.fPlane < 0 || c.fPlane >= NPLANES ) {
	Error( here, "Illegal plane number = %d in cluster"
	       "Should never happen. Call expert.", c.fPlane );
	simEvent->Print("clust");
	return HED_FATAL;
      }
      Int_t nback = GetNBackTracks();
      Int_t ib = 0;
      for( ; ib < nback; ++ib ) {
	TSolSimBackTrack* theTrack = GetBackTrack(ib);
	if( theTrack && theTrack->GetType() == c.fType ) {
	  if( theTrack->TestHitBit(c.fPlane) ) {
	    Warning( here, "Event %d: Multiple hits of primary particle "
		     "in plane %d\nShould never happen. Call expert.",
		     event_num, c.fPlane );
	  } else {
	    if( !theTrack->TestHitBit(0) ) {
	      // If no plane 0 position information is yet recorded,
	      // see if the position information can be improved.
	      if( c.fPlane == 0 ) {
		// Jackpot; we got plane 0 now
		if( theTrack->Update(c) )
		  return HED_FATAL;
	      }
	      else {
		// Got a plane closer to 0 than what's been seen so far?
		Int_t ibit = 1;
		for( ; ibit < NPLANES; ++ibit ) {
		  if( theTrack->TestHitBit(ibit) )
		    break;
		}
		if( ibit == NPLANES ) {
		  Error( here, "No hit bit set in existing back track. "
			 "Should never happen. Call expert." );
		  theTrack->Print();
		  return HED_FATAL;
		}
		if( c.fPlane < ibit )
		  if( theTrack->Update(c) )
		    return HED_FATAL;
	      }
	    }
	    theTrack->SetHitBit( c.fPlane );
	  }
	  break;
	}
      }
      if( ib == nback )
	new( (*fBackTracks)[nback] ) TSolSimBackTrack(c);
    }
  }
  
  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNTracks() << endl;
  //fTracks.Print();

  return HED_OK;
}

//-----------------------------------------------------------------------------
TSolSimGEMHit::TSolSimGEMHit( const TSolSimEvent::GEMCluster& c )
  : fID(c.fID), fSector(c.fSector), fPlane(c.fPlane), fType(c.fType),
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
  : fType(c.fType), fPID(c.fPID), fSector(c.fSector), fHitBits(0)
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

  // z position of first tracker plane.
  // FIXME: Get this from plane object
  static double z0 = 1.536914;

  if( fType != c.fType || fPID != c.fPID || fSector != c.fSector ) {
    Error( here, "Updating with inconsistent GEMCluster data: "
	   "type = %d/%d, pid = %d/%d, sector = %d/%d.\n"
	   "Should never happen. Call expert.",
	   fType, c.fType, fPID, c.fPID, fSector, c.fSector );
    return -1;
  }
  if( c.fPlane > 0 ) {
    Double_t dz = c.fMCpos.Z() - z0;
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
