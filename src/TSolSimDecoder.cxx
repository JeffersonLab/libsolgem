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
#include "TSolSimEvent.h"
#include "THaCrateMap.h"
#include "THaBenchmark.h"
#include "VarDef.h"
#include "TClonesArray.h"

using namespace std;

// Prefix of our own gloabl variables (true MC data)
static const char* const MC_PREFIX = "MC.";

//-----------------------------------------------------------------------------
TSolSimDecoder::TSolSimDecoder() : fIsSetup(false)
{
  // Constructor

  DefineVariables();
}

//-----------------------------------------------------------------------------
TSolSimDecoder::~TSolSimDecoder() {

  DefineVariables( THaAnalysisObject::kDelete );
}

//-----------------------------------------------------------------------------
Int_t TSolSimDecoder::DefineVariables( THaAnalysisObject::EMode mode )
{
  const char* const here = "TSolSimDecoder::DefineVariables";

  if( mode == THaAnalysisObject::kDefine && fIsSetup ) 
    return THaAnalysisObject::kOK;
  fIsSetup = ( mode == THaAnalysisObject::kDefine );

  RVarDef vars[] = {
    { "tr.n",    "Number of tracks",     "GetNTracks()" },
    { "tr.x",    "Track origin x (m)",   "fTracks.TSolSimTrack.TX()" },
    { "tr.y",    "Track origin y (m)",   "fTracks.TSolSimTrack.TY()" },
    { "tr.th",   "Track tan(theta)",     "fTracks.TSolSimTrack.TTheta()" },
    { "tr.ph",   "Track tan(phi)",       "fTracks.TSolSimTrack.TPhi()" },
    { "tr.p",    "Track momentum (GeV)", "fTracks.TSolSimTrack.P() "},
    { "tr.hit1", "Chamber # of 1st hit", "fTracks.TSolSimTrack.fHit1" },
    { "tr.pid",  "Track PID (PDG)",      "fTracks.TSolSimTrack.fPID" },
    { "tr.t0",   "Track time offset",    "fTracks.TSolSimTrack.fTimeOffset" },
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

  // Never delete the tracks, only clear the list. The tracks are deleted
  // in TSolSimFile::ReadEvent() by the call to fEvent->Clear().
  fTracks.Clear("nodelete");
}

//-----------------------------------------------------------------------------
int TSolSimDecoder::LoadEvent(const int* evbuffer, THaCrateMap* map)
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'

  fMap = map;

  // Local copy of evbuffer pointer - any good use for it?
  buffer = evbuffer;

  // Cast the evbuffer pointer back to exactly the event type that is present
  // in the input file (in TSolSimFile). The pointer-to-integer is there
  // just for compatibility with the standard decoder.
  // Note: simEvent can't be constant - ROOT does not like to iterate
  // over const TList.
  TSolSimEvent* simEvent = (TSolSimEvent*)(evbuffer);

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

  // Decode the digitized data.  Populate crateslot array.
  for( Int_t i = 0; i < simEvent->fNStrips; i++) {

    Int_t roc  = simEvent->fStpGEM[i];
    Int_t slot = simEvent->fStpNum[i] / CHAN_PER_SLOT;
    Int_t chan = simEvent->fStpNum[i] % CHAN_PER_SLOT;

    for( Int_t k = 0; k < MC_MAXSAMP; k++ ) {
      Int_t raw = simEvent->fStpADC[k][i];
      if (crateslot[idx(roc,slot)]->loadData("adc",chan,raw,raw)
	  == SD_ERR) return HED_ERR;
    }
  }

  // Extract MC track info, so we can access it via global variables
  // The list of tracks is already part of the event - no need to generate
  // our own tracks here. 
  // FIXME: However, we have to copy the list because the global variable system 
  // cannot handle variable pointers.

  TClonesArray* tracks = simEvent->fMCTracks;
  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    fTracks.Add( tracks->UncheckedAt(i) );
  }

  if( fDoBench ) fBench->Stop("physics_decode");

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNTracks() << endl;
  //fTracks.Print();

  return HED_OK;
}

//-----------------------------------------------------------------------------
ClassImp(TSolSimDecoder)
