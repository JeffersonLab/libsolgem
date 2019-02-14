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
#include "TSBSDBManager.h"
#include "ha_compiledata.h"

#include "TError.h"
#include "TSystem.h"
#include "TMath.h"
#include "TDatabasePDG.h"
//#include "TRandom3.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace Podd;

//EFuchey: 2016/12/10: it is necessary to declare the TSBSDBManager as a static instance here 
// (and not make it a member) because it is used by functions whic are defined as "static inline".
static TSBSDBManager* fManager = TSBSDBManager::GetInstance();
//static TRandom3* Rdec = new TRandom3(0);
static const Int_t kPrimaryType = 1, kPrimarySource = 0;
// Projection types must match the definitions in TreeSearch
enum EProjType { kUPlane = 0, kVPlane =1, kXPlane = 2, kYPlane = 3};

typedef vector<int>::size_type vsiz_t;

// TODO: Have dbconvert write out MAXSLOT (and possibly other parameters)
//  to the database to allow client to understand the generated detector maps.
// FIXME: The number 30 is hardcoded in dbconvert
static const Int_t SIM_MAXSLOT = TMath::Min(Decoder::MAXSLOT,30);

// Hard coded stuff for the time being...
// what we will need is to switch all that stuff to SBS-offline anyway
static const double Mp = 0.938272; //GeV
static const double Lx_scint_CDET = 0.51; //m
static const double Ly_scint_CDET = 0.005; //5 mm
static const double mu_p = 2.793; // nuclear magneton
static const double PI = TMath::Pi();
static const double SBS_tracker_pitch=5.0*PI/180.0; //5 degrees

static const double ECAL_phe_per_GeV=300.0;
static const double ECAL_max_cell_size = 0.042; //meters
static const double X0_ECAL = 0.0274; //radiation length
static const double Ec_ECAL = 0.015; //GeV
static const double yoff_ECAL = 0.0077; //meters, average y deflection in SBS fringe field, to be SUBTRACTED from recconstructed shower coordinate!
static const double sigx_ECAL = 0.008; //meters, shower x coordinate resolution
static const double sigy_ECAL = 0.006; //meters, shower y coordinate resolution
static const double Ltgt = 0.4; //meters
static const double sigy_CDET = 0.003;

// Bin width for vertex z "filtering": 
static const double sig_vz = 0.0064;//in m
static const double vz_bin_width = 3.0*0.0064;
static const double vz_scan_stepsize = 0.3;

static const double dE_E_MPV = 0.13; //Most probable electron energy loss before reaching ECAL:


//-----------------------------------------------------------------------------
TSBSSimDecoder::TSBSSimDecoder()
{
  // Constructor

  fMCHits     = new TClonesArray( "TSBSSimGEMHit",    200 );
  fMCTracks   = new TClonesArray( "TSBSSimTrack",       1 );
  fBackTracks = new TClonesArray( "TSBSSimBackTrack",   5 );
  
  DefineVariables();

  gSystem->Load("libEG.so");  // for TDatabasePDG
  
  for (int i=0; i<fManager->GetNSigParticle(); i++){
    fSignalInfo.push_back(SignalInfo(fManager->GetSigPID(i),
                                     fManager->GetSigTID(i)));
  }
  if(fManager->Getg4sbsDetectorType()==3 && fManager->DoCalo()){
    cout << "Det type = " << fManager->Getg4sbsDetectorType() << " and manager = " << fManager->DoCalo() << endl;
    if(load_shower_profiles("ECAL_shower_profiles.txt"))
      cout << "ECal shower profiles successfully loaded" << endl;
  }
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
    //{ "tr.n",      "Number of tracks",      "GetNMCTracks()" },  // already defined in Podd::SimDecoder
    { "tr.vx",     "Track origin x (m)",    "fMCTracks.TSBSSimTrack.targetX()" },
    { "tr.vy",     "Track origin y (m)",    "fMCTracks.TSBSSimTrack.targetY()" },
    { "tr.vz",     "Track origin z (m)",    "fMCTracks.TSBSSimTrack.targetZ()" },
    { "tr.p",      "Track momentum (GeV)",  "fMCTracks.TSBSSimTrack.targetP() "},
    { "tr.theta",  "Track theta_p (rad)",   "fMCTracks.TSBSSimTrack.targetPTheta()" },
    { "tr.phi",    "Track phi_p (rad)",     "fMCTracks.TSBSSimTrack.targetPPhi()" },
    { "tr.ptarx", "Track px(GeV)",         "fMCTracks.TSBSSimTrack.targetPX()" },
    { "tr.ptary", "Track py(GeV)",         "fMCTracks.TSBSSimTrack.targetPY()" },
    { "tr.ptarz", "Track pz(GeV)",         "fMCTracks.TSBSSimTrack.targetPZ()" },
    //variables in GEM detector coordinates
    { "tr.x",     "Track x in Transport",    "fMCTracks.TSBSSimTrack.VX()" },
    { "tr.y",     "Track y in Transport",    "fMCTracks.TSBSSimTrack.VY()" },
    { "tr.p_transport",   "Track momentum in Transport", "fMCTracks.TSBSSimTrack.P()"},
    { "tr.px",   "Track x_momentum in Transport", "fMCTracks.TSBSSimTrack.PX()"},
    { "tr.py",   "Track y_momentum in Transport", "fMCTracks.TSBSSimTrack.PY()"},
    { "tr.pz",   "Track z_momentum in Transport", "fMCTracks.TSBSSimTrack.PZ()"},

    { "tr.pid",    "Track PID (PDG)",       "fMCTracks.TSBSSimTrack.fPID" },
    { "tr.num",    "GEANT track number",    "fMCTracks.TSBSSimTrack.fNumber" },
    { "tr.planes", "Bitpattern of planes hit", "fMCTracks.TSBSSimTrack.fHitBits" },
    { "tr.nhits",  "Number of tracker hits","fMCTracks.TSBSSimTrack.fNHits" },
    { "tr.nfound", "Number of hits found",  "fMCTracks.TSBSSimTrack.fNHitsFound" },
    { "tr.flags",  "Reconstruction status", "fMCTracks.TSBSSimTrack.fReconFlags" },

    // Results of fit to MC points - measures multiple scattering
    // Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
    // refer to comment in TSBSSimEvent.h l. 30-32
    // { "tr.mcfit.r",     "Track x from MC fit [m]", "fMCTracks.TSBSSimTrack.MCFitR()" },
    // { "tr.mcfit.phi",   "Track phi from MC fit [rad]", "fMCTracks.TSBSSimTrack.MCFitPhi()" },
    // { "tr.mcfit.thdir", "Track dir theta from MC fit [rad]", "fMCTracks.TSBSSimTrack.MCFitThetaDir()" },
    // { "tr.mcfit.phdir", "Track x from MC fit [rad]", "fMCTracks.TSBSSimTrack.MCFitPhiDir()" },
    // { "tr.mcfit.x",     "Track x from MC fit [m]",       "fMCTracks.TSBSSimTrack.MCFitX_print()" },
    { "tr.mcfit.x",     "Track x from MC fit [m]",       "fMCTracks.TSBSSimTrack.fMCFitPar[0]" },
    { "tr.mcfit.xdir",  "Track dir x from MC fit [rad]", "fMCTracks.TSBSSimTrack.fMCFitPar[1]" },
    { "tr.mcfit.y",     "Track y from MC fit [rad]",     "fMCTracks.TSBSSimTrack.fMCFitPar[2]" },
    { "tr.mcfit.ydir",  "Track dir y from MC fit [rad]", "fMCTracks.TSBSSimTrack.fMCFitPar[3]" },
    { "tr.mcfit.chi2",  "Chi2 of MC fit",                "fMCTracks.TSBSSimTrack.fMCFitPar[4]" },
    { "tr.mcfit.ndof",  "NDoF of MC fit",                "fMCTracks.TSBSSimTrack.fMCFitPar[5]" },
    { "tr.mcfit.vx",    "Vertex x from MC fit [m]",      "fMCTracks.TSBSSimTrack.fMCFitPar[6]" },
    { "tr.mcfit.vy",    "Vertex y from MC fit [m]",      "fMCTracks.TSBSSimTrack.fMCFitPar[7]" },
    { "tr.mcfit.vz",    "Vertex z from MC fit [m]",      "fMCTracks.TSBSSimTrack.fMCFitPar[8]" },

    // Results of fit to reconstructed MC hits - checks hit resolution effects
    // independent of track finding
    // Those below are not useful for SBS, which needs X, Y, Xdir, Ydir (unless otherwise demonstrated)
    // refer to comment in TSBSSimEvent.h l. 30-32
    // { "tr.fit.r",     "Track x from rec hit fit [m]", "fMCTracks.TSBSSimTrack.RcFitR()" },
    // { "tr.fit.phi",   "Track phi from rec hit fit [rad]", "fMCTracks.TSBSSimTrack.RcFitPhi()" },
    // { "tr.fit.thdir", "Track dir theta from rec hit fit [rad]", "fMCTracks.TSBSSimTrack.RcFitThetaDir()" },
    // { "tr.fit.phdir", "Track x from rec hit fit [rad]", "fMCTracks.TSBSSimTrack.RcFitPhiDir()" },
    { "tr.fit.x",     "Track x from rec hit fit [m]",       "fMCTracks.TSBSSimTrack.fRcFitPar[0]" },
    { "tr.fit.xdir",  "Track dir x from rec hit fit [rad]", "fMCTracks.TSBSSimTrack.fRcFitPar[1]" },
    { "tr.fit.y",     "Track y from rec hit fit [rad]",     "fMCTracks.TSBSSimTrack.fRcFitPar[2]" },
    { "tr.fit.ydir",  "Track dir y from rec hit fit [rad]", "fMCTracks.TSBSSimTrack.fRcFitPar[3]" },
    { "tr.fit.chi2",  "Chi2 of rec hit fit",                "fMCTracks.TSBSSimTrack.fRcFitPar[4]" },
    { "tr.fit.ndof",  "NDoF of rec hit fit",                "fMCTracks.TSBSSimTrack.fRcFitPar[5]" },
    { "tr.fit.vx",    "Vertex x from rec hit fit [m]",      "fMCTracks.TSBSSimTrack.fRcFitPar[6]" },
    { "tr.fit.vy",    "Vertex y from rec hit fit [m]",      "fMCTracks.TSBSSimTrack.fRcFitPar[7]" },
    { "tr.fit.vz",    "Vertex z from rec hit fit [m]",      "fMCTracks.TSBSSimTrack.fRcFitPar[8]" },

    // "Back tracks": hits of the primary particle in the first tracker plane
    { "btr.n",     "Number of back tracks",     "GetNBackTracks()" },
    { "btr.pid",   "Track PID (PDG)",           "fBackTracks.TSBSSimBackTrack.fPID" },
    { "btr.num",   "GEANT particle number",     "fBackTracks.TSBSSimBackTrack.fType" },
    { "btr.planes","Bitpattern of planes hit",  "fBackTracks.TSBSSimBackTrack.fHitBits" },
    { "btr.ufail", "Undigitized u planes",      "fBackTracks.TSBSSimBackTrack.fUfailBits" },
    { "btr.vfail", "Undigitized v planes",      "fBackTracks.TSBSSimBackTrack.fVfailBits" },
    { "btr.sect",  "Sector number",             "fBackTracks.TSBSSimBackTrack.fSector" },
    { "btr.p",     "Track momentum (GeV)",      "fBackTracks.TSBSSimBackTrack.P() "},
    // Track position in Cartesian/TRANSPORT coordinates, optimal for SBS, not for SoLID
    { "btr.x",     "Track pos lab x [m]",       "fBackTracks.TSBSSimBackTrack.X()" },
    { "btr.y",     "Track pos lab y [m]",       "fBackTracks.TSBSSimBackTrack.Y()" },
    { "btr.th",    "Track dir tan(theta)",      "fBackTracks.TSBSSimBackTrack.ThetaT()" },
    { "btr.ph",    "Track dir tan(phi)",        "fBackTracks.TSBSSimBackTrack.PhiT()" },
    // Track position and direction in cylindrical coordinates, good for SoLID
    // { "btr.r",     "Track pos lab r_trans (m)", "fBackTracks.TSBSSimBackTrack.R()" },
    // { "btr.theta", "Track pos lab theta [rad]", "fBackTracks.TSBSSimBackTrack.Theta()" },
    // { "btr.phi",   "Track pos lab phi [rad]",   "fBackTracks.TSBSSimBackTrack.Phi()" },
    // { "btr.thdir", "Track dir theta [rad]",     "fBackTracks.TSBSSimBackTrack.ThetaDir()" },
    // { "btr.phdir", "Track dir phi [rad]",       "fBackTracks.TSBSSimBackTrack.PhiDir()" },
    // Hit coordinates in first tracker plane, relative to plane origin
    { "btr.hx",    "Track pos plane x [m]",     "fBackTracks.TSBSSimBackTrack.HX()" },
    { "btr.hy",    "Track pos plane y [m]",     "fBackTracks.TSBSSimBackTrack.HY()" },

    // Digitized hits registered in the GEMs
    //    { "hit.n",     "Number of MC hits",          "GetNMCHits()" },
    { "hit.id",    "MC hit number",              "fMCHits.TSBSSimGEMHit.fID" },
    { "hit.sect",  "MC hit sector",              "fMCHits.TSBSSimGEMHit.fSector" },
    { "hit.rsect", "MC hit non-mapped sector",   "fMCHits.TSBSSimGEMHit.fRealSector" },
    { "hit.plane", "MC hit plane",               "fMCHits.TSBSSimGEMHit.fPlane" },
    { "hit.module","MC hit module",              "fMCHits.TSBSSimGEMHit.fModule" },
    { "hit.src",   "MC data set source",         "fMCHits.TSBSSimGEMHit.fSource" },
    { "hit.type",  "MC hit GEANT counter",       "fMCHits.TSBSSimGEMHit.fType" },
    { "hit.pid",   "MC hit PID (PDG)",           "fMCHits.TSBSSimGEMHit.fPID" },
    { "hit.p",     "MC hit particle mom [GeV]",  "fMCHits.TSBSSimGEMHit.P()" },
    { "hit.x",     "MC hit lab x position [m]",  "fMCHits.TSBSSimGEMHit.X()" },
    { "hit.y",     "MC hit lab y position [m]",  "fMCHits.TSBSSimGEMHit.Y()" },
    { "hit.z",     "MC hit lab z position [m]",  "fMCHits.TSBSSimGEMHit.Z()" },
    // Hit position in cylindrical/spherical coordinates, good for SoLID
    // { "hit.r",     "MC hit lab r [m]",           "fMCHits.TSBSSimGEMHit.R()" },
    // { "hit.theta", "MC hit lab theta [rad]",     "fMCHits.TSBSSimGEMHit.Theta()" },
    // { "hit.phi",   "MC hit lab phi [rad]",       "fMCHits.TSBSSimGEMHit.Phi()" },
    { "hit.charge","MC hit cluster charge",      "fMCHits.TSBSSimGEMHit.fCharge" },
    { "hit.time",  "MC hit time offset [s]",     "fMCHits.TSBSSimGEMHit.fTime" },
    { "hit.usz",   "MC hit u cluster size",      "fMCHits.TSBSSimGEMHit.fUSize" },
    { "hit.ustart","MC hit u cluster 1st strip", "fMCHits.TSBSSimGEMHit.fUStart" },
    { "hit.upos",  "MC hit u cluster center [m]","fMCHits.TSBSSimGEMHit.fUPos" },
    { "hit.vsz",   "MC hit v cluster size",      "fMCHits.TSBSSimGEMHit.fVSize" },
    { "hit.vstart","MC hit v cluster 1st strip", "fMCHits.TSBSSimGEMHit.fVStart" },
    { "hit.vpos",  "MC hit v cluster center [m]","fMCHits.TSBSSimGEMHit.fVPos" },
    
    { "pt.fmctrk", "MC point track number",      "fMCPoints.Podd::MCTrackPoint.fMCTrack" },
    
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
  
  // cout << "Chan per slot ? " << fManager->GetChanPerSlot() << endl;
  // cout << "Module per readout ? " << fManager->GetModulesPerReadOut() << endl;
  // cout << "N readout ? " << fManager->GetNReadOut() << ", N Chambers ? " << fManager->GetNChamber() << endl;
  // cout << "Chambers per crate ? " << fManager->GetChambersPerCrate() << endl;
  // cout << "Module per readout ? " << fManager->GetModulesPerChamber() << endl;
  
  div_t d = div( s_chan, fManager->GetChanPerSlot() );
  Int_t module = d.quot;
  chan = d.rem;
  Int_t ix = module +
    fManager->GetModulesPerReadOut()*( s_proj + fManager->GetNReadOut()*( s_plane + fManager->GetNChamber()*s_sector ));
  
  //cout << "StripToROC: module " << module << ", ix " << ix << endl;
  
  d = div( ix, fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
void StripToROCD( Int_t s_plane, Int_t s_module, Int_t s_proj,
		 Int_t s_chan,
		 Int_t& crate, Int_t& slot, Int_t& chan )
{
  div_t d = div( s_chan, fManager->GetChanPerSlot() );
  //  Int_t module = d.quot;
  chan = d.rem;
  //total slot id
  Int_t ix = s_proj + 2*( s_module + fManager->GetNModule(s_plane-1)*s_plane );
  
  //  cout << "StripToROC: module " << module << ", ix " << ix << Decoder::MAXSLOT<<endl;
  
  d = div( ix, SIM_MAXSLOT);//fManager->GetChambersPerCrate()*fManager->GetModulesPerChamber() );
  crate = d.quot;
  slot  = d.rem;
}

//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan +
    fManager->GetChanPerSlot()*( slot + SIM_MAXSLOT*crate );
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
std::vector<std::vector<Double_t>> TSBSSimDecoder::GetAllMCHits() const
{
  std::vector<std::vector<Double_t>> hits;
  std::vector<Double_t> vtemp = {0,0,0,0,0,0};//v[0]--posx, v[1]--posy, v[2]--charge, v[3]--planeID v[4]--moduleID v[5]time_zero
  assert( buffer );       // Must still have the event buffer
  const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);
  for(size_t i=0;i<simEvent->fGEMClust.size();i++){
    const TSBSSimEvent::GEMCluster& clust = simEvent->fGEMClust[i];
    if(clust.fSource!=0){continue;}
    vtemp[0] = clust.fMCpos.X();
    vtemp[1] = clust.fMCpos.Y();
    vtemp[2] = clust.fCharge;
    vtemp[3] = clust.fPlane;
    vtemp[4] = clust.fModule;
    vtemp[5] = clust.fTime;
    // cout<<"######## "<<vtemp[1]<<endl;
    hits.push_back(vtemp);
  }

  return hits;
}

//-----------------------------------------------------------------------------
TSBSMCHitInfo TSBSSimDecoder::GetSBSMCHitInfo( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Get MC truth info for the given hardware channel
//  const char* const here = "TSBSSimDecoder::GetSBSMCHitInfo";

  Int_t istrip = StripFromROC( crate, slot, chan );
  assert( istrip >= 0 );  // else logic error in caller or bad fStripMap
  
  
  assert( buffer );       // Must still have the event buffer
  const TSBSSimEvent* simEvent = reinterpret_cast<const TSBSSimEvent*>(buffer);
  
  assert( static_cast<vsiz_t>(istrip) < simEvent->fGEMStrips.size() );
  const TSBSSimEvent::DigiGEMStrip& strip = simEvent->fGEMStrips[istrip];
  assert( strip.fProj >= 0 && strip.fProj < fManager->GetNReadOut() );
  
  TSBSMCHitInfo mc;
  mc.fSigType = strip.fSigType;
  //if(strip.fSigType==1)cout << "TSBSSimDecoder::GetSBSMCHitInfo mc.fSigType " << mc.fSigType << " strip ? " << istrip << " crate " << crate << " slot " << slot << " chan " << chan << endl;
  // cout<<kPrimaryStrip<<" "<<kSecondaryStrip<<" "<<kInducedStrip<<endl;getchar();0 1 2
    // if(strip.fProj==0 && strip.fPlane==4 && strip.fTime1>50.0)
  //   printf("%f \n", strip.fTime1);
  
  //for cross talk
  if (TESTBIT(strip.fSigType, kInducedStrip) && !TESTBIT(strip.fSigType, kPrimaryStrip) &&
      !TESTBIT(strip.fSigType, kSecondaryStrip) ){
    mc.fMCTrack = 0;
    mc.fMCPos = fManager->GetPosFromModuleStrip(strip.fProj, strip.fPlane, strip.fModule, strip.fChan);
    mc.fMCTime = strip.fTime1;
    
    //cout << "strip = " << strip.fChan << ", time = " << mc.fMCTime << ", pos = " <<  mc.fMCPos << endl;
    return mc;
  }

  mc.fMCCharge = strip.fCharge;
  // cout<<"strip: "<<istrip<<" cc: "<<mc.fMCCharge<<endl;
  Double_t nOverlapSignal = 0.;
  // cout<<strip.fClusters.GetSize()<<" # "<<strip.fClusterRatio[0].GetSize()<<" # "<<strip.fClusterRatio[1].GetSize()<<" # "<<strip.fClusterRatio[2].GetSize()<<" # "<<strip.fClusterRatio[3].GetSize()<<endl;getchar();
  for( Int_t i = 0; i<strip.fClusters.GetSize(); ++i ) {
    Int_t iclust = strip.fClusters[i] - 1;  // yeah, array index = clusterID - 1

   
    //   cout<<mc.vClusterID.size()<<" : "<<mc.vClusterADC[1].size()<<endl;getchar(); 

    assert( iclust >= 0 && static_cast<vsiz_t>(iclust) < simEvent->fGEMClust.size() );
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[iclust];
    mc.vClusterID.push_back(iclust);

    //work here! remove after finish!
    // add cluster type info here to tell whether a cluster is a primary or background using "fSource" in simevent.h !!
    // then go to get b_overtotal in gemplane.cxx and so u can tell whether a cluster has background in it. and how much! great!
    mc.vClusterType.push_back(c.fSource);




    //
    mc.vClusterPeakTime.push_back(c.fTime);
    mc.vClusterPos.push_back(c.fHitpos);
    mc.vClusterCharge.push_back(c.fCharge);

    mc.vClusterStripWeight.push_back(strip.fStripWeightInCluster[i]);
    // cout<<strip.fStripWeightInCluster[i]<<endl;getchar();
    for(Int_t its=0;its<6;its++)
      {
	mc.vClusterADC[its].push_back(strip.fClusterRatio[its][i]);
      }

    assert( c.fID == iclust+1 );
    assert( strip.fPlane == c.fPlane && strip.fSector == c.fSector );
    Int_t signalID = -1;
    for (unsigned int ii = 0; ii<fSignalInfo.size(); ii++){
      //cout << "c.fType " << c.fType << " fSignalInfo.at(ii).tid " << fSignalInfo.at(ii).tid << ", c.fPID " << c.fPID << " fSignalInfo.at(ii).pid " << fSignalInfo.at(ii).pid << endl;
      //if (c.fType == fSignalInfo.at(ii).tid && c.fPID == fSignalInfo.at(ii).pid) // cluster_type(primary or secondary) == type_requested(primary) && particle == partical_requested 
      // cluster type has _nothing_ to do with signal track ID !!!!!!!! 
      if (c.fType == 1 && c.fPID == fSignalInfo.at(ii).pid) // cluster_type(primary or secondary) == type_requested(primary) && particle == partical_requested
	signalID = ii;
    }
    // cout << "Plane " << strip.fPlane << ", proj (x: 0, y: 1) " << strip.fProj 
    //  	 << ": pos[proj] = "  << c.fXProj[strip.fProj] << endl;
    //cout << "signalID " << signalID << "; " << c.fSource << " ==? " << kPrimarySource << endl;
    if( signalID >= 0 && c.fSource == kPrimarySource ) {
      if( mc.fMCTrack > 0 ) {
        //this means that there two signal hits overlapping
        //for now I keep the fMCTrack to the first one, by average the fMCPos nad fMCTime
        //Weizhi Xiong
        //assert(manager->GetNSigParticle() > 1); //otherwise should not happen
        
        mc.fMCPos += c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
        mc.fMCTime += c.fTime; 
      }else{
        // Strip contains a contribution from a primary particle hit :)
        mc.fMCTrack = fSignalInfo.at(signalID).tid; 
        mc.fMCPos   = c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
        mc.fMCTime  = c.fTime;
      }
      nOverlapSignal++;
    } else {
      ++mc.fContam;
      if( mc.fMCTrack == 0 ) {
	mc.fMCPos += c.fXProj[strip.fProj]+(1-strip.fProj)*fManager->GetXOffset(c.fPlane,c.fModule);
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
  }else{
    mc.fMCPos /= nOverlapSignal;
    mc.fMCTime /= nOverlapSignal;
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
  // cout<<"Tsbs sim decoder DoLoadEvent"<<endl;
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
 
  /* 
  cout<<simEvent->fGEMClust.size()<<endl;getchar();
  for(Int_t i=0;i<simEvent->fGEMClust.size();i++){
    cout<<"Plane: "<<simEvent->fGEMClust[i].fPlane<<" type: "<<simEvent->fGEMClust[i].fType<<"   x: "<<simEvent->fGEMClust[i].fMCpos.X()
	<<"   y: "<<simEvent->fGEMClust[i].fMCpos.Y()<<endl;
  }
  */

  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    if( (ret = init_cmap()) != HED_OK )
      return ret;
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
    if( (ret = init_slotdata()) != HED_OK)
#else
    if( (ret = init_slotdata(fMap)) != HED_OK)
#endif
      return ret;
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

  // Event weight
  fWeight = simEvent->fWeight;

  //
  if( fDoBench ) fBench->Begin("physics_decode");

  // Decode the digitized strip data.  Populate crateslot array.
  for( vector<TSBSSimEvent::DigiGEMStrip>::size_type i = 0;
       i < simEvent->fGEMStrips.size(); i++) {
    //cout << "i " << i << endl;
    const TSBSSimEvent::DigiGEMStrip& s = simEvent->fGEMStrips[i];
    Int_t crate, slot, chan;
    //cout << "striptoroc: " << endl;
    //StripToROC( s.fPlane, s.fSector, s.fProj, s.fChan, crate, slot, chan );
    StripToROCD( s.fPlane, s.fModule, s.fProj, s.fChan, crate, slot, chan );
    //cout<<"Plane: "<<s.fPlane<<" Module:  "<<s.fModule<<" projection:  "<<s.fProj<<" channel:  "<<s.fChan<<" sigType: "<<s.fSigType<<endl;
    //cout << "crate = " << crate << ", slot = " << slot << ", chan " << chan << endl;getchar();
    //cout << "samples: " << endl;
    for( Int_t k = 0; k < s.fNsamp; k++ ) { 
      Int_t raw = s.fADC[k];
      //cout << raw << " ### ";
      
      if( crateslot[idx(crate,slot)]->loadData("adc",chan,raw,raw) == SD_ERR )
	return HED_ERR;
    }
    //    cout << endl<<endl;
    //cout << "stripmap : " << endl;
    // Build map from ROC address to strip index. This is needed to extract
    // the MC truth info later in the tracking detector decoder via GetMCChanInfo.
#ifndef NDEBUG
    pair<StripMap_t::const_iterator,bool> ins =
#endif
      fStripMap.insert( make_pair( MakeROCKey(crate,slot,chan), i ) );
    //cout<<crate<<" "<<slot<<" "<<chan<<endl;
    // cout<<MakeROCKey(crate,slot,chan)<<endl;
    // getchar();
    
    // cout << "ROC key inserted in strip map " << endl;
    //  cout << "ins.second ? " << ins.second << endl;
    assert( ins.second );
  }
  
  // Create lists of two types of tracks:
  // 1) Physics tracks, as generated at the target
  // 2) "Back tracks": hits in any GEM plane from the primary particle

  // Physics tracks. We need to copy them here so we can export them as global
  // variables.
  TClonesArray* tracks = simEvent->fMCTracks;
  assert( tracks );

  double trkProjCaloX, trkProjCaloY;

  Int_t itrack = 0;

  double E_evt = 0;//-2*fManager->GetCaloThreshold();
  //
  //if(fManager->DoCalo()){
  //E_evt = 0;
  if(simEvent->fECalClusters.size()>0){
    const TSBSECalCluster& eCalHit = simEvent -> fECalClusters[0];
    E_evt = eCalHit.GetEnergy();
  }
  
  if(fManager->DoCalo() && fabs(E_evt)<fManager->GetCaloThreshold())return HED_ERR;

  for( Int_t i = 0; i < tracks->GetLast()+1; i++ ) {
    TSBSSimTrack* trk = static_cast<TSBSSimTrack*>(tracks->UncheckedAt(i));
    
    trkProjCaloX = trk->fMomentum.X()/trk->P()*(fManager->GetCaloZ()-0.8) + trk->fOrigin.X()/1000;
    trkProjCaloY = trk->fMomentum.Y()/trk->P()*(fManager->GetCaloZ()-0.8) + trk->fOrigin.Y()/1000;
    //cout<<"MC trk projected calo position: "<<trkProjCaloX<<" : "<<trkProjCaloY<<endl;
    
    trkProjCaloX = trk->fMomentum.X()/trk->P()*(1.3-0.8) + trk->fOrigin.X()/1000;
    trkProjCaloY = trk->fMomentum.Y()/trk->P()*(1.3-0.8) + trk->fOrigin.Y()/1000;
    //cout<<"MC trk projected 4th  GEM plane position: "<<trkProjCaloX<<" : "<<trkProjCaloY<<endl;
    
    trkProjCaloX = trk->fMomentum.X()/trk->P()*(2.33-0.8) + trk->fOrigin.X()/1000;
    trkProjCaloY = trk->fMomentum.Y()/trk->P()*(2.33-0.8) + trk->fOrigin.Y()/1000;
    //cout<<"MC trk projected last GEM plane position: "<<trkProjCaloX<<" : "<<trkProjCaloY<<endl;
    
    new( (*fMCTracks)[i] ) TSBSSimTrack(*trk);
    if(i==0)itrack = trk->fNumber;//even better
  }
  assert( GetNMCTracks() > 0 );
  



  // MC hit data ("clusters") and "back tracks"
  Int_t best_primary = -1, best_primary_plane = fManager->GetNChamber(), primary_sector = -1;
  UInt_t primary_hitbits = 0, ufail = 0, vfail = 0;
  for( vector<TSBSSimEvent::GEMCluster>::size_type i = 0;
       i < simEvent->fGEMClust.size(); ++i ) {
    const TSBSSimEvent::GEMCluster& c = simEvent->fGEMClust[i];

    if( c.fPlane < 0 || c.fPlane >= fManager->GetNChamber() ) {
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
      //Int_t itrack = 1;//dont' use the by default value...
      //Int_t itrack = fManager->GetSigTID(0);//better
      //cout << "itrack ??? " << itrack << endl;
      primary_sector = c.fSector;
      MCTrackPoint* upt = // kUPlane changed to kXPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kXPlane,
							  c.fMCpos, c.fP );
      upt->fMCTime = c.fTime;
      MCTrackPoint* vpt =// kVPlane changed to kYPlane: necessary to match TreeSearch EProjType
	new( (*fMCPoints)[GetNMCPoints()] ) MCTrackPoint( itrack,
							  c.fPlane, kYPlane,
							  c.fMCpos, c.fP );
      vpt->fMCTime = c.fTime;
      
      // //debug...
      // cout << "TSBSSimDecoder.cxx: Print MC points " << endl;
      // cout << "kXplane ? " << kXPlane << endl;
      // upt->Print("");
      // cout << "kVYlane ? " << kYPlane << endl;
      // vpt->Print("");
      
      // Keep bitpattern of planes crossed by this primary
      SETBIT(primary_hitbits,c.fPlane);

      //cout << "Plane number " << c.fPlane << ", primary hitbits " << primary_hitbits << endl; 
      
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

    //cout << "Backtrack primary hitbits " << primary_hitbits << endl;
    
    btr->SetHitBits(primary_hitbits);
    btr->SetUfailBits(ufail);
    btr->SetVfailBits(vfail);

    // Use the back track to emulate calorimeter hits.
    // Assumptions:
    // - Only tracks crossing all fManager->GetNChamber() GEMs (points in all planes)
    //   make a calorimeter hit. This is a crude model for the trigger.
    // - The track propagates without deflection from the last GEM plane
    //   to the front of the emulated calorimeter.
    // - The measured calorimeter position is independent of the incident
    //   track angle.
    
    if( fManager->DoCalo() ){// && trk->fNHits == 2*fManager->GetNChamber() ) {
	double tempCaloX;
	double tempCaloY;
	Int_t flag;
	double E;
	double x_S, y_S;
	
	//all the stuff for ECal analysis
	double xmom, ymom;
	double xECal, yECal;
	
	//double tempScintX;
	//double tempScintY;
	//Int_t tempPlane;
	double kpx_xc, kpy_xc, kpz_xc;
	double  kp, thetak, phik;
	double  kp_xc, thetak_xc, phik_xc;
	double  pp_xc, thetap_xc, phip_xc;
	//double ppx, ppy, ppz;
	double pp, thetap, phip;
		
	double xtar, ytar, xptar, yptar;
	double fterm;
	double xfp, yfp, xpfp, ypfp;
	double x_S_2, y_S_2;

	vector<double> xfinal,yfinal,zfinal,wxfinal,wyfinal;
	double eX, eY, eXp, eYp;
	
	double vz;
	
	//this is becoming very dirty... :/
	const double k0 = 11.0;
	const double th_earm = 29.0*TMath::DegToRad();
	const double z_earm[3] = {4.5, 4.08, 4.13};
	//{5.02, 4.08, 4.13};
	//double theta[3];
	//double phi[3];
	TVector3 ECAL_zaxis(sin(th_earm),0,cos(th_earm));
	TVector3 ECAL_yaxis(0,1,0);
	TVector3 ECAL_xaxis = ECAL_yaxis.Cross(ECAL_zaxis).Unit();

	TVector3 ehat_final_ECAL;
	TVector3 ehat_final_global;
	
	const double th_sbs = 16.9*TMath::DegToRad();
	TVector3 pvect, pvect_SBS;
	TVector3 vertex, vertex_SBS, vertex_ECAL;
	TVector3 SBS_zaxis( -sin(th_sbs), 0, cos(th_sbs) );
	TVector3 SBS_xaxis(0,-1,0);
	TVector3 SBS_yaxis = (SBS_zaxis.Cross(SBS_xaxis)).Unit();
	
	TVector3 sigvtx_global(0.001,0.001,0.007);
	
	TVector3 sigvtx_ECAL( sigvtx_global.Dot( ECAL_xaxis ),
			      sigvtx_global.Dot( ECAL_yaxis ),
			      sigvtx_global.Dot( ECAL_zaxis ) );
	
	for(UInt_t i=0; i< simEvent->fECalClusters.size(); i++){
	  const TSBSECalCluster& eCalHit = simEvent -> fECalClusters[i];
	  // TODO: correct for correlations between "projected" and reconstructed cluster position
	  tempCaloX = eCalHit.GetXPos();
	  tempCaloY = eCalHit.GetYPos();
	  flag = eCalHit.GetDetFlag();
	  union FloatIntUnion {
	    Float_t f;
	    Int_t   i;
	  } datx, daty, datx_2, daty_2;
	  
	  //cout<<"Reconstructed ECal pos X: "<<tempCaloX<<"  Y: "<<tempCaloY<<endl;
	  Int_t crate, slot, chan;
	  
	  //Recalculate all kinematics... a bit dumb
	  E = eCalHit.GetEnergy();
	  if(flag==0){
	    /*
	      datx.f = static_cast<Float_t>(tempCaloX);
	      daty.f = static_cast<Float_t>(tempCaloY);
	      
	      crate = 4;
	      slot  = 0;
	      chan  = 2;
	      
	      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	      }
	      crate = 4;
	      slot  = 0;
	      chan  = 3;
	      if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	      }
	    */
	  }
	  
	  if(flag==10 && E>=fManager->GetCaloThreshold()){
	    datx.f = static_cast<Float_t>(tempCaloX);
	    daty.f = static_cast<Float_t>(tempCaloY);
	    
	    crate = 2;
	    slot  = 0;
	    chan  = 0;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	    }
	    crate = 2;
	    slot  = 1;
	    chan  = 0;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	    }
	  }
	  
	  if(flag==12 && E>=fManager->GetCaloThreshold()){
	    //cout << "Event "<< simEvent->fEvtID << endl;
	    
    	    xfp = yfp = xpfp = ypfp = 0;
	    thetak = phik = 0;

	    
	    xmom = (tempCaloX -  eCalHit.GetXMax())/ECAL_max_cell_size;
	    ymom = (tempCaloY -  eCalHit.GetYMax())/ECAL_max_cell_size;
	    
	    Calc_Shower_Coordinates( xmom, ymom, eCalHit.GetXMax(), eCalHit.GetYMax(), E, z_earm[0], xECal, yECal);
	    //test pour voir...
	    // if(simEvent->fEvtID==2){
	    //   xECal = -0.0898902;
	    //   yECal = -0.666955;
	    // }
	    // if(simEvent->fEvtID==7){
	    //   xECal = 0.292138;
	    //   yECal = 0.371915;
	    // }
	    
	    /*
	    cout << "Event " <<  simEvent->fEvtID << ": " << endl;
	    cout << " X_rec " << tempCaloX << " X_max " << eCalHit.GetXMax() << ", Y_rec " << tempCaloY << ", Y_max " << eCalHit.GetYMax() << endl;
	    cout << " X_mom " << xmom << " X_ECal " << xECal << ", Y_mom " << ymom << ", Y_ECal " << yECal << endl;
	    */
	    
	    //vz = 0;
	    //vz = trk->vertex_target.Z();//
	    vz = find_vertex_bin(trk->vertex_target.Z());//
	    vertex = TVector3(0.0, 0.0, vz);
	    
	    kpx_xc = z_earm[0]*sin(th_earm)+xECal*cos(th_earm);
	    kpy_xc = yECal - yoff_ECAL;
	    kpz_xc = vz+z_earm[0]*cos(th_earm)-xECal*sin(th_earm);
	    kp_xc = sqrt(kpx_xc*kpx_xc+kpy_xc*kpy_xc+kpz_xc*kpz_xc);
	    
	    thetak_xc = acos(kpz_xc/kp_xc);
	    phik_xc = atan2(kpy_xc, kpx_xc);
	    
	    kp_xc = k0*Mp/(k0*(1-cos(thetak_xc))+Mp);
	    pp_xc = sqrt(pow(k0-kp_xc+Mp, 2)-Mp*Mp);
	    thetap_xc = asin(kp_xc*sin(thetak_xc)/pp_xc);
	    
	    vertex_ECAL = TVector3( vertex.Dot(ECAL_xaxis),
				    vertex.Dot(ECAL_yaxis),
				    vertex.Dot(ECAL_zaxis) );
	    
	    xfinal.push_back( vertex_ECAL.X() );
	    yfinal.push_back( vertex_ECAL.Y() );
	    zfinal.push_back( vertex_ECAL.Z() );
	    
	    // wxfinal.push_back( pow( sigvtx_ECAL.X(), -2 ) );
	    // wyfinal.push_back( pow( sigvtx_ECAL.Y(), -2 ) );
	    wxfinal.push_back( pow( 0.001, -2 ) );
	    wyfinal.push_back( pow( 0.001, -2 ) );
	    
	    for(UInt_t j=0; j< simEvent->fScintClusters.size(); j++){
	      const TSBSScintCluster& SciHit = simEvent -> fScintClusters[j];
	      if(SciHit.GetDetFlag()!=31 || SciHit.GetEnergy()<=0)return HED_ERR;
	      //kill the event if we don't have two properly reconstructed CDet hits.
	      //tempScintX = SciHit.GetXPos();
	      //tempScintY = SciHit.GetYPos();
	      //tempPlane = SciHit.GetPlane();
	      //cout << "CDet plane " << SciHit.GetPlane() << ": X = " << SciHit.GetXPos() << ", Y = " << SciHit.GetYPos() << endl;
	      
	      xfinal.push_back( SciHit.GetXPos() );
	      yfinal.push_back( SciHit.GetYPos() - yoff_ECAL );
	      zfinal.push_back( z_earm[SciHit.GetPlane()] );
	      
	      wxfinal.push_back( pow( Lx_scint_CDET, -2 ) );
	      wyfinal.push_back( pow( sigy_CDET, -2 ) );
	    }
	    xfinal.push_back( xECal );
	    yfinal.push_back( yECal - yoff_ECAL );
	    zfinal.push_back( z_earm[0] );
	    
	    wxfinal.push_back( pow( sigx_ECAL, -2 ) );
	    wyfinal.push_back( pow( sigy_ECAL, -2 ) );
	    
	    Fit_3D_Track( xfinal, yfinal, zfinal, wxfinal, wyfinal,
			  eX, eY, eXp, eYp );
	    
	    
	    ehat_final_ECAL = TVector3( eXp, eYp, 1.0 );
	    ehat_final_ECAL = ehat_final_ECAL.Unit();
	    
	    ehat_final_global =
	      ehat_final_ECAL.X() * ECAL_xaxis +
	      ehat_final_ECAL.Y() * ECAL_yaxis +
	      ehat_final_ECAL.Z() * ECAL_zaxis;
			    
			    
			    
	    thetak = acos( ehat_final_global.Z() );
	    phik = atan2( ehat_final_global.Y(), ehat_final_global.X() );
	    
	    // double Eprime_etheta_recon = Ebeam / (1.0 + Ebeam/Mp*(1.-cos(etheta_recon)));
	    // double nu_etheta_recon = Ebeam - Eprime_etheta_recon;
	    
	    // double pp_etheta_recon = sqrt(pow(nu_etheta_recon,2) + 2.*Mp*nu_etheta_recon);
	    
	    // double ptheta_etheta_recon = acos( (Ebeam-Eprime_etheta_recon*cos(etheta_recon))/pp_etheta_recon);
	    
	    phip = phik+TMath::Pi();
	    kp = k0*Mp/(k0*(1-cos(thetak))+Mp);
	    pp = sqrt(pow(k0-kp+Mp, 2)-Mp*Mp);
	    thetap = asin(kp*sin(thetak)/pp);
	    
	    pvect = TVector3(pp*sin(thetap)*cos(phip), pp*sin(thetap)*sin(phip), pp*cos(thetap));
	    pvect_SBS = TVector3( pvect.Dot(SBS_xaxis), pvect.Dot(SBS_yaxis), pvect.Dot(SBS_zaxis) );

	    xptar = pvect_SBS.X()/pvect_SBS.Z();
	    yptar = pvect_SBS.Y()/pvect_SBS.Z();
	    vertex_SBS = TVector3( vertex.Dot(SBS_xaxis), vertex.Dot(SBS_yaxis), vertex.Dot(SBS_zaxis) );
	    ytar = vertex_SBS.Y() - yptar * vertex_SBS.Z();

	    xtar = 0;
	    
	    /*
	    cout << "electron: k = " << kp << " GeV, theta " << thetak << ", phi " << phik << endl;
	    cout << "electron X-check: k = " << kp_xc << " GeV, theta " << thetak_xc << ", phi " << phik_xc << endl;
	    cout << "proton X-check: k = " << pp_xc << " GeV, theta " << thetap_xc << ", phi " << phip_xc << endl;
	    	    
	    cout << "xptar " << xptar << " yptar " << yptar << " ytar " << ytar << " 1/p " << 1.0/pp << " xtar " << xtar << endl;
	    
	    cout << "MC proton vector at target: " << endl;
	    trk->momentum_target.Print();
	    cout << "reconstructed proton vector at target: " << endl;
	    pvect.Print();
	    
	    cout << "MC proton vector in spectrometer: " << endl;
	    trk->fMomentum.Print();
	    cout << "reconstructed proton vector in spectrometer: " << endl;
	    pvect_SBS.Print();
	    */
	    
	    //xtar assumed to be 0... 
	    //Calculate the optics
	    for( int i_par=0; i_par<fManager->GetNOpticsTerms(); i_par++ ){
	      // cout << i_par << " " << fManager->GetOpticsCoeff(0, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(1, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(2, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(3, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(8, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(7, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(6, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(5, i_par) << " "
	      // 	   << fManager->GetOpticsCoeff(4, i_par) << " ";
	      fterm = 
		pow(xptar,fManager->GetOpticsCoeff(8, i_par))*
		pow(yptar,fManager->GetOpticsCoeff(7, i_par))*
		pow(ytar,fManager->GetOpticsCoeff(6, i_par))*
		pow(1.0/pp,fManager->GetOpticsCoeff(5, i_par))*
		pow(xtar,fManager->GetOpticsCoeff(4, i_par));
	      //cout << fterm << endl;
	      
	      xfp += fterm * fManager->GetOpticsCoeff(0, i_par);//b_xfp
	      yfp += fterm * fManager->GetOpticsCoeff(1, i_par);//b_yfp;
	      xpfp += fterm * fManager->GetOpticsCoeff(2, i_par);//b_xpfp;
	      ypfp += fterm * fManager->GetOpticsCoeff(3, i_par);//b_ypfp;
	    }
	    
	    /*
	    cout << "xfp_true " << trk->fOrigin.X()*1.e-3+1.819555*trk->fMomentum.X()/trk->fMomentum.Z()  
		 << " yfp_true " << trk->fOrigin.Y()*1.e-3+1.819555*trk->fMomentum.Y()/trk->fMomentum.Z()
		 << " xpfp_true " << trk->fMomentum.X()/trk->fMomentum.Z() 
		 << " ypfp_true " << trk->fMomentum.Y()/trk->fMomentum.Z() << endl;
	    cout << "xfp_rec " << xfp << " yfp_rec " << yfp << " xpfp_rec " << xpfp << " ypfp_rec " << ypfp << endl;
	    */
	    
	    x_S = xfp;//x_S = xfp-0.015955*xpfp;
	    y_S = yfp;//y_S = yfp-0.015955*ypfp;
	    
	    //filouterie... pour voir...
	    // x_S+= -0.003;
	    // y_S+= -0.003;
	    //cout << "x_S = " << x_S << ", y_S = " << y_S << endl;
	    
	    //x_S = 0.0337-tempCaloX*0.4377;//Ad hoc parameters (probably still much better than above)
	    // y_S = (0.0838-x_S*tempCaloY/tempCaloX)/0.879;
	    datx.f = static_cast<Float_t>(x_S);
	    daty.f = static_cast<Float_t>(y_S);
	    crate = 4;
	    slot  = 0;
	    chan  = 0;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	    }
	    crate = 4;
	    slot  = 0;
	    chan  = 1;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx.i,daty.i) == SD_ERR ){
	      return HED_ERR;
	    }
	    
	    x_S_2 = xfp+0.48191*xpfp;
	    y_S_2 = yfp+0.48191*ypfp;
	    //filouterie... pour voir...
	    // x_S_2+= -0.003;
	    // y_S_2+= -0.003;
	    
	    datx_2.f = static_cast<Float_t>(x_S_2);
	    daty_2.f = static_cast<Float_t>(y_S_2);
	    crate = 4;
	    slot  = 0;
	    chan  = 2;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx_2.i,daty_2.i) == SD_ERR ){
	      return HED_ERR;
	    }
	    crate = 4;
	    slot  = 0;
	    chan  = 3;
	    if( crateslot[idx(crate,slot)]->loadData("adc",chan,datx_2.i,daty_2.i) == SD_ERR ){
	      return HED_ERR;
	    }
	    
	    /*
	    //Recalculate all kinematics... a bit dumb??? =>yes...
	    E = eCalHit.GetEnergy();
	    if(E<fManager->GetCaloThreshold())return HED_ERR;
	    kpx = 5.02*sin(0.506)+tempCaloY*cos(0.506);
	    kpy = -tempCaloX;
	    kpz = 5.02*cos(0.506)-tempCaloY*sin(0.506);
	    
	    kp = sqrt(kpx*kpx+kpy*kpy+kpz*kpz);
	    kpx*= E/kp;
	    kpy*= E/kp;
	    kpz*= E/kp;
	    kp = E;
	    
	    pp = 11.0-kp;
	    thetap = asin(kp/pp*sin(acos(kpz/kp)));
	    phip = atan2(-kpy, -kpx);
	    
	    ppx = pp*cos(thetap);
	    ppy = pp*cos(thetap);
	    ppz = pp*cos(thetap);
	    
	    alpha = atan2(ppx,ppz);
	    tanbeta = ppy/sqrt(ppx*ppx+ppz*ppz);
	    //Double_t v = vz + vx*ppx/ppz;
	    
	    y_S = (3.448)*tan(alpha-16.9);
	    x_S = -(3.448)*tanbeta/cos(alpha-16.9);
	    //cout << x_S << " " << y_S << endl;	    
	    */ 
	  }
	  /*
	  switch(flag){
	  case(0):
	    break;
	  case(10):
	    break;
	  case(12):
	    break;
	  default:
	    break;
	  }
	  */
	  
	}
	/*The previous predicting method below seems to be quite off from the reconstructed
	  ECal position and projected position from the MCTrack, thus switching to reconstructed 
	  ECal position method above

	  // Retrieve last MC track point
	  assert( GetNMCPoints() == 2*fManager->GetNChamber() );
	  MCTrackPoint* pt =
	  static_cast<MCTrackPoint*>( fMCPoints->UncheckedAt(2*fManager->GetNChamber()-1) );
	  assert( pt );
	  const TVector3& pos = pt->fMCPoint;
	  TVector3 dir = pt->fMCP.Unit();
	  if( fManager->GetCaloZ() <= pos.Z() ) {
	  Error( here, "Calorimeter z = %lf less than z of last GEM plane = %lf. "
	  "Set correct value with SetCaloZ() or turn off calo emulation.",
	  fManager->GetCaloZ(), pos.Z() );
	  return HED_FATAL;
	  }
	  if( TMath::Abs(dir.Z()) < 1e-6 ) {
	  Error( here, "Illegal primary track direction (%lf,%lf,%lf). "
	  "Should never happen. Call expert.", dir.X(), dir.Y(), dir.Z() );
	  return HED_ERR;
	  }
	  dir *= 1.0/dir.Z();  // Make dir a transport vector
	  TVector3 hitpos = pos + (fManager->GetCaloZ()-pos.Z()) * dir;

	  // Smear the position with the given resolution
	  // Assumes z-axis normal to calorimeter plane. Otherwise we would have to
	  // get the plane's fXax and fYax
	  TVector3 res( gRandom->Gaus(0.0, fManager->GetCaloRes()),
	  gRandom->Gaus(0.0, fManager->GetCaloRes()), 0.0 );
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
	*/
    }
  }

  // DEBUG:
  //cout << "SimDecoder: nTracks = " << GetNMCTracks() << endl;
  //fMCTracks.Print();
  // cout<<" \n\n###\n\n"<<endl;
  return HED_OK;
}

//-----------------------------------------------------------------------------
TSBSSimGEMHit::TSBSSimGEMHit( const TSBSSimEvent::GEMCluster& c )
  : fID(c.fID), fSector(c.fSector), fPlane(c.fPlane),
    //fRealSector(c.fRealSector), 
    fModule(c.fModule), fSource(c.fSource), fType(c.fType),
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
    Double_t dz = c.fMCpos.Z() - fManager->GetZ0();
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
// stuff for ECal analysis
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void TSBSSimDecoder::Calc_Shower_Coordinates( double xmom, double ymom, double xmax, double ymax, double Eclust, double Rcal, double &xclust, double &yclust )//, double &xf, double &yf
{
  //calculate longitudinal depth of max. shower energy deposition
  double tmax = TMath::Max(0.0, X0_ECAL * (log( Eclust/Ec_ECAL ) - 0.5) );

  int binx = int( (100.0*xmax-profx_xmin)/(profx_xmax-profx_xmin)*profx_nbins );

  shower_profile_t proftemp;
  if( binx >= 0 && binx < profx_nbins ){
    proftemp = profx[binx];
  } else {
    proftemp = profxdefault;
  }
  
  int binxmom = int( (xmom - proftemp.xmin)/(proftemp.xmax-proftemp.xmin) * proftemp.nbins_mom );

  //set some sensible default value:
  xclust = xmax + xmom*ECAL_max_cell_size; 
  
  //do linear interpolation within the bin:
  if( binxmom < 0 ){
    xclust = xmax - 0.5 * ECAL_max_cell_size;
    //xf = -0.5*ECAL_max_cell_size;
  } else if( binxmom < proftemp.nbins_mom ){
    double binwidth = (proftemp.xmax-proftemp.xmin)/double(proftemp.nbins_mom);
    double flo = proftemp.frac[binxmom];
    double fhi = proftemp.frac[binxmom+1];
    
    double xlo = proftemp.xmin + binxmom * binwidth;
    //double xhi = xlo + binwidth;
    //linear interpolation within the bin:
    double f = flo + (fhi-flo) * (xmom-xlo)/binwidth;
    xclust = xmax + (f - 0.5)*ECAL_max_cell_size;

    // xf = xclust-xmax;
    //cout << "xmom, binxmom, f = " << xmom << ", " << binxmom << ", " << f << endl;
  } else {
    xclust = xmax + 0.5 * ECAL_max_cell_size;
    
  }

  //before incident-angle correction, record x position within cell:
  //xf = xclust-xmax;
  
  int biny = int( (100.0*ymax-profy_ymin)/(profy_ymax-profy_ymin)*profy_nbins );

  if( biny >= 0 && biny < profy_nbins ){
    proftemp = profy[biny];
  } else {
    proftemp = profydefault;
  }

  int binymom = int( (ymom - proftemp.xmin)/(proftemp.xmax-proftemp.xmin) * proftemp.nbins_mom );

  yclust = ymax + ymom*ECAL_max_cell_size;

  if( binymom < 0 ){
    yclust = ymax - 0.5 * ECAL_max_cell_size;
  } else if( binymom < proftemp.nbins_mom ){
    double binwidth = (proftemp.xmax-proftemp.xmin)/double(proftemp.nbins_mom);
    double flo = proftemp.frac[binymom];
    double fhi = proftemp.frac[binymom+1];
    
    double ylo = proftemp.xmin + binymom * binwidth;
    //double xhi = xlo + binwidth;
    //linear interpolation within the bin:
    double f = flo + (fhi-flo) * (ymom-ylo)/binwidth;
    yclust = ymax + (f - 0.5)*ECAL_max_cell_size;

    //cout << "ymom, biny, binymom, f = " << ymom << ", " << biny << "," <<  binymom << ", " << f << endl;
  } else {
    yclust = ymax + 0.5 * ECAL_max_cell_size;
  }

  //yf = yclust-ymax;
  
  //Apply incident-angle correction under the assumption that track starts at (x,y,z) = (0,0,0)
  double xptemp = xclust/Rcal;
  double yptemp = yclust/Rcal;

  double dz = tmax / sqrt(1.0 + pow(xptemp,2) + pow(yptemp,2) );

  xclust -= dz * xptemp;
  yclust -= dz * yptemp;
  
  return;
}


bool TSBSSimDecoder::load_shower_profiles( const char *filename ){
  ifstream fshower(filename);

  TString currentfile;
  
  currentfile.ReadFile( fshower );
  
  //Now the whole file is loaded into the TString:

  int istart, istop, length;

  std::istringstream sstream_temp;
  
  TString stemp,substrtemp;
  stemp = "ECAL_shower_profile_x_default";

  //Read to the next newline character:
  istart = currentfile.Index(stemp);
  istop = currentfile.Index( '\n', istart ); 

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  //cout << "'" << substrtemp << "'" << endl;

  substrtemp.ReplaceAll(stemp,"");

  //cout << "'" << substrtemp << "'" << endl;

  sstream_temp = std::istringstream(substrtemp.Data());

  //extract number of bins and min and max moment value for default x profile:
  sstream_temp >> profxdefault.nbins_mom >> profxdefault.xmin >> profxdefault.xmax;

  istart = istop;
  istop = currentfile.Index( "ECAL_shower_profile", istart ); //returns position of next header

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  //Grab profile data:
  substrtemp = TString(currentfile(istart+1,istop-istart-1));

  //cout << "'" << substrtemp << "'" << endl;
  //for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
  sstream_temp = std::istringstream(substrtemp.Data());

  profxdefault.frac.resize(profxdefault.nbins_mom+1);

  //extract bin contents from profile data:
  for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
    sstream_temp >> profxdefault.frac[bin];
    //cout << "bin, frac = " << bin << ", " << profxdefault.frac[bin] << endl;
  }

  //Do the same for the default y profile:
  stemp = "ECAL_shower_profile_y_default";

  
  
  istart = currentfile.Index(stemp);
  istop = currentfile.Index( '\n', istart ); 

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  //cout << "'" << substrtemp << "'" << endl;

  substrtemp.ReplaceAll(stemp,"");

  //cout << "'" << substrtemp << "'" << endl;

  sstream_temp = std::istringstream(substrtemp.Data());

  sstream_temp >> profydefault.nbins_mom >> profydefault.xmin >> profydefault.xmax;

  istart = istop;
  istop = currentfile.Index( "ECAL_shower_profile", istart ); //returns position of next header

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart+1,istop-istart-1));

  //cout << "'" << substrtemp << "'" << endl;
  //for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
  sstream_temp = std::istringstream(substrtemp.Data());

  profydefault.frac.resize(profydefault.nbins_mom+1);
    
  for( int bin=0; bin<=profydefault.nbins_mom; bin++ ){
    sstream_temp >> profydefault.frac[bin];
    //cout << "bin, frac = " << bin << ", " << profxdefault.frac[bin] << endl;
  }

  //Now read x-dependent x profiles and y-dependent y profiles:
  stemp = "ECAL_shower_profiles_x";

  istart = currentfile.Index( stemp );
  istop = currentfile.Index('\n',istart);

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  substrtemp.ReplaceAll(stemp,"");

  sstream_temp = std::istringstream(substrtemp.Data());

  sstream_temp >> profx_nbins >> profx_xmin >> profx_xmax;

  profx.resize( profx_nbins );
  
  for( int xbin=0; xbin<profx_nbins; xbin++ ){
    stemp.Form("ECAL_shower_profile_x_bin %d",xbin+1);
    istart = currentfile.Index(stemp);
    istop = currentfile.Index('\n',istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart,istop-istart));
    substrtemp.ReplaceAll(stemp,"");

    sstream_temp = std::istringstream(substrtemp.Data());

    sstream_temp >> profx[xbin].nbins_mom >> profx[xbin].xmin >> profx[xbin].xmax;

    istart = istop;
    istop = currentfile.Index("ECAL",istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart+1,istop-istart-1));

    sstream_temp = std::istringstream(substrtemp.Data());

    profx[xbin].frac.resize(profx[xbin].nbins_mom+1);

    for( int bin=0; bin<=profx[xbin].nbins_mom; bin++ ){
      sstream_temp >> profx[xbin].frac[bin];

      //cout << "xbin, bin, frac = " << xbin << ", " << bin << ", " << profx[xbin].frac[bin] << endl;
    }
  }

  ///y-dependent y profiles:
  stemp = "ECAL_shower_profiles_y";

  istart = currentfile.Index( stemp );
  istop = currentfile.Index('\n',istart);

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  substrtemp.ReplaceAll(stemp,"");

  sstream_temp = std::istringstream(substrtemp.Data());

  sstream_temp >> profy_nbins >> profy_ymin >> profy_ymax;

  profy.resize( profy_nbins );
  
  for( int ybin=0; ybin<profy_nbins; ybin++ ){
    stemp.Form("ECAL_shower_profile_y_bin %d",ybin+1);
    istart = currentfile.Index(stemp);
    istop = currentfile.Index('\n',istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart,istop-istart));
    substrtemp.ReplaceAll(stemp,"");

    sstream_temp = std::istringstream(substrtemp.Data());

    sstream_temp >> profy[ybin].nbins_mom >> profy[ybin].xmin >> profy[ybin].xmax;

    istart = istop;
    istop = currentfile.Index("ECAL",istart);

    //cout << "ybin, istop = " << ybin << ", " << istop << endl;

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart+1,istop-istart-1));

    sstream_temp = std::istringstream(substrtemp.Data());

    profy[ybin].frac.resize(profy[ybin].nbins_mom+1);

    for( int bin=0; bin<=profy[ybin].nbins_mom; bin++ ){
      sstream_temp >> profy[ybin].frac[bin];

      //cout << "ybin, bin, frac = " << ybin << ", " << bin << ", " << profy[ybin].frac[bin] << endl;
    }
  }
  
  return true;
}

/*
void SBS_tgt_reconstruct( double xfp, double yfp, double xpfp, double ypfp, double xtar, double &xptar, double &yptar, double &ytar, double &p ){
  double sum_xptar=0.0;
  double sum_yptar=0.0;
  double sum_ytar=0.0;
  double sum_ptheta=0.0;

  for( int i=0; i<nterms_optics; i++ ){
    double term = pow(xfp,Cexpon[i][0]) * pow(yfp,Cexpon[i][1]) * pow(xpfp,Cexpon[i][2])
      * pow(ypfp,Cexpon[i][3])*pow(xtar,Cexpon[i][4]);
    sum_xptar += Cxptar[i]*term;
    sum_yptar += Cyptar[i]*term;
    sum_ytar += Cytar[i]*term;
    sum_ptheta += Cpthetabend[i]*term;
  }

  TVector3 zaxis_FP(-sin(SBS_tracker_pitch),0,cos(SBS_tracker_pitch));
  TVector3 yaxis_FP(0,1,0);
  TVector3 xaxis_FP = yaxis_FP.Cross(zaxis_FP).Unit();

  TVector3 nhat_fp(xpfp,ypfp,1.0);
  nhat_fp = nhat_fp.Unit();

  TVector3 nhat_fp_global = nhat_fp.X() * xaxis_FP + nhat_fp.Y() * yaxis_FP + nhat_fp.Z() * zaxis_FP;
  TVector3 nhat_tgt(xptar,yptar,1.0);
  nhat_tgt = nhat_tgt.Unit();

  double thetabend = acos( nhat_fp_global.Dot( nhat_tgt ) );

  p = sum_ptheta/thetabend;

  xptar = sum_xptar;
  yptar = sum_yptar;
  ytar = sum_ytar;
}

void SBS_fp_reconstruct( double xtar, double ytar, double xptar, double yptar, double p, double &xfp, double &yfp, double &xpfp, double &ypfp ){
  double sum_xfp = 0.0;
  double sum_yfp = 0.0;
  double sum_xpfp = 0.0;
  double sum_ypfp = 0.0;

  for( int i=0; i<nterms_foptics; i++ ){
    double term = pow(xptar, Cfexpon[i][0]) * pow(yptar,Cfexpon[i][1]) * pow(ytar,Cfexpon[i][2]) * pow(1.0/p,Cfexpon[i][3]) * pow(xtar, Cfexpon[i][4] );
    sum_xfp += Cxfp[i]*term;
    sum_yfp += Cyfp[i]*term;
    sum_xpfp += Cxpfp[i]*term;
    sum_ypfp += Cypfp[i]*term;
  }
  xfp = sum_xfp;
  yfp = sum_yfp;
  xpfp = sum_xpfp;
  ypfp = sum_ypfp;
}

void find_ECAL_clusters( gep_tree_with_spin *T, int &nclust, vector<double> &xclust, vector<double> &yclust, vector<double> &Eclust, vector<int> &nhitclust, vector<vector<int> > &hitlist_clust ){ //This method is only intended to group ECAL hits together in clusters and do "crude" coordinate reconstuction (shower center-of-gravity)

  nclust = 0;
  xclust.clear();
  yclust.clear();
  Eclust.clear();
  nhitclust.clear();
  //nxclust.clear();
  //nyclust.clear();
  hitlist_clust.clear();

  double Eclustmin = 0.5*T->ev_ep; //50% of elastic

  bool foundclust = false;

  int nhitstot = T->Earm_ECalTF1_hit_nhits;

  //cout << "ECAL cluster finding, nhits = " << nhitstot << endl;

  set<int> unused_hits;
  //vector<bool> hitused(nhitstot);
  vector<double> Ehit_recon(nhitstot);
  
  for( int hit=0; hit<nhitstot; hit++ ){
    //hitused[hit] = false;
    unused_hits.insert(hit);

    double Etemp = (*(T->Earm_ECalTF1_hit_sumedep))[hit];
    double nphemean = Etemp*ECAL_phe_per_GeV;
    
    double nphe = random_generator->PoissonD(nphemean);

    Ehit_recon[hit] = nphe/ECAL_phe_per_GeV;

    //cout << "(ihit, Etrue, Erecon)=(" << hit << ", " << Etemp << ", " << Ehit_recon[hit] << ")" << endl;
    
  }
  
  //repeat cluster search until no new clusters found:
  do {
    foundclust = false;
    //Step 1: loop over all unused hits; find maximum.
    double Ehitmax = 0.0;
    int ihitmax = -1; //position in hit array of hit with largest energy.

    for( set<int>::iterator hit=unused_hits.begin(); hit != unused_hits.end(); ++hit ){
      int ihit = *hit;

      double Ehit = Ehit_recon[ihit];

      ihitmax = Ehit > Ehitmax ? ihit : ihitmax;
      Ehitmax = Ehit > Ehitmax ? Ehit : Ehitmax;
      
    }

    //cout << "(ihitmax, Ehitmax)=(" << ihitmax << ", " << Ehitmax << ")" << endl;

    if( ihitmax < 0 ) {
      foundclust = false;
    } else { //found a new maximum: start a cluster around this maximum
      //unused_hits.erase( ihitmax );

      int ncellclust_temp = 1;
      vector<int> hitlist_temp;
      hitlist_temp.push_back( ihitmax );

      double sum_logE = log(Ehitmax);
      double Eclust_temp = Ehitmax;
      double xclust_temp = (*(T->Earm_ECalTF1_hit_xcell))[ihitmax]*Ehitmax;
      double yclust_temp = (*(T->Earm_ECalTF1_hit_ycell))[ihitmax]*Ehitmax;
      //int nxclust_temp = 1;
      //int nyclust_temp = 1;
      
      int icelltemp = 0;
      //Next, we want a do while loop over all hits in the cluster to add "nearest neighbors"
      while ( icelltemp < hitlist_temp.size() ) {

	// unused_hits.erase( hitlist_temp[icelltemp] );
	// at the beginning of each iteration of finding nearest-neighbors, mark all hits
	// associated with this cluster as used, starting with icelltemp
	// everything preceding icelltemp should have already been removed!
	for( int jhit=icelltemp; jhit<hitlist_temp.size(); jhit++ ){
	  unused_hits.erase( hitlist_temp[jhit] ); 
	}
	//grab information about the current cell:
	int rowtemp = (*(T->Earm_ECalTF1_hit_row))[hitlist_temp[icelltemp]];
	int coltemp = (*(T->Earm_ECalTF1_hit_col))[hitlist_temp[icelltemp]];
	double xcelltemp = (*(T->Earm_ECalTF1_hit_xcell))[hitlist_temp[icelltemp]];
	double ycelltemp = (*(T->Earm_ECalTF1_hit_ycell))[hitlist_temp[icelltemp]];
	double Ecelltemp = Ehit_recon[hitlist_temp[icelltemp]];
	//loop over the list of unused hits, adding any nearest-neighbor hits found:
	//note that when we start this loop, the central maximum of the cluster has already been removed from the
	//set of unused hits!

	//Do NOT modify the set while iterating through it!
	for( set<int>::iterator hit=unused_hits.begin(); hit != unused_hits.end(); ++hit ){
	  int ihit = *hit;

	  //check for same row and column +/- 1, or row +/- 1 and (x - xcell)

	  double xhit = (*(T->Earm_ECalTF1_hit_xcell))[ihit];
	  double yhit = (*(T->Earm_ECalTF1_hit_ycell))[ihit];
	  double Ehit = Ehit_recon[ihit];

	  if( pow(xhit-xcelltemp,2)+pow(yhit-ycelltemp,2)<=pow(1.5*ECAL_max_cell_size,2) ){ //nearest-neighbor:
	    hitlist_temp.push_back( ihit );
	    xclust_temp += xhit*Ehit;
	    yclust_temp += yhit*Ehit;
	    Eclust_temp += Ehit;
	    sum_logE += log(Ehit);
	  }
	}

	icelltemp++; //after first iteration, which always happens, icelltemp = 1.
	// if any nearest-neighbors were found, hitlist_temp.size() > 1
	// we keep looping over all hits in the list of hits added to this cluster until
	// no more unused nearest-neighbor hits are found!
	// As long as at least one new nearest neighbor was found, we keep going!
	// As soon as we don't find any new neighbors, the condition below fails and we exit the loop!
	
      } 

      //Add this cluster to the output arrays:

      if( Eclust_temp >= Eclustmin ){
      
	xclust.push_back( xclust_temp/Eclust_temp );
	yclust.push_back( yclust_temp/Eclust_temp );
	Eclust.push_back( Eclust_temp );
	nhitclust.push_back( hitlist_temp.size() );
	hitlist_clust.push_back( hitlist_temp );
      
	foundclust = true;

	// cout << "found cluster E, x, y = " << Eclust_temp << ", " << xclust_temp/Eclust_temp
	//      << ", " << yclust_temp/Eclust_temp << endl;
      }
    }
  } while( foundclust );
  
  nclust = xclust.size();

}
*/

//-----------------------------------------------------------------------------
void TSBSSimDecoder::Fit_3D_Track( vector<double> xpoints, vector<double> ypoints, vector<double> zpoints, vector<double> wx, vector<double> wy,
		   double &X, double &Y, double &Xp, double &Yp ){

  //Setting up the matrices for the fit:
  
  TMatrixD A(4,4);
  TVectorD b(4);

  for( int i=0; i<4; i++ ){
    for( int j=0; j<4; j++ ){
      A(i,j) = 0.0;
    }
    b(i) = 0.0;
  }

  int npoints = xpoints.size();

  if( npoints<2 ) return;

  int ndf=0;

  //For a 3D fit to a straight-line:
  // chi^2 = sum_i wxi * (xi- (X + Xp*zi))^2 + wyi*(y - (Y+Yp*zi))^2
  // dchi^2/dX = -2 * (xi - (X+Xp*zi))* wxi = 0
  // dchi^2/dY = -2 * (yi - (Y+Yp*zi))* wyi = 0
  // dchi^2/dXp = -2 * (xi - (X+Xp*zi))*zi * wxi = 0
  // dchi^2/dYp = -2 * (yi - (Y+Yp*zi))*zi * wyi = 0
  for( int i=0; i<npoints; i++ ){
    b(0) += wx[i]*xpoints[i];
    b(1) += wy[i]*ypoints[i];
    b(2) += wx[i]*xpoints[i]*zpoints[i];
    b(3) += wy[i]*ypoints[i]*zpoints[i];

    A(0,0) += wx[i];
    A(0,1) += 0.0;
    A(0,2) += wx[i]*zpoints[i];
    A(0,3) += 0.0;

    A(1,0) += 0.0;
    A(1,1) += wy[i];
    A(1,2) += 0.0;
    A(1,3) += wy[i]*zpoints[i];

    A(2,0) += wx[i]*zpoints[i];
    A(2,1) += 0.0;
    A(2,2) += wx[i]*pow(zpoints[i],2);
    A(2,3) += 0.0;

    A(3,0) += 0.0;
    A(3,1) += wy[i]*zpoints[i];
    A(3,2) += 0.0;
    A(3,3) += wy[i]*pow(zpoints[i],2);
    
  }

  TVectorD solution = A.Invert() * b;

  X = solution(0);
  Y = solution(1);
  Xp = solution(2);
  Yp = solution(3);
}

double TSBSSimDecoder::find_vertex_bin(double vz_true)
{
  int Ntgtbins = TMath::CeilNint(Ltgt/vz_bin_width);
  double minbin = -vz_bin_width*Ntgtbins/2.0;
  for(int i = 0; i<Ntgtbins; i++){
    if( minbin+i*vz_bin_width<=vz_true && vz_true<=minbin+(i+1)*vz_bin_width )
      return(minbin+(i+0.5)*vz_bin_width);
  }
  cout << "haven't found a vertex bin for value " << vz_true << endl;
  return vz_true;
}

//-----------------------------------------------------------------------------
ClassImp(TSBSSimDecoder)
ClassImp(TSBSSimGEMHit)
ClassImp(TSBSSimBackTrack)
