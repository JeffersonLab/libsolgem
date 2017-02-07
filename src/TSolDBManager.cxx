#include "TSolDBManager.h"
#include "TSolSimDecoder.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"

TSolDBManager * TSolDBManager::fManager = NULL;

TSolDBManager::TSolDBManager() 
: fErrID(-999), fErrVal(-999.)
{
}
//______________________________________________________________
TSolDBManager::~TSolDBManager()
{
}
//______________________________________________________________
void TSolDBManager::LoadGeneralInfo(const string& fileName)
{
    ifstream input(fileName.c_str());
    if (!input.is_open()){
        cout<<"cannot find general information file "<<fileName
            <<". Exiting the program"<<endl;
        exit(0);
    }
    const string prefix = "generalinfo.";
    DBRequest request[] = {
        {"do_map_sector",       &fDoMapSector         , kInt,    0, 1},
        {"self_define_sector",  &fDoSelfDefinedSector , kInt,    0, 1},
        {"sector_mapped",       &fMappedSector        , kInt,    0, 1},
	{"ngsector",            &fNgSector            , kInt,    0, 1},
	{"nchamber1",           &fNChamber1           , kInt,    0, 1},
        {"nsector1",            &fNSector1            , kInt,    0, 1},
        {"nchamber2",           &fNChamber2           , kInt,    0, 1},
        {"nsector2",            &fNSector2            , kInt,    0, 1},
        {"nreadout",            &fNReadOut            , kInt,    0, 1},
        {"gem_drift_id",        &fGEMDriftID          , kInt,    0, 1},
        {"gem_copper_front_id", &fGEMCopperFrontID    , kInt,    0, 1},
        {"gem_copper_back_id",  &fGEMCopperBackID     , kInt,    0, 1},
        {"gem_strip_id",        &fGEMStripID          , kInt,    0, 1},
        {"faec_id",             &fFAECID              , kInt,    0, 1},
        {"laec_id",             &fLAECID              , kInt,    0, 1},
        {"nsignal",             &fNSigParticle        , kInt,    0, 1},
        {"chan_per_slot",       &fChanPerSlot         , kInt,    0, 1},
        {"modules_per_readout", &fModulesPerReadOut   , kInt,    0, 1},
        {"self_defined_sector", &fDoSelfDefinedSector , kInt,    0, 1},
	{"g4sbs_detectortype",  &fg4sbsDetectorType   , kInt,    0, 1},
	{"g4sbs_z_specoffset",  &fg4sbsZSpecOffset    , kDouble, 0, 1},
	{"z0",                  &fgZ0                 , kDouble, 0, 1},
	{"calo_z",              &fgCaloZ              , kDouble, 0, 1},
	{"calo_res",            &fgCaloRes            , kDouble, 0, 1},
	{"docalo",              &fgDoCalo             , kInt,    0, 1},
        { 0 }
    };
    int pid, tid;
    DBRequest signalRequest[] = {
        {"pid",                 &pid,                   kInt, 0, 1},
        {"tid",                 &tid,                   kInt, 0, 1},
        { 0 }
    };
    int err = LoadDB( input, request,  prefix);
   
    if( err ) exit(2); 
    
    for (int i=0; i<fNSigParticle; i++){
        ostringstream signal_prefix(prefix, ios_base::ate);
        signal_prefix<<"signal"<<i+1<<".";
        
        err = LoadDB(input, signalRequest, signal_prefix.str());
        
        fSigPID.push_back(pid);
        fSigTID.push_back(tid);
	
	if( err ) exit(2); 
    }
    
    for (int i=0; i<GetNChamber(); i++){
         vector<GeoInfo> thisInfo;
         thisInfo.clear();
         fGeoInfo[i] = thisInfo;
    }
    
    fModulesPerChamber = fModulesPerReadOut * fNReadOut;
    
    fChambersPerCrate  = 
      (TSolSimDecoder::GetMAXSLOT()/fModulesPerChamber/(fNChamber1+fNChamber2)) * (fNChamber1+fNChamber2);
}

//______________________________________________________________
void TSolDBManager::LoadGeoInfo(const string& prefix)
{
  const string& fileName = "db_"+prefix+".dat";
    
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find geometry file "<<fileName
	<<". Exiting the program"<<endl;
    exit(0);
  }
  
  //const string prefix = "gemc.";
  
  GeoInfo thisGeo;
  
  DBRequest request[] = {
    {"dmag",        &thisGeo.dmag,         kDouble, 0, 1},
    {"d0",          &thisGeo.d0,           kDouble, 0, 1},
    {"xoffset",     &thisGeo.xoffset,      kDouble, 0, 1},
    {"dx",          &thisGeo.dx,           kDouble, 0, 1},
    {"dy",          &thisGeo.dy,           kDouble, 0, 1},
    {"thetaH",      &thisGeo.thetaH,       kDouble, 0, 1},
    {"thetaV",      &thisGeo.thetaV,       kDouble, 0, 1},
    {"depth",       &thisGeo.depth,        kDouble, 0, 1},
    { 0 }
  };
  
  DBRequest plane_request[] = {
    { "x.stripangle",     &thisGeo.stripangle_u,   kDouble, 0, 1},
    { "x.pitch",          &thisGeo.pitch_u,        kDouble, 0, 1},
    { "y.stripangle",     &thisGeo.stripangle_v,   kDouble, 0, 1},
    { "y.pitch",          &thisGeo.pitch_v,        kDouble, 0, 1},
    { 0 }
  };
  
  for (int i=0; i<fNChamber2; i++){
    map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
    if (it == fGeoInfo.end()) { cout<<"unexpected chamber id "<<i<<endl; }
    
    for (int j=0; j<fNSector2; j++){
      ostringstream sector_prefix(prefix, ios_base::ate);
      int idx = i*fNSector2 + j;
      sector_prefix<<".gem"<<idx<<".";
      
      int err = LoadDB(input, request, sector_prefix.str());
      if( err ) exit(2);
      
      sector_prefix<<"gem"<<idx;
      err = LoadDB(input, plane_request, sector_prefix.str());
      if (err) exit(2);
      
      fGeoInfo[i].push_back(thisGeo);
    }
  }
  for (int i=0; i<fNChamber1; i++){
    map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
    if (it == fGeoInfo.end()) { cout<<"unexpected chamber id "<<i<<endl; }
    
    for (int j=0; j<fNSector1; j++){
      ostringstream sector_prefix(prefix, ios_base::ate);
      int idx = fNChamber2*fNSector2 +i*fNSector1 + j;
      sector_prefix<<".gem"<<idx<<".";
      
      int err = LoadDB(input, request, sector_prefix.str());
      if( err ) exit(2);
      
      sector_prefix<<"gem"<<idx;
      err = LoadDB(input, plane_request, sector_prefix.str());
      if (err) exit(2);
      
      fGeoInfo[fNChamber2+i].push_back(thisGeo);
    }
  }
  
  // cout << "fGeo size: " << fGeoInfo.size() << endl;
  // for(uint i = 0; i<fGeoInfo.size(); i++){
  //   cout << fGeoInfo[i].size() << endl;
  // }
}

//______________________________________________________________
string TSolDBManager::FindKey( ifstream& inp, const string& key )
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}
//_________________________________________________________________________
bool TSolDBManager::CheckIndex(int i, int j, int k)
{
    if (i >= fNChamber1+fNChamber2 || i < 0){
        cout<<"invalid chamber ID requested: "<<i<<endl;
        return false;
    }
    else if(i<fNChamber2 && j>=fNSector2){
      cout<<"invalid sector id requested: "<<j<<endl;
      return false;
    }
    else if(i>=fNChamber2 && j>=fNSector1){
      cout<<"invalid sector id requested: "<<j<<endl;
      return false;
    }
    else if (k >= fNReadOut || k < 0){
        cout<<"invalid readout id requested: "<<k<<endl;
    }
    return true;
}
//_________________________________________________________________
int TSolDBManager::LoadDB( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = FindKey(inp,key);
    if( !val.empty() ) {
      istringstream sv(val);
      switch(item->type){
        case kDouble:
          sv >> *((double*)item->var);
          break;
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        default:
          return 1;
        break;
      }
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}
//_____________________________________________________________________
const int & TSolDBManager::GetSigPID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigPID[i];
}
//______________________________________________________________________
const int & TSolDBManager::GetSigTID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigTID[i];
}

//const 
int TSolDBManager::GetNChamber()
{ 
  const int NchamberTot = fNChamber1+fNChamber2;
  return NchamberTot;
}

//______________________________________________________________________
const double & TSolDBManager::GetDMag(int i, int j)
{
  // cout << "D0: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << fGeoInfo[i].size() << endl;
  return fGeoInfo[i].at(j).dmag;
}
//______________________________________________________________________
const double & TSolDBManager::GetD0(int i, int j)
{
  // cout << "D0: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << fGeoInfo[i].size() << endl;
  return fGeoInfo[i].at(j).d0;
}
//______________________________________________________________________
const double & TSolDBManager::GetXOffset(int i, int j)
{
  // cout << "XOff: i, j " << i << " " << j << " Geo size, Geo[i] size "  << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << fGeoInfo[i].size() << endl;
  return fGeoInfo[i].at(j).xoffset;
}
//______________________________________________________________________
const double & TSolDBManager::GetDX(int i, int j)
{
  //cout << "DX: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size();
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << " " << fGeoInfo[i].size() << endl;
  return fGeoInfo[i].at(j).dx;
}
//______________________________________________________________________
const double & TSolDBManager::GetDY(int i, int j)
{
  // cout << "DY: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size();
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << " " << fGeoInfo[i].size() << endl;
  return fGeoInfo[i].at(j).dy;
}
//______________________________________________________________________
const double & TSolDBManager::GetThetaH(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).thetaH;
}
//______________________________________________________________________
const double & TSolDBManager::GetThetaV(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).thetaV;
}
//_________________________________________________________________________
const double & TSolDBManager::GetStripAngle(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[i].at(j).stripangle_u;
    else return fGeoInfo[i].at(j).stripangle_u;
}
//_________________________________________________________________________
const double & TSolDBManager::GetPitch(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[i].at(j).pitch_u;
    else return fGeoInfo[i].at(j).pitch_u;
}

//__________________________________________________________________________
int TSolDBManager::GetSectorIDFromPos(int ichamber, double x, double y)
{
  if (!CheckIndex(ichamber)) return fErrVal;

  double sector = -1;
  for(int k = 0; k<fGeoInfo[ichamber].size(); k++){
    if(fGeoInfo[ichamber].at(k).xoffset-fGeoInfo[ichamber].at(k).dx/2.0<x && 
       x<fGeoInfo[ichamber].at(k).xoffset+fGeoInfo[ichamber].at(k).dx/2.0){
      sector = k;
    }
  }
  return sector;
}








