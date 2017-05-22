#include "TSolDBManager.h"
#include "TSBSSimDecoder.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"

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
  // EFuchey: 2017/02/09: Since this date, the reading, digitization, etc. of the data 
  // has been splitted for Forward Tracker and Focal Plane Polarimeter.
  // This simplifies the number of parameters for General info.
  // On the other hand I also changed the ordering of the GEMs reading: 
  // It used to be: chamber 0: loop on sectors, chmaber 1: loop on sectors...
  // now it is: sector 0: loop on chambers, sector 1: loop on chambers, etc...
  // That implied some modifications of the functions retrieving 
  // a plane parameter as a function of the plane and sector number, 
  // as well as the function which retrives the sector number from the position
  
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
	{"nchamber",            &fNChamber            , kInt,    0, 1},
        {"nsector",             &fNSector             , kInt,    0, 1},
	// // see comment TSolDBManager.h, l. 92-94.
	// {"nchamber2",            &fNChamber2            , kInt,    0, 1},
        // {"nsector2",             &fNSector2             , kInt,    0, 1},
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
    
    for (int i=0; i<GetNSector(); i++){
      vector<GeoInfo> thisInfo;
      thisInfo.clear();
      fGeoInfo[i] = thisInfo;
    }
    
    fModulesPerChamber = fModulesPerReadOut * fNReadOut;
    
    fChambersPerCrate = 
      (TSBSSimDecoder::GetMAXSLOT()/fModulesPerChamber/fNChamber) * fNChamber;
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
  
  for (int i=0; i<fNSector; i++){
    map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
    if (it == fGeoInfo.end()) { cout<<"unexpected chamber id "<<i<<endl; }
    
    for (int j=0; j<fNChamber; j++){
      ostringstream sector_prefix(prefix, ios_base::ate);
      int idx = i*fNChamber + j;
      sector_prefix<<".gem"<<idx<<".";
      
      int err = LoadDB(input, request, sector_prefix.str());
      if( err ) exit(2);
      
      sector_prefix<<"gem"<<idx;
      err = LoadDB(input, plane_request, sector_prefix.str());
      if (err) exit(2);
      
      fGeoInfo[i].push_back(thisGeo);
    }
  }
  /*
  for (int i=0; i<fNSector2; i++){
    map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
    if (it == fGeoInfo.end()) { cout<<"unexpected chamber id "<<i<<endl; }
    
    for (int j=0; j<fNChamber2; j++){
      ostringstream sector_prefix(prefix, ios_base::ate);
      int idx = i*fNChamber2 + j;
      sector_prefix<<".gem"<<idx<<".";
      
      int err = LoadDB(input, request, sector_prefix.str());
      if( err ) exit(2);
      
      sector_prefix<<"gem"<<idx;
      err = LoadDB(input, plane_request, sector_prefix.str());
      if (err) exit(2);
      
      fGeoInfo[i].push_back(thisGeo);
    }
  }
  */
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
    if (i >= fNChamber || i < 0){
        cout<<"invalid chamber ID requested: "<<i<<endl;
        return false;
    }
    else if(j>=fNSector){
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

//______________________________________________________________________
const double & TSolDBManager::GetDMag(int i, int j)
{
  // cout << "D0: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << fGeoInfo[j].size() << endl;
  return fGeoInfo[j].at(i).dmag;
}
//______________________________________________________________________
const double & TSolDBManager::GetD0(int i, int j)
{
  //cout << "D0: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  //cout << fGeoInfo[j].size() << endl;
  return fGeoInfo[j].at(i).d0;
}
//______________________________________________________________________
const double & TSolDBManager::GetXOffset(int i, int j)
{
  //cout << "XOff: i, j " << i << " " << j << " Geo size, Geo[i] size "  << fGeoInfo.size() << " ";
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << fGeoInfo[j].size() << endl;
  return fGeoInfo[j].at(i).xoffset;
}
//______________________________________________________________________
const double & TSolDBManager::GetDX(int i, int j)
{
  //cout << "DX: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size();
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << " " << fGeoInfo[j].size() << endl;
  return fGeoInfo[j].at(i).dx;
}
//______________________________________________________________________
const double & TSolDBManager::GetDY(int i, int j)
{
  // cout << "DY: i, j " << i << " " << j << " Geo size, Geo[i] size " << fGeoInfo.size();
  if (!CheckIndex(i, j)) return fErrVal;
  // cout << " " << fGeoInfo[j].size() << endl;
  return fGeoInfo[j].at(i).dy;
}
//______________________________________________________________________
const double & TSolDBManager::GetThetaH(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[j].at(i).thetaH;
}
//______________________________________________________________________
const double & TSolDBManager::GetThetaV(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[j].at(i).thetaV;
}
//_________________________________________________________________________
const double & TSolDBManager::GetStripAngle(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[j].at(i).stripangle_u;
    else return fGeoInfo[j].at(i).stripangle_u;
}
//_________________________________________________________________________
const double & TSolDBManager::GetPitch(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[j].at(i).pitch_u;
    else return fGeoInfo[j].at(i).pitch_u;
}

//__________________________________________________________________________
int TSolDBManager::GetSectorIDFromPos(int ichamber, double x, double y)
{
  if (!CheckIndex(ichamber)) return fErrVal;
  
  //int sector = -1;
  std::vector<int> sector;//yet something else...
  //printf("chamber %d, x = %1.6f\n", ichamber, x);
  for(int k = 0; k<fGeoInfo.size(); k++){
    //printf("%d: %1.2f, %1.2f \n", k, 
    //fGeoInfo[k].at(ichamber).xoffset-fGeoInfo[k].at(ichamber).dx/2.0,
    //fGeoInfo[k].at(ichamber).xoffset+fGeoInfo[k].at(ichamber).dx/2.0);
    if(fGeoInfo[k].at(ichamber).xoffset-fGeoInfo[k].at(ichamber).dx/2.0<x && 
       x<fGeoInfo[k].at(ichamber).xoffset+fGeoInfo[k].at(ichamber).dx/2.0){
      //if(abs(sector)==1)
      //sector = k;// EFuchey, 2017/05/16: if the sector value is -1 then it is not being over ridden
      // if the sector value is 1, then it is being over ridden, 
      // but it should happen only in the case of BB GEMs, plane 5.
      // I want that to happen because I divided plane 5 in 3 "sectors, each 1m long in x, 
      // and centered on -0.5, 0, 0.5, i.e. sector 1 overlaps completely with the 2 others.
      // I didn't figure any simpler way to integrate plane 5(UVA GEM) along with with the other 4 (INFN GEMs).
      sector.push_back(k);
    }
  }

  TRandom3 R(0);
  double P;
  //return sector;
  switch(sector.size()){
  case(0):
    return -1;
    break;
  case(1):
    return sector.at(0);
    break;
  case(2):// should actually happen **only** in the case of BB GEMs, plane 5.
    //return -1;cos(2*TMath::Pi()*x)/2+0.5
    P = R.Uniform(0, 1);
    if(P<= cos(2*TMath::Pi()*x)/2+0.5 ){
      return 1;
    }else{
      for(int i = 0; i<sector.size(); i++){
	if(sector.at(i)!=1)return sector.at(i);
      }
    }
    break;
  default:// >2 (should never happen)
    cout << "Sector size = " << sector.size() << "; something wrong, needs to be fixed !" << endl;
    return -1;
    break;
  }
  
}








