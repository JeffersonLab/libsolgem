#include "TSBSGeant4File.h"
#include "TSBSDBManager.h"
//#include "g4sbs_types.h"
#include "gemc_types.h"
#include "TRandom3.h"
#include "fstream"

using namespace std;

#ifndef __CINT__

// Set following variables to 1 (and recompile) t get some useful printouts
#ifndef D_FLAG
#define D_FLAG 0 //0: nothing; 1: warning; 2: debug;
#endif

TSBSGeant4File::TSBSGeant4File() : fFile(0), fTree(0), fSource(0), fEvNum(0), fManager(0) {
  fFilename[0] = '\0';
}

TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  //TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  SetFilename(f);
  fManager = TSBSDBManager::GetInstance();
  //InitMiscParam(filedbpath);
  
  //cout << " Gas data file container: -> " << &gasdatafilename << " <- " << endl;
  //cout << " Gas data file name: -> " << gasdatafilename.c_str() << " <- " << endl;
  //ReadGasData(gasdatafilename.c_str());
  //ReadGasData("gasErange.txt"); // NB: See comment lines 128-129 of TSBSGeant4File.h 
  
  //Filling the table that will be used to calculate the low energy electron range in the gas. 
#if D_FLAG>1
  cout << "Initialization completed" << endl;
#endif
}

/* // NB: See comment lines 138-141 of TSBSGeant4File.h
void TSBSGeant4File::ReadGasData(const char* filename){
  double D_gas;
  double T, p, R;
  
  ifstream in(filename);
  
  if(in.is_open()){
    in.ignore(100,';');
    in >> D_gas;
    in.ignore(50,';');
    
    p = -10.0;
    while(T<1000.0){//Kinetic energy cut to 0.1 MeV
      in >> T >> R;
      if(!in.good())break;
      
      p = T*sqrt(1.0+2.0*0.511/T)*1.0e-3;// in GeV
      feMom.push_back(p);
      fgasErange.push_back(R/D_gas*1.0e-2);// in m...
    }
  }else{
#if D_FLAG>1
    cout << "TSBSGeant4File Warning: file " << filename << " does not exist, using defaut values" << endl;
#endif
    D_gas = 1.662E-03;
    double eMom[81] = {
      1.000E-02, 1.250E-02, 1.500E-02, 1.750E-02, 2.000E-02, 2.500E-02, 3.000E-02, 3.500E-02, 4.000E-02, 
      4.500E-02, 5.000E-02, 5.500E-02, 6.000E-02, 7.000E-02, 8.000E-02, 9.000E-02, 1.000E-01, 1.250E-01, 
      1.500E-01, 1.750E-01, 2.000E-01, 2.500E-01, 3.000E-01, 3.500E-01, 4.000E-01, 4.500E-01, 5.000E-01, 
      5.500E-01, 6.000E-01, 7.000E-01, 8.000E-01, 9.000E-01, 1.000E+00, 1.250E+00, 1.500E+00, 1.750E+00, 
      2.000E+00, 2.500E+00, 3.000E+00, 3.500E+00, 4.000E+00, 4.500E+00, 5.000E+00, 5.500E+00, 6.000E+00, 
      7.000E+00, 8.000E+00, 9.000E+00, 1.000E+01, 1.250E+01, 1.500E+01, 1.750E+01, 2.000E+01, 2.500E+01, 
      3.000E+01, 3.500E+01, 4.000E+01, 4.500E+01, 5.000E+01, 5.500E+01, 6.000E+01, 7.000E+01, 8.000E+01, 
      9.000E+01, 1.000E+02, 1.250E+02, 1.500E+02, 1.750E+02, 2.000E+02, 2.500E+02, 3.000E+02, 3.500E+02, 
      4.000E+02, 4.500E+02, 5.000E+02, 5.500E+02, 6.000E+02, 7.000E+02, 8.000E+02, 9.000E+02, 1.000E+03
    };
    double gasErange[81] = {
      3.921E-04, 5.740E-04, 7.849E-04, 1.024E-03, 1.289E-03, 1.896E-03, 2.599E-03, 3.394E-03, 4.276E-03, 
      5.240E-03, 6.283E-03, 7.402E-03, 8.594E-03, 1.118E-02, 1.403E-02, 1.712E-02, 2.042E-02, 2.958E-02, 
      3.985E-02, 5.106E-02, 6.309E-02, 8.920E-02, 1.175E-01, 1.474E-01, 1.787E-01, 2.109E-01, 2.440E-01, 
      2.777E-01, 3.119E-01, 3.814E-01, 4.518E-01, 5.227E-01, 5.939E-01, 7.718E-01, 9.483E-01, 1.123E+00, 
      1.295E+00, 1.632E+00, 1.959E+00, 2.278E+00, 2.589E+00, 2.892E+00, 3.189E+00, 3.479E+00, 3.763E+00, 
      4.315E+00, 4.848E+00, 5.363E+00, 5.861E+00, 7.045E+00, 8.151E+00, 9.192E+00, 1.018E+01, 1.200E+01, 
      1.365E+01, 1.518E+01, 1.659E+01, 1.790E+01, 1.913E+01, 2.029E+01, 2.139E+01, 2.341E+01, 2.525E+01, 
      2.692E+01, 2.847E+01, 3.188E+01, 3.479E+01, 3.732E+01, 3.957E+01, 4.341E+01, 4.663E+01, 4.940E+01, 
      5.182E+01, 5.398E+01, 5.592E+01, 5.769E+01, 5.932E+01, 6.221E+01, 6.473E+01, 6.696E+01, 6.897E+01
    };
    for(int i = 0; i<81; i++){
      feMom.push_back(eMom[i]*sqrt(1.0+2.0*0.511/eMom[i])*1.0e-3);
      fgasErange.push_back(gasErange[i]/D_gas*1.0e-2);
    }
  }
}
*/ 

TSBSGeant4File::~TSBSGeant4File() {
  Clear();
  delete fFile;
}

void TSBSGeant4File::SetFilename( const char *f ){
  if( !f ) return;
  strcpy( fFilename, f );
}

Int_t TSBSGeant4File::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    delete fFile;
    fFile = new TFile(fFilename);
    
    if( !fFile->IsOpen() ){ 
      fprintf(stderr, "%s: File could not be made\n",__PRETTY_FUNCTION__);
      return 0; 
    }
    
    TChain* C1 = (TChain*)fFile->Get("T");//Get the tree from the file
    
#if D_FLAG>1 
    cout << "Detector option " << fManager->Getg4sbsDetectorType() << endl;
#endif
    
    // EFuchey: 2017/02/09: Since this date, the reading, digitization, etc. of the data 
    // has been splitted for Forward Tracker and Focal Plane Polarimeter.
    // This will make the reconstruction step easier to organize.
    // This implies an additional option for the detector type flag.
    // See the printout content below.
    if(fManager->Getg4sbsDetectorType()<1 && fManager->Getg4sbsDetectorType()>4){
      cout << "Invalid detector option: Set correct option in db_generalinfo.dat" << endl;
      cout << "(remider: 1 - BB GEMs; 2 - SIDIS SBS GEMs; 3 - FT; 4 - FPP)" << endl;
      return 0;
    }
    
    fTree = new g4sbs_tree(C1, fManager->Getg4sbsDetectorType());
    
    // g4sbs_tree declare all variables, branches, etc... 
    // to read, event by event, the varaibles stored in the tree. 
    // See comments in g4sbs_tree for more details...
    
    fEvNum = -1;
 
    return 1;
}

Int_t TSBSGeant4File::Close(){
    // Return 0 on fail, 1 on success
    Int_t ret = 1;
    
    if( !fFile->IsOpen() ){ return 0; }
    
    Clear();
    fFile->Close();
    
    delete fTree; 
    delete fFile; fFile = 0;
    return ret;
}



Int_t TSBSGeant4File::ReadNextEvent(int d_flag){
  // Return 1 on success
  
  // Channel not open
  if( !fFile->IsOpen() ){ 
    fprintf(stderr, "%s %s line %d Channel not open\n",
	    __FILE__,__PRETTY_FUNCTION__,__LINE__ );
    return 0; 
  }
  
  Clear();
  
  int n_hits = 0;//total number of hits at the end of the event
  int n_gen = 0;//total number of tracks at the end of the event
  // bool newtrk, dupli;// These variables help avoid store many times the same MC track info
  bool res = false;
  
  fEvNum++;

  //cout << "Read Next Event: Evt " << fEvNum << endl;
  
  res = fTree->GetEntry(fEvNum);
  //Test that the next entry exist
  if( !res ){
    // Don't need to print this out.  Not really an error
    if(d_flag>1){
      fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
	      __FILE__, __FUNCTION__, __LINE__ );
    } //DEBUG
    return 0;
  }
  
  double weight = fTree->ev_solang*fTree->ev_sigma; 
  
  //  int det_id;//2017/02/09: now corresponds to fManager->Getg4sbsDetectorType()
  
  //These are now local variables within each case statement below
  //as they logically should be
  // int pid;
  // int trid;
  // int type;
  // int plane;
  //FIXME: there seems to be confusion as to the meanings of "module" and "sector"
  // in the cases below. Keep them global for now, even though I strongly believe
  // they should also be local, like the other variables here
  int module = 0;
  int sector = 0;
  // double edep;
  // double tmin;
  // double tmax;
    
  // TVector3 Mom;
  // double pz;

  // TVector3 X_in;
  // TVector3 X_out;
  // TVector3 X_RO;
    
  // TVector3 Vtx;
  
  //FIXME: are these ever reset between modules?
  double hit_data_temp[24];
  double gen_data_temp[15];
  
  // NB: See comment lines 138-141 of TSBSGeant4File.h
  // double eRangeSlope;
  // double eRangeGas;
  // double temp;
  
  double pmax = 0.0;
  std::vector<int> trid_hits;
  double gen_data_temp_max[9];

  //variables for clustering
  TRandom3* R = new TRandom3(0);
  const double Npe_Edep_PS = 1.550e3;
  const double sigma_Npe_Edep_PS = 1.58e2;
  const double Npe_Edep_SH = 2.335e3;
  const double sigma_Npe_Edep_SH = 1.47e2;
  const int clustercrown_size = 2;
  const double BBECalBlock_size = 0.085;
  double edep_cal = 0;
  double npe = 0;
  
  double Edep_PS = 0;
  double Edep_SH = 0;
  
  double Edep_PS_max = 0;
  double Edep_SH_max = 0;
  double EdepXmax = 0;
  double EdepYmax = 0;
  
  double E_rec = 0;
  double X_rec = 0;
  double Y_rec = 0;
  double E_rec_SH = 0;
  
  std::vector<double> bbps_edep;
  //std::vector<double> bbps_xcell;
  std::vector<double> bbps_ycell;
  std::vector<double> bbsh_edep;
  std::vector<double> bbsh_xcell;
  std::vector<double> bbsh_ycell;
  
  switch(fManager->Getg4sbsDetectorType()){

  case(1)://BB GEMs
    if(d_flag>1){
      cout << "Number of Hits: " << fTree->Earm_BBGEM_hit_nhits << endl;
    } // DEBUG
    for(int i = 0; i<fTree->Earm_BBGEM_hit_nhits; i++){
      int det_id = 1;
      int pid = fTree->Earm_BBGEM_hit_pid->at(i);
      int trid = fTree->Earm_BBGEM_hit_trid->at(i);// track ID: particle counter
      int type = fTree->Earm_BBGEM_hit_mid->at(i)+1;//=1 if primary, >1 if secondary...
      int plane = fTree->Earm_BBGEM_hit_plane->at(i)-1;
      double edep = fTree->Earm_BBGEM_hit_edep->at(i)*1.0e3;
      double tmin = fTree->Earm_BBGEM_hit_tmin->at(i);
      double tmax = fTree->Earm_BBGEM_hit_tmax->at(i);
      
      trid_hits.push_back(trid);
      
      module = fManager->GetModuleIDFromPos(plane, fTree->Earm_BBGEM_hit_tx->at(i));
      if(module==-1)continue;
      
      double pz = sqrt( pow(fTree->Earm_BBGEM_hit_p->at(i), 2)/
		 ( pow(fTree->Earm_BBGEM_hit_txp->at(i), 2) + 
		   pow(fTree->Earm_BBGEM_hit_typ->at(i), 2) + 1.0) );
     
      TVector3 Mom = TVector3(fTree->Earm_BBGEM_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Earm_BBGEM_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      TVector3 X_in = TVector3((fTree->Earm_BBGEM_hit_xin->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Earm_BBGEM_hit_yin->at(i)*1.0e3, // in mm
		      (fTree->Earm_BBGEM_hit_zin->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
      
      TVector3 X_out = TVector3((fTree->Earm_BBGEM_hit_xout->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		       fTree->Earm_BBGEM_hit_yout->at(i)*1.0e3, // in mm
		       (fTree->Earm_BBGEM_hit_zout->at(i)+fManager->Getg4sbsZSpecOffset()-
			fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(fTree->Earm_BBGEM_hit_xout->at(i)-fTree->Earm_BBGEM_hit_xin->at(i))*9.185/3.0 ,
		      // in mm
		      X_in.Y()+(fTree->Earm_BBGEM_hit_yout->at(i)-fTree->Earm_BBGEM_hit_yin->at(i))*9.185/3.0 ,
		      // in mm
		      +7.685);// in mm
      /* // see comment lines 138-141 of TSBSGeant4File.h
      X_in = TVector3((fTree->Earm_BBGEM_hit_tx->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Earm_BBGEM_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Earm_BBGEM_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
       
      //cout << fTree->Earm_BBGEM_hit_z->at(i) << "+" << fManager->Getg4sbsZSpecOffset() << "-" << fManager->GetD0(plane, sector) << " = " << X_in.z() << endl;
      
      X_out = TVector3(X_in.X()+3.0*fTree->Earm_BBGEM_hit_txp->at(i), // in mm 
		       X_in.Y()+3.0*fTree->Earm_BBGEM_hit_typ->at(i), // in mm
		       +1.5);//-1.6825);// in mm
      
      //if(X_out.X()-X_in.X()>0.01)
      
      X_RO = TVector3(X_in.X()+9.185*fTree->Earm_BBGEM_hit_txp->at(i), // in mm 
		      X_in.Y()+9.185*fTree->Earm_BBGEM_hit_typ->at(i), // in mm 
		      +7.685);//4.5025);// in mm      
      
      //cout << "BBGEM: momentum: " << fTree->Earm_BBGEM_hit_p->at(i) << " < ? " << feMom.back() << endl;
      if(fabs(fTree->Earm_BBGEM_hit_pid->at(i))==11 && fTree->Earm_BBGEM_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Earm_BBGEM_hit_txp->at(i), 2)+pow(fTree->Earm_BBGEM_hit_typ->at(i), 2))*3.0e-3;//m
	if(eRangeSlope>=5.0e-3 && fTree->Earm_BBGEM_hit_p->at(i)>=1.0e-2)cout << "e Range Slope " << eRangeSlope << " p (MeV) = " << fTree->Earm_BBGEM_hit_p->at(i)*1000 << " hit edep = " << edep << endl;
	eRangeGas = FindGasRange(fTree->Earm_BBGEM_hit_p->at(i));//m
	if(d_flag>1)cout << " range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(X_in.X()+3.0*fTree->Earm_BBGEM_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(X_in.Y()+3.0*fTree->Earm_BBGEM_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(X_out.X());
	  X_RO.SetY(X_out.Y());
       	}
      }
      */
      
      if(fabs(X_out.X())>=fManager->GetDX(plane, module)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.X " << X_out.X() << " outside BBGEM plane " << plane
	       << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, module)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << fTree->Earm_BBGEM_hit_tx->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }  
      if(fabs(X_out.Y())>=fManager->GetDY(plane, module)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.Y " << X_out.Y() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, module)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << fTree->Earm_BBGEM_hit_ty->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }
      
      TVector3 Vtx = TVector3(fTree->Earm_BBGEM_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Earm_BBGEM_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Earm_BBGEM_hit_vz->at(i)*1.0e3);// in mm
      
      //Filling hit_data temporary array...
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      hit_data_temp[19] = module;
      hit_data_temp[23] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  24));//data_size(__GEM_TAG)));

      // ... to copy it in the actual g4sbsHitData structure.
      for(int j = 0; j<24; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
      

            
      if(d_flag>1){
	cout << "Hit number: " << i << " BBGEM: X_global : " 
	     << fTree->Earm_BBGEM_hit_xg->at(i) << ", " 
	     << fTree->Earm_BBGEM_hit_yg->at(i) << ", " 
	     << fTree->Earm_BBGEM_hit_zg->at(i) << endl;
	cout << "X_local (g4sbs): " << fTree->Earm_BBGEM_hit_tx->at(i) << ", " 
	     << fTree->Earm_BBGEM_hit_ty->at(i) << ", " 
	     << fTree->Earm_BBGEM_hit_z->at(i) << endl;
	cout << "detector ID: " << det_id << ", plane: " << plane << ", sector: " << sector << endl
	     << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	     << "energy deposit (eV): " << edep << endl;
	cout << "Momentum (MeV): ";
	for(int k = 0; k<3; k++){
	  cout << Mom[k] << ", ";
	}
	cout << " norm " << fTree->Earm_BBGEM_hit_p->at(i) << endl;
	cout << "hit position at drift entrance (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_in[k] << ", ";
	}
	cout << " time : " << tmin << endl;
	cout << "hit position at drift exit (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_out[k] << " ";
	}
	cout << " time : " << tmax << endl;
	cout << "hit position at readout (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_RO[k] << ", ";
	}
	cout << endl;
	cout << "Vertex position (mm): ";
	for(int k = 0; k<3; k++){
	  cout << Vtx[k] << ", ";
	}
	cout << endl;
      } //DEBUG      
    }//end loop on hits
    
    
    //Filling gen_data temporary array...
    
    
    // ... to copy it in the actual g4sbsGenData structure.
    // only store signal, primary MC tracks
    if(fSource==0 && n_gen==0 && fTree->Earm_BBGEM_Track_ntracks>0  && fTree->Earm_BBGEM_Track_MID->at(0)==0){
      
      TVector3 trackOrigin = TVector3(fTree->Earm_BBGEM_Track_X->at(0)*1.0e3, // in mm
				      fTree->Earm_BBGEM_Track_Y->at(0)*1.0e3, // in mm
				      0);
      double pz = fTree->Earm_BBGEM_Track_P->at(0) / sqrt(1 + pow(fTree->Earm_BBGEM_Track_Xp->at(0),2) + pow(fTree->Earm_BBGEM_Track_Yp->at(0),2));
      TVector3 Mom = TVector3(pz * fTree->Earm_BBGEM_Track_Xp->at(0),
			      pz * fTree->Earm_BBGEM_Track_Yp->at(0),
			      pz);
      TVector3 trackVertex = TVector3(fTree->ev_vx, fTree->ev_vy, fTree->ev_vz);
      TVector3 trackVertexMom = TVector3(fTree->ev_epx, fTree->ev_epy, fTree->ev_epz);
      
      
      gen_data_temp[0] = fTree->Earm_BBGEM_Track_TID->at(0);
      gen_data_temp[1] = fTree->Earm_BBGEM_Track_PID->at(0);
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];  //average momentum in tracker                in GeV
	gen_data_temp[k+5] = trackOrigin[k];//interception at focal plane          in mm
	gen_data_temp[k+9] = trackVertexMom[k];//momentum at target in global axis    in GeV
	gen_data_temp[k+12]= trackVertex[k];//vertex in global axis,                  in m
      }
      gen_data_temp[8] = weight;
      
      fg4sbsGenData.push_back(new g4sbsgendata());
      
      int jdd=0;
      for(; jdd<15; jdd++){
	fg4sbsGenData[n_gen]->SetData(jdd, gen_data_temp[jdd]);
      }
      n_gen++;
    }
    
    
      
    // ---------------------------------
    // Add calorimeter clustering here.
    // ---------------------------------
    bbps_edep.clear();
    bbps_ycell.clear();
    bbsh_edep.clear();
    bbsh_xcell.clear();
    bbsh_ycell.clear();
    //first, loop on PS hits
    for(Int_t i = 0; i<fTree->Earm_BBPSTF1_hit_nhits; i++){
      // evaluate energy deposit including smearing from the photoelectron yield:
      // convert energy deposit in pe yield with Npe_Edep_PS conversion coefficiency, 
      // and smearing sigma_Npe_Edep_PS. 
      // Then reconvert into enegry deposit using Npe_Edep_PS
      edep_cal = fTree->Earm_BBPSTF1_hit_sumedep->at(i);
      npe = R->Gaus(Npe_Edep_PS*edep_cal, sigma_Npe_Edep_PS*edep_cal);
      //Npe_fEdep_PS(R, fTree->Earm_BBPSTF1_hit_sumedep->at(i));
      edep_cal = npe/Npe_Edep_PS;
      //Edep_fNpe_PS(npe);
      Edep_PS+= edep_cal;
	
      // if edep is the maximal energy deposit, store it in the "max" variables
      if(edep_cal>Edep_PS_max){
	Edep_PS_max = edep_cal;
	EdepYmax = fTree->Earm_BBPSTF1_hit_ycell->at(i);
      }
	
      // record all smeared energy deposits no matter what, + photon coordinates 
      bbps_edep.push_back(edep_cal);
      // X coordinates (in calorimeter) are not relevant for PS.
      //bbps_xcell.push_back(fTree->Earm_BBPSTF1_hit_xcell->at(i));
      bbps_ycell.push_back(fTree->Earm_BBPSTF1_hit_ycell->at(i));
    }

    //then, loop on SH hits
    // same story as for PS, except we also store X coordinates
    for(Int_t i = 0; i<fTree->Earm_BBSHTF1_hit_nhits; i++){
      edep_cal = fTree->Earm_BBSHTF1_hit_sumedep->at(i);
      npe = R->Gaus(Npe_Edep_SH*edep_cal, sigma_Npe_Edep_SH*edep_cal);
      //Npe_fEdep_SH(R, fTree->Earm_BBSHTF1_hit_sumedep->at(i));
      edep_cal = npe/Npe_Edep_SH;
      //Edep_fNpe_SH(npe);
      Edep_SH+= edep_cal;
	
      if(edep_cal>Edep_SH_max){
	Edep_SH_max = edep_cal;
	EdepXmax = fTree->Earm_BBSHTF1_hit_xcell->at(i);
	if(Edep_SH_max>Edep_PS_max){
	  EdepYmax = fTree->Earm_BBSHTF1_hit_ycell->at(i);
	}
      }
	
      bbsh_edep.push_back(edep_cal);
      bbsh_xcell.push_back(fTree->Earm_BBSHTF1_hit_xcell->at(i));
      bbsh_ycell.push_back(fTree->Earm_BBSHTF1_hit_ycell->at(i));
    }
      
    // calculate reconstructed energy: 5x5 blocks around the max.
    for(size_t l = 0; l<bbps_edep.size(); l++){
      if(fabs(bbps_ycell[l]-EdepYmax)<=BBECalBlock_size*clustercrown_size+1.0e-2){
	E_rec+= bbps_edep[l];
	//X_rec+= bbps_edep[l]*(-1)*bbps_ycell[l];
      }
    }
      
    // calculate reconstructed energy: 5x5 blocks around the max.
    for(size_t l = 0; l<bbsh_edep.size(); l++){
      if(fabs(bbsh_xcell[l]-EdepXmax)<=BBECalBlock_size*clustercrown_size+1.0e-2 &&
	 fabs(bbsh_ycell[l]-EdepYmax)<=BBECalBlock_size*clustercrown_size+1.0e-2){
	E_rec+= bbsh_edep[l];
	E_rec_SH+= bbsh_edep[l];
	//transforming calo coordinates as stored in g4sbs output to transport coordinates.
	X_rec+= bbsh_edep[l]*(-1)*bbsh_ycell[l];
	Y_rec+= bbsh_edep[l]*bbsh_xcell[l];
      }
    }
      
    // calculate reconstructed position: 
    // mean of position of all cluster blocks weighted with energy deposit.
    //X_rec = X_rec/E_rec;
    X_rec = X_rec/E_rec_SH;
    Y_rec = Y_rec/E_rec_SH;
    //TSBSECalCluster clus(E_rec, X_rec, Y_rec);

    //hardcode threshold for the moment, will add in the DB later.
    if(E_rec>fManager->GetCaloThreshold()){
      fECalClusters.push_back(new TSBSECalCluster(E_rec, X_rec, Y_rec));
    }
    // -------------------------------
    // end: calorimeter reconstruction
    // -------------------------------
      
    /*
    // Rescue ngen data here...
    // ngen data is rescued if: 
    // * the file read is signal; 
    // * there is no gen data stored already; 
    // * and there are at least 3 hits associated to the rescued particle
    if(fSource==0 && n_gen==0 && n_hits>=3){
    int n_hits_max = 0;
    for(UInt_t k = 0; k<trid_hits.size(); k++){
    if(gen_data_temp_max[0]==trid_hits[k])n_hits_max++;
    }
    if(n_hits_max>=3){
    fg4sbsGenData.push_back(new g4sbsgendata());
    for(int j = 0; j<9; j++){
    fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
    }
    n_gen++;
    }
    }
    */
    break;
  case(2)://SIDIS SBS GEMs
    for(int i = 0; i<fTree->Harm_SBSGEM_hit_nhits; i++){
      int det_id = 2;

      int pid = fTree->Harm_SBSGEM_hit_pid->at(i);
      int trid = fTree->Harm_SBSGEM_hit_trid->at(i);
      int type = fTree->Harm_SBSGEM_hit_mid->at(i)+1;//=1 if primary, >1 if secondary...
      int plane = fTree->Harm_SBSGEM_hit_plane->at(i)-1;
      int edep = fTree->Harm_SBSGEM_hit_edep->at(i)*1.0e3;
      int tmin = fTree->Harm_SBSGEM_hit_tmin->at(i);
      int tmax = fTree->Harm_SBSGEM_hit_tmax->at(i);
      //FIXME: note how 'sector' is set here and 'module' isn't, but then 'module'
      // is used in the Get... calls right below, holding the value left over from
      // a different case
      sector = fManager->GetModuleIDFromPos(plane, fTree->Harm_SBSGEM_hit_tx->at(i));

      trid_hits.push_back(trid);

      double pz = sqrt( pow(fTree->Harm_SBSGEM_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_SBSGEM_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_SBSGEM_hit_typ->at(i), 2) + 1.0) );
      
      TVector3 Mom = TVector3(fTree->Harm_SBSGEM_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_SBSGEM_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      TVector3 X_in = TVector3((fTree->Harm_SBSGEM_hit_xin->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Harm_SBSGEM_hit_yin->at(i)*1.0e3, // in mm
		      (fTree->Harm_SBSGEM_hit_zin->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
      
      TVector3 X_out = TVector3((fTree->Harm_SBSGEM_hit_xout->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		       fTree->Harm_SBSGEM_hit_yout->at(i)*1.0e3, // in mm
		       (fTree->Harm_SBSGEM_hit_zout->at(i)+fManager->Getg4sbsZSpecOffset()-
			fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(fTree->Harm_SBSGEM_hit_xout->at(i)-fTree->Harm_SBSGEM_hit_xin->at(i))*9.185/3.0 ,
		      // in mm
		      X_in.Y()+(fTree->Harm_SBSGEM_hit_yout->at(i)-fTree->Harm_SBSGEM_hit_yin->at(i))*9.185/3.0 ,
		      // in mm
		      +7.685);// in mm
      
      /* // see comment lines 138-141 of TSBSGeant4File.h
      X_in = TVector3((fTree->Harm_SBSGEM_hit_tx->at(i)-fManager->GetXOffset(plane, sector))*1.0e3, // in mm
		      fTree->Harm_SBSGEM_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_SBSGEM_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, sector))*1.0e3);// in mm
      
      X_out = TVector3(X_in.X()+3.0*fTree->Harm_SBSGEM_hit_txp->at(i), // in mm 
		       X_in.Y()+3.0*fTree->Harm_SBSGEM_hit_typ->at(i), // in mm
		       +1.5);//-1.6825);// in mm

      X_RO = TVector3(X_in.X()+9.185*fTree->Harm_SBSGEM_hit_txp->at(i), // in mm 
		      X_in.Y()+9.185*fTree->Harm_SBSGEM_hit_typ->at(i), // in mm 
		      +7.685);//4.5025);// in mm      
      
      if(fabs(fTree->Harm_SBSGEM_hit_pid->at(i))==11 && fTree->Harm_SBSGEM_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_SBSGEM_hit_txp->at(i), 2)+pow(fTree->Harm_SBSGEM_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_SBSGEM_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
       	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(X_in.X()+3.0*fTree->Harm_SBSGEM_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(X_in.Y()+3.0*fTree->Harm_SBSGEM_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(X_out.X());
	  X_RO.SetY(X_out.Y());
       	}
      }
      */
      if(fabs(X_out.X())>=fManager->GetDX(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.X " << X_out.X() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << fTree->Harm_SBSGEM_hit_tx->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=fManager->GetDY(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.Y " << X_out.Y() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << fTree->Harm_SBSGEM_hit_ty->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }
      
      TVector3 Vtx = TVector3(fTree->Harm_SBSGEM_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_SBSGEM_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_SBSGEM_hit_vz->at(i)*1.0e3);// in mm

      //Filling hit_data temporary array...
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      hit_data_temp[19] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  data_size(__GEM_TAG)));

      // ... to copy it in the actual g4sbsHitData structure.
      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
      
      //Filling gen_data temporary array...
      gen_data_temp[0] = trid;
      gen_data_temp[1] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];
	gen_data_temp[k+5] = Vtx[k];
      }
      gen_data_temp[8] = weight;
      
      //store information to rescue, if necessary, the generated info
      if(fTree->Harm_SBSGEM_hit_p->at(i)>pmax && pid==fManager->GetSigPID(0)){
	pmax = fTree->Harm_SBSGEM_hit_p->at(i);
	for(int k = 0; k<9; k++){
	  gen_data_temp_max[k] = gen_data_temp[k];
	}
      }
      
      // ... to copy it in the actual g4sbsGenData structure.
      // only store signal, primary MC tracks
      if(fSource==0 && n_gen==0 && type==1){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      /*else{// this determines if the track is new or not
	newtrk = true; 
	for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(1)!=gen_data_temp[1]){
	dupli=false;
	}else{
	for(int j = 5; j<8; j++){
	if(fg4sbsGenData[z]->GetData(j)!=gen_data_temp[j]){
	dupli=false;
	break;
	}
	}
	}
	if(dupli){
	newtrk = false;
	break;
	}
	}
	
	if(newtrk){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
	}
	}
      */
    }// endl loop on hits
    
    // Rescue ngen data here...
    // ngen data is rescued if: 
    // * the file read is signal; 
    // * there is no gen data stored already; 
    // * and there are at least 3 hits associated to the rescued particle
    if(fSource==0 && n_gen==0 && n_hits>=3){
      int n_hits_max = 0;
      for(UInt_t k = 0; k<trid_hits.size(); k++){
	if(gen_data_temp_max[0]==trid_hits[k])n_hits_max++;
      }
      if(n_hits_max>=3){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
    }
    break;
      
  case(3)://FT
    //Loop on the Forward Tracker detector hits:
    for(int i = 0; i<fTree->Harm_FT_hit_nhits; i++){
      int det_id = 3;
      int pid = fTree->Harm_FT_hit_pid->at(i);
      int trid = fTree->Harm_FT_hit_trid->at(i);
      int type = fTree->Harm_FT_hit_mid->at(i)+1;//=1 if primary, >1 if secondary...
      int plane = fTree->Harm_FT_hit_plane->at(i)-1;
      double edep = fTree->Harm_FT_hit_edep->at(i)*1.0e3;
      double tmin = fTree->Harm_FT_hit_tmin->at(i);
      double tmax = fTree->Harm_FT_hit_tmax->at(i);
      //cout << "pid ??" << pid << endl;
      module = fManager->GetModuleIDFromPos(plane, fTree->Harm_FT_hit_tx->at(i));
      //cout << fTree->Harm_FT_hit_xin->at(i) << " " << fTree->Harm_FT_hit_xout->at(i) << endl;
      
      trid_hits.push_back(trid);

      double pz = sqrt( pow(fTree->Harm_FT_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FT_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FT_hit_typ->at(i), 2) + 1.0) );
      
      TVector3 Mom = TVector3(fTree->Harm_FT_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FT_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      TVector3 X_in = TVector3((fTree->Harm_FT_hit_xin->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Harm_FT_hit_yin->at(i)*1.0e3, // in mm
		      (fTree->Harm_FT_hit_zin->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
      
      TVector3 X_out = TVector3((fTree->Harm_FT_hit_xout->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		       fTree->Harm_FT_hit_yout->at(i)*1.0e3, // in mm
		       (fTree->Harm_FT_hit_zout->at(i)+fManager->Getg4sbsZSpecOffset()-
			fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(fTree->Harm_FT_hit_xout->at(i)-fTree->Harm_FT_hit_xin->at(i))*9.185/3.0 ,
		      // in mm
		      X_in.Y()+(fTree->Harm_FT_hit_yout->at(i)-fTree->Harm_FT_hit_yin->at(i))*9.185/3.0 ,
		      // in mm
		      +7.685);// in mm

      /* // see comment lines 138-141 of TSBSGeant4File.h
      X_in = TVector3((fTree->Harm_FT_hit_tx->at(i)-fManager->GetXOffset(plane, sector))*1.0e3, // in mm
		      fTree->Harm_FT_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_FT_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, sector))*1.0e3);// in mm
      
      X_out = TVector3(X_in.X()+3.0*fTree->Harm_FT_hit_txp->at(i), // in mm 
		       X_in.Y()+3.0*fTree->Harm_FT_hit_typ->at(i), // in mm
		       +1.5);//-1.6825);// in mm
      
      X_RO = TVector3(X_in.X()+9.185*fTree->Harm_FT_hit_txp->at(i), // in mm 
		      X_in.Y()+9.185*fTree->Harm_FT_hit_typ->at(i), // in mm 
		      +7.685);//4.5025);// in mm      
      
      //Calculation of very low momentum electrons range ingas.
      if(fabs(fTree->Harm_FT_hit_pid->at(i))==11 && fTree->Harm_FT_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FT_hit_txp->at(i), 2)+pow(fTree->Harm_FT_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FT_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(X_in.X()+3.0*fTree->Harm_FT_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(X_in.Y()+3.0*fTree->Harm_FT_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(X_out.X());
	  X_RO.SetY(X_out.Y());
	  //cout << "Coucou ! FT" << endl;
       	}
      }
      */
      if(fabs(X_out.X())>=fManager->GetDX(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.X " << X_out.X() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << fTree->Harm_FT_hit_tx->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=fManager->GetDY(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << i 
	       << ": X_out.Y " << X_out.Y() << " outside FT plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << fTree->Harm_FT_hit_ty->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }
      
      TVector3 Vtx = TVector3(fTree->Harm_FT_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vz->at(i)*1.0e3);// in mm

      //Filling hit_data temporary array...
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      //cout << "pid ??" << pid << " hit_data_temp[18] ?? " << hit_data_temp[18] << endl;
      hit_data_temp[19] = module;
      //cout << "module !!!!!!!!!!! " << module << " " << hit_data_temp[19] << endl;
      hit_data_temp[23] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  24));

      // ... to copy it in the actual g4sbsHitData structure.
      for(int j = 0; j<24; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
 
      /*
      //Filling gen_data temporary array...
      gen_data_temp[0] = trid;
      gen_data_temp[1] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];
	gen_data_temp[k+5] = Vtx[k];
      }
      gen_data_temp[8] = weight;
      
      //store information to rescue, if necessary, the generated info
      if(fTree->Harm_FT_hit_p->at(i)>pmax && pid==fManager->GetSigPID(0)){
	pmax = fTree->Harm_FT_hit_p->at(i);
	for(int k = 0; k<9; k++){
	  gen_data_temp_max[k] = gen_data_temp[k];
	}
      }
      
      // ... to copy it in the actual g4sbsGenData structure.
      // only store signal, primary MC tracks
      if(fSource==0 && n_gen==0 && type==1){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      */
      /*else{// this determines if the track is new or not
	newtrk = true; 
	for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(1)!=gen_data_temp[1]){
	dupli=false;
	}else{
	for(int j = 5; j<8; j++){
	if(fg4sbsGenData[z]->GetData(j)!=gen_data_temp[j]){
	dupli=false;
	break;
	}
	}
	}
	if(dupli){
	newtrk = false;
	break;
	}
	}
	
	if(newtrk){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
	}
	}
      */
      // Print out block
      if(d_flag>1){
	cout << "Hit number: " << i << " FT: X_global : " 
	     << fTree->Harm_FT_hit_xg->at(i) << ", " 
	     << fTree->Harm_FT_hit_yg->at(i) << ", " 
	     << fTree->Harm_FT_hit_zg->at(i) << endl;
	cout << "detector ID: " << det_id << ", plane: " << plane << endl
	     << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	     << "energy deposit (eV): " << edep << endl;
	cout << "Momentum (MeV): ";
	for(int k = 0; k<3; k++){
	  cout << Mom[k] << ", ";
	}
	cout << " norm " << fTree->Harm_FT_hit_p->at(i) << endl
	     << "dpx/dpz = " << fTree->Harm_FT_hit_txp->at(i)
	     << ", dpx/dpz = " << fTree->Harm_FT_hit_typ->at(i)
	     << endl;
	cout << "hit position at drift entrance (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_in[k] << ", ";
	}
	cout << " time : " << tmin << endl;
	cout << "hit position at drift exit (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_out[k] << " ";
	}
	cout << " time : " << tmax << endl;
	cout << "hit position at readout (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_RO[k] << ", ";
	}
	cout << endl;
	cout << "Vertex position (mm): ";
	for(int k = 0; k<3; k++){
	  cout << Vtx[k] << ", ";
	}
	cout << endl;
      } //DEBUG
    }
    
    // only store signal, primary MC tracks
    if(fSource==0 && n_gen==0 && fTree->Harm_FT_Track_ntracks>0  && fTree->Harm_FT_Track_MID->at(0)==0){
      
      TVector3 trackOrigin = TVector3(fTree->Harm_FT_Track_X->at(0)*1.0e3, // in mm
				      fTree->Harm_FT_Track_Y->at(0)*1.0e3, // in mm
				      0);
      double pz = fTree->Harm_FT_Track_P->at(0) / sqrt(1 + pow(fTree->Harm_FT_Track_Xp->at(0),2) + pow(fTree->Harm_FT_Track_Yp->at(0),2));
      TVector3 Mom = TVector3(pz * fTree->Harm_FT_Track_Xp->at(0),
			      pz * fTree->Harm_FT_Track_Yp->at(0),
			      pz);
      TVector3 trackVertex = TVector3(fTree->ev_vx, fTree->ev_vy, fTree->ev_vz);
      TVector3 trackVertexMom = TVector3(fTree->ev_epx, fTree->ev_epy, fTree->ev_epz);
      
      
      gen_data_temp[0] = fTree->Harm_FT_Track_TID->at(0);
      gen_data_temp[1] = fTree->Harm_FT_Track_PID->at(0);
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];  //average momentum in tracker                in GeV
	gen_data_temp[k+5] = trackOrigin[k];//interception at focal plane          in mm
	gen_data_temp[k+9] = trackVertexMom[k];//momentum at target in global axis    in GeV
	gen_data_temp[k+12]= trackVertex[k];//vertex in global axis,                  in m
      }
      gen_data_temp[8] = weight;
      
      fg4sbsGenData.push_back(new g4sbsgendata());
      
      int jdd=0;
      for(; jdd<15; jdd++){
	fg4sbsGenData[n_gen]->SetData(jdd, gen_data_temp[jdd]);
      }
      n_gen++;
    }
        
    /*
    // Rescue ngen data here...
    // ngen data is rescued if: 
    // * the file read is signal; 
    // * there is no gen data stored already; 
    // * and there are at least 3 hits associated to the rescued particle
    if(fSource==0 && n_gen==0 && n_hits>=3){
      int n_hits_max = 0;
      for(UInt_t k = 0; k<trid_hits.size(); k++){
	if(gen_data_temp_max[0]==trid_hits[k])n_hits_max++;
      }
      if(n_hits_max>=3){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      
    }
    */
    break;// end case(3)
    
  case(4):
    //Loop on the Focal Plane Polarimeter 1 hits: detectors 0 to 4
    // This block is not well commented, 
    // as it is very similar to the previous block of instructions
    // where Forward Tracker data are unfolded.
    for(int i = 0; i<fTree->Harm_FPP1_hit_nhits; i++){
      int det_id = 4;

      int pid = fTree->Harm_FPP1_hit_pid->at(i);
      int trid = fTree->Harm_FPP1_hit_trid->at(i);
      int type = fTree->Harm_FPP1_hit_mid->at(i)+1;//=1 if primary, >1 if secondary...
      int plane = fTree->Harm_FPP1_hit_plane->at(i)-1;
      double edep = fTree->Harm_FPP1_hit_edep->at(i)*1.0e3;
      double tmin = fTree->Harm_FPP1_hit_tmin->at(i);
      double tmax = fTree->Harm_FPP1_hit_tmax->at(i);
      
      if(d_flag>1)cout << plane << " " << fTree->Harm_FPP1_hit_tx->at(i) << endl; 
      
      module = fManager->GetModuleIDFromPos(plane, fTree->Harm_FPP1_hit_tx->at(i));

      if(d_flag>1)cout << plane << ",  " << sector << ",  " << fTree->Harm_FPP1_hit_tx->at(i) << endl; 
      
      trid_hits.push_back(trid);
      
      double pz = sqrt( pow(fTree->Harm_FPP1_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP1_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP1_hit_typ->at(i), 2) + 1.0) );
      
      TVector3 Mom = TVector3(fTree->Harm_FPP1_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP1_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      TVector3 X_in = TVector3((fTree->Harm_FPP1_hit_xin->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Harm_FPP1_hit_yin->at(i)*1.0e3, // in mm
		      (fTree->Harm_FPP1_hit_zin->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
      
      TVector3 X_out = TVector3((fTree->Harm_FPP1_hit_xout->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		       fTree->Harm_FPP1_hit_yout->at(i)*1.0e3, // in mm
		       (fTree->Harm_FPP1_hit_zout->at(i)+fManager->Getg4sbsZSpecOffset()-
			fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(fTree->Harm_FPP1_hit_xout->at(i)-fTree->Harm_FPP1_hit_xin->at(i))*9.185/3.0 ,
		      // in mm
		      X_in.Y()+(fTree->Harm_FPP1_hit_yout->at(i)-fTree->Harm_FPP1_hit_yin->at(i))*9.185/3.0 ,
		      // in mm
		      +7.685);// in mm

      /* // see comment lines 138-141 of TSBSGeant4File.h
      X_in = TVector3((fTree->Harm_FPP1_hit_tx->at(i)-fManager->GetXOffset(plane, sector))*1.0e3, // in mm
		      fTree->Harm_FPP1_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_FPP1_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, sector))*1.0e3);// in mm
      
      X_out = TVector3(X_in.X()+3.0*fTree->Harm_FPP1_hit_txp->at(i), // in mm 
		       X_in.Y()+3.0*fTree->Harm_FPP1_hit_typ->at(i), // in mm
		       -1.6825);// in mm
      
      X_RO = TVector3(X_in.X()+9.185*fTree->Harm_FPP1_hit_txp->at(i), // in mm 
		      X_in.Y()+9.185*fTree->Harm_FPP1_hit_typ->at(i), // in mm 
		      4.5025);// in mm      
       
      if(d_flag>1)cout << "FPP1: momentum: " << fTree->Harm_FPP1_hit_p->at(i) << " < ? " << feMom.back() << endl;
      if(fabs(fTree->Harm_FPP1_hit_pid->at(i))==11 && fTree->Harm_FPP1_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FPP1_hit_txp->at(i), 2)+pow(fTree->Harm_FPP1_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FPP1_hit_p->at(i));//m
	if(d_flag>1)cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
       	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(X_in.X()+3.0*fTree->Harm_FPP1_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(X_in.Y()+3.0*fTree->Harm_FPP1_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(X_out.X());
	  X_RO.SetY(X_out.Y());
	  if(d_flag>1)cout << "Coucou ! FPP1 " << endl;
       	}
      }
      */
      
      if(fabs(X_out.X())>=fManager->GetDX(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << fTree->Harm_FT_hit_nhits+i 
	       << ": X_out.X " << X_out.X() << " outside FPP1 plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << fTree->Harm_FPP1_hit_tx->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=fManager->GetDY(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " << fTree->Harm_FT_hit_nhits+i 
	       << ": X_out.Y " << X_out.Y() << " outside FPP1 plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << fTree->Harm_FPP1_hit_ty->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }
      
      TVector3 Vtx = TVector3(fTree->Harm_FPP1_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vz->at(i)*1.0e3);// in mm
      
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      hit_data_temp[19] = module;
      hit_data_temp[23] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id, 24));

      for(int j = 0; j<24; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
      
      /*
      //Filling gen_data temporary array...
      gen_data_temp[0] = trid;
      gen_data_temp[1] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];
	gen_data_temp[k+5] = Vtx[k];
      }
      gen_data_temp[8] = weight;
      
      //store information to rescue, if necessary, the generated info
      if(fTree->Harm_FPP1_hit_p->at(i)>pmax && pid==fManager->GetSigPID(0)){
	pmax = fTree->Harm_FPP1_hit_p->at(i);
	for(int k = 0; k<9; k++){
	  gen_data_temp_max[k] = gen_data_temp[k];
	}
      }
      
      // ... to copy it in the actual g4sbsGenData structure.
      // only store signal, primary MC tracks
      if(fSource==0 && n_gen==0 && type==1){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      */
      
      /*else{// this determines if the track is new or not
	newtrk = true; 
	for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(1)!=gen_data_temp[1]){
	dupli=false;
	}else{
	for(int j = 5; j<8; j++){
	if(fg4sbsGenData[z]->GetData(j)!=gen_data_temp[j]){
	dupli=false;
	break;
	}
	}
	}
	if(dupli){
	newtrk = false;
	break;
	}
	}
	
	if(newtrk){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
	}
	}
      */
      if(d_flag>1){
	cout << "Hit number: " << fTree->Harm_FT_hit_nhits+i << " FPP1: X_global : " 
	     << fTree->Harm_FPP1_hit_xg->at(i) << ", " 
	     << fTree->Harm_FPP1_hit_yg->at(i) << ", " 
	     << fTree->Harm_FPP1_hit_zg->at(i) << endl;
	cout << "detector ID: " << det_id << ", plane: " << plane << endl
	     << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	     << "energy deposit (MeV): " << edep << endl;
	cout << "Momentum (MeV): ";
	for(int k = 0; k<3; k++){
	  cout << Mom[k] << ", ";
	}
	cout << " norm " << fTree->Harm_FPP1_hit_p->at(i) << endl;
	cout << "hit position at drift entrance (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_in[k] << ", ";
	}
	cout << " time : " << tmin << endl;
	cout << "hit position at drift exit (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_out[k] << " ";
	}
	cout << " time : " << tmax << endl;
	cout << "hit position at readout (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_RO[k] << ", ";
	}
	cout << endl;
	cout << "Vertex position (mm): ";
	for(int k = 0; k<3; k++){
	  cout << Vtx[k] << ", ";
	}
	cout << endl;
      } //DEBUG 
      
      if(fSource==0 && n_gen==0 && fTree->Harm_FPP1_Track_ntracks>0  && fTree->Harm_FPP1_Track_PID->at(0)==2212){
	
	TVector3 trackOrigin = TVector3(fTree->Harm_FPP1_Track_X->at(0)*1.0e3, // in mm
					fTree->Harm_FPP1_Track_Y->at(0)*1.0e3, // in mm
					0);
	double pz = fTree->Harm_FPP1_Track_P->at(0) / sqrt(1 + pow(fTree->Harm_FPP1_Track_Xp->at(0),2) + pow(fTree->Harm_FPP1_Track_Yp->at(0),2));
	TVector3 Mom = TVector3(pz * fTree->Harm_FPP1_Track_Xp->at(0),
				pz * fTree->Harm_FPP1_Track_Yp->at(0),
				pz);
	TVector3 trackVertex = TVector3(fTree->ev_vx, fTree->ev_vy, fTree->ev_vz);
	TVector3 trackVertexMom = TVector3(fTree->ev_epx, fTree->ev_epy, fTree->ev_epz);
	
	
	gen_data_temp[0] = fTree->Harm_FPP1_Track_TID->at(0);
	gen_data_temp[1] = fTree->Harm_FPP1_Track_PID->at(0);
	for(int k = 0; k<3; k++){
	  gen_data_temp[k+2] = Mom[k];  //average momentum in tracker                in GeV
	  gen_data_temp[k+5] = trackOrigin[k];//interception at focal plane          in mm
	  gen_data_temp[k+9] = trackVertexMom[k];//momentum at target in global axis    in GeV
	  gen_data_temp[k+12]= trackVertex[k];//vertex in global axis,                  in m
	}
	gen_data_temp[8] = weight;
	
	fg4sbsGenData.push_back(new g4sbsgendata());
	
	int jdd=0;
	for(; jdd<15; jdd++){
	  fg4sbsGenData[n_gen]->SetData(jdd, gen_data_temp[jdd]);
	}
	n_gen++;
      }
      
    }
    
    //Loop on the Focal Plane Polarimeter 2 hits: detectors 5 to 9
    // This block is not well commented, 
    // as it is very similar to the previous block of instructions
    // where Forward Tracker data are unfolded.
    if(d_flag>0)cout << "read FPP2" << endl;
    for(int i = 0; i<fTree->Harm_FPP2_hit_nhits; i++){
      int det_id = 4;

      int pid = fTree->Harm_FPP2_hit_pid->at(i);
      int trid = fTree->Harm_FPP2_hit_trid->at(i);
      int type = fTree->Harm_FPP2_hit_mid->at(i)+1;//=1 if primary, >1 if secondary...
      int plane = fManager->GetNChamber()/2+fTree->Harm_FPP2_hit_plane->at(i)-1;
      double edep = fTree->Harm_FPP2_hit_edep->at(i)*1.0e3;
      double tmin = fTree->Harm_FPP2_hit_tmin->at(i);
      double tmax = fTree->Harm_FPP2_hit_tmax->at(i);
      module = fManager->GetModuleIDFromPos(plane, fTree->Harm_FPP2_hit_tx->at(i));
      
      trid_hits.push_back(trid);

      double pz = sqrt( pow(fTree->Harm_FPP2_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP2_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP2_hit_typ->at(i), 2) + 1.0) );
      
      TVector3 Mom = TVector3(fTree->Harm_FPP2_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP2_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      TVector3 X_in = TVector3((fTree->Harm_FPP2_hit_xin->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_yin->at(i)*1.0e3, // in mm
		      (fTree->Harm_FPP2_hit_zin->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, module))*1.0e3);// in mm
      
      TVector3 X_out = TVector3((fTree->Harm_FPP2_hit_xout->at(i)-fManager->GetXOffset(plane, module))*1.0e3, // in mm
		       fTree->Harm_FPP2_hit_yout->at(i)*1.0e3, // in mm
		       (fTree->Harm_FPP2_hit_zout->at(i)+fManager->Getg4sbsZSpecOffset()-
			fManager->GetD0(plane, module))*1.0e3);// in mm

      // we use X_in and X_out to extrapolate X_RO; 
      // not very clean, but since we don't use it in the digitization at all, it does not really matter...
      TVector3 X_RO = TVector3(X_in.X()+(fTree->Harm_FPP2_hit_xout->at(i)-fTree->Harm_FPP2_hit_xin->at(i))*9.185/3.0 ,
		      // in mm
		      X_in.Y()+(fTree->Harm_FPP2_hit_yout->at(i)-fTree->Harm_FPP2_hit_yin->at(i))*9.185/3.0 ,
		      // in mm
		      +7.685);// in mm

      /* // see comment lines 138-141 of TSBSGeant4File.h
      X_in = TVector3((fTree->Harm_FPP2_hit_tx->at(i)-fManager->GetXOffset(plane, sector))*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_FPP2_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()-
		       fManager->GetD0(plane, sector))*1.0e3);// in mm
      
      X_out = TVector3(X_in.X()+3.0*fTree->Harm_FPP2_hit_txp->at(i), // in mm 
		       X_in.Y()+3.0*fTree->Harm_FPP2_hit_typ->at(i), // in mm
		       -1.6825);// in mm
      
      X_RO = TVector3(X_in.X()+9.185*fTree->Harm_FPP2_hit_txp->at(i), // in mm 
		      X_in.Y()+9.185*fTree->Harm_FPP2_hit_typ->at(i), // in mm 
		      4.5025);// in mm      
      
      if(fabs(fTree->Harm_FPP2_hit_pid->at(i))==11 && fTree->Harm_FPP2_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FPP2_hit_txp->at(i), 2)+pow(fTree->Harm_FPP2_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FPP2_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
       	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(X_in.X()+3.0*fTree->Harm_FPP2_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(X_in.Y()+3.0*fTree->Harm_FPP2_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(X_out.X());
	  X_RO.SetY(X_out.Y());
       	}
      }
      */
      
      if(fabs(X_out.X())>=fManager->GetDX(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " 
	       << fTree->Harm_FPP1_hit_nhits+fTree->Harm_FT_hit_nhits+i 
	       << ": X_out.X " << X_out.X() << " outside FPP2 plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.X());
	X_out[0]*=fManager->GetDX(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.X() << " mm " << endl;
	  cout << "(X_in.X = " << X_in.X() << ",  " << fTree->Harm_FPP2_hit_tx->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=fManager->GetDY(plane, sector)*5.0e2){
	if(d_flag>0){
	  cout << "Warning: Evt " << fEvNum << ", hit " 
	       << fTree->Harm_FPP1_hit_nhits+fTree->Harm_FT_hit_nhits+i 
	       << ": X_out.Y " << X_out.Y() << " outside FPP2 plane " << plane << " sector " << sector;
	} //D_FLAG
	double temp = fabs(X_out.Y());
	X_out[1]*=fManager->GetDY(plane, sector)*5.0e2/temp;
	if(d_flag>0){
	  cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
	  cout << "(X_in.Y = " << X_in.Y() << ",  " << fTree->Harm_FPP2_hit_ty->at(i)*1.0e3 << " mm)" << endl;
	} //D_FLAG
	X_RO.SetY(X_out.Y());
      }

      TVector3 Vtx = TVector3(fTree->Harm_FPP2_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vz->at(i)*1.0e3);// in mm

      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = (double)type;
      hit_data_temp[17] = (double)trid;
      hit_data_temp[18] = (double)pid;
      hit_data_temp[19] = module;
      hit_data_temp[23] = sector;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id, 24));
      
      for(int j = 0; j<24; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;

      /*
      //Filling gen_data temporary array...
      gen_data_temp[0] = trid;
      gen_data_temp[1] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+2] = Mom[k];
	gen_data_temp[k+5] = Vtx[k];
      }
      gen_data_temp[8] = weight;
      
      //store information to rescue, if necessary, the generated info
      if(fTree->Harm_FPP2_hit_p->at(i)>pmax && pid==fManager->GetSigPID(0)){
	pmax = fTree->Harm_FPP2_hit_p->at(i);
	for(int k = 0; k<9; k++){
	  gen_data_temp_max[k] = gen_data_temp[k];
	}
      }
      
      // ... to copy it in the actual g4sbsGenData structure.
      // only store signal, primary MC tracks
      if(fSource==0 && n_gen==0 && type==1){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
	}
      */
      /*else{// this determines if the track is new or not
	newtrk = true; 
	for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(1)!=gen_data_temp[1]){
	dupli=false;
	}else{
	for(int j = 5; j<8; j++){
	if(fg4sbsGenData[z]->GetData(j)!=gen_data_temp[j]){
	dupli=false;
	break;
	}
	}
	}
	if(dupli){
	newtrk = false;
	break;
	}
	}
	
	if(newtrk){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
	}
	}
      */
      if(d_flag>1){
	cout << "Hit number: " << fTree->Harm_FPP1_hit_nhits+fTree->Harm_FT_hit_nhits+i << " FPP2: X_global : " 
	     << fTree->Harm_FPP2_hit_xg->at(i) << ", " 
	     << fTree->Harm_FPP2_hit_yg->at(i) << ", " 
	     << fTree->Harm_FPP2_hit_zg->at(i) << endl;
	cout << "detector ID: " << det_id << ", plane: " << plane << endl
	     << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	     << "energy deposit (eV): " << edep << endl;
	cout << "Momentum (MeV): ";
	for(int k = 0; k<3; k++){
	  cout << Mom[k] << ", ";
	}
	cout << " norm " << fTree->Harm_FPP2_hit_p->at(i) << endl;
	cout << "hit position at drift entrance (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_in[k] << ", ";
	}
	cout << " time : " << tmin << endl;
	cout << "hit position at drift exit (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_out[k] << " ";
	}
	cout << " time : " << tmax << endl;
	cout << "hit position at readout (mm): ";
	for(int k = 0; k<3; k++){
	  cout << X_RO[k] << ", ";
	}
	cout << endl;
	cout << "Vertex position (mm): ";
	for(int k = 0; k<3; k++){
	  cout << Vtx[k] << ", ";
	}
	cout << endl;
      } //DEBUG          
    }
    
    /*
    // Rescue ngen data here...
    // ngen data is rescued if: 
    // * the file read is signal; 
    // * there is no gen data stored already; 
    // * and there are at least 3 hits associated to the rescued particle
    if(fSource==0 && n_gen==0 && n_hits>=3){
      int n_hits_max = 0;
      for(UInt_t k = 0; k<trid_hits.size(); k++){
	if(gen_data_temp_max[0]==trid_hits[k])n_hits_max++;
      }
      if(n_hits_max>=3){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<9; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
    }
    */
    break;// end case(4)
      
   }//end switch(fManager->Getg4sbsDetectorType)
    
  delete R;
  return 1;
}


void TSBSGeant4File::Clear(){
  // Clear out hit and generated data

#if D_FLAG>1
  fprintf(stderr, "%s %s line %d: Deleting hits\n",
	  __FILE__, __FUNCTION__, __LINE__);
#endif //DEBUG

  unsigned int i;
  for( i = 0; i < fg4sbsHitData.size(); i++ ){
    delete fg4sbsHitData[i];
  }

  for( i = 0; i < fg4sbsGenData.size(); i++ ){
    delete fg4sbsGenData[i];
  }

  fg4sbsHitData.clear();
  fg4sbsGenData.clear();
  
  fECalClusters.clear();
  
#if D_FLAG>1
  fprintf(stderr, "%s %s line %d: Hits deleted\n",
	  __FILE__, __FUNCTION__, __LINE__);
#endif //DEBUG

  return;
}

/* // NB: See comment lines 138-141 of TSBSGeant4File.h
double TSBSGeant4File::FindGasRange(double p)
{
  //find the electron range in the gas. Useful for very low energy electrons. 
  //momentum cut set rather high (~1.1 MeV) should be more than enough
  double res = -10.0;

  if(p<feMom.at(0)){
    //cout << "Momentum " << p << " < " << feMom.at(0) << endl;
    res = fgasErange.at(0)*exp(log(p/feMom.at(0))*log(fgasErange.at(1)/fgasErange.at(0))/log(feMom.at(1)/feMom.at(0)));
    return (res);
  }
  
  for(uint i = 1; i<feMom.size(); i++){
    if(feMom.at(i-1)<= p && p<feMom.at(i)){
      res = fgasErange.at(i-1)*exp(log(p/feMom.at(i-1))*log(fgasErange.at(i)/fgasErange.at(i-1))/log(feMom.at(i)/feMom.at(i-1)));
      //cout << "Momentum: " << feMom.at(i-1) << " < " << p << " < " << feMom.at(i) << endl;
      //cout << "Range: " << fgasErange.at(i-1) << " <? " << res << " <? " << fgasErange.at(i) << endl;
      return(res);
    }
  }
  return(res);
}
*/

TSBSGEMSimHitData* TSBSGeant4File::GetGEMData()
{
  // Return TSBSGEMSimHitData object filled with GEM data of present event.
  // The returned object pointer must be deleted by the caller!

  TSBSGEMSimHitData* gd = new TSBSGEMSimHitData();

  GetGEMData(gd);
  return gd;
}

void TSBSGeant4File::GetGEMData(TSBSGEMSimHitData* gd)
{
  // Pack data into TSBSGEMSimHitData
   
  //    printf("NEXT EVENT ---------------------------\n");

  if( !gd ) return;
  gd->ClearEvent();
  gd->SetSource(fSource);
  gd->SetEvent(fEvNum);
  
  if (GetNData() == 0) {
    return;
  }
  gd->InitEvent(GetNData());

  g4sbshitdata *h;//, *hs;
  //bool matchedstrip;
  unsigned int i, ngdata = 0;// j,
  for(i=0; i<GetNData(); i++){
    h = GetHitData(i);

    if( h->GetData(1)>0.0 ){
	
      // Chamber IDs are numbered as 
      	
      //if( h->GetDetID()%100 == __GEM_DRIFT_ID &&  h->GetData(1)>0.0 ){
      // Vector information
      TVector3 p(h->GetData(20), h->GetData(21), h->GetData(22));
      gd->SetMomentum(ngdata, p);
	
      TVector3 li(h->GetData(5), h->GetData(6), h->GetData(7));
      gd->SetHitEntrance(ngdata, li);
	
      TVector3 lo(h->GetData(9), h->GetData(10), h->GetData(11));
      gd->SetHitExit(ngdata, lo);
	
      // Average over entrance and exit time
      gd->SetHitTime(ngdata, (h->GetData(8)+h->GetData(12))/2.0);
	
      TVector3 vert(h->GetData(14), h->GetData(15), h->GetData(16));
      gd->SetVertex(ngdata, vert);
	
      TVector3 lr(h->GetData(2), h->GetData(3), h->GetData(4));
      gd->SetHitReadout(ngdata, lr);
      // printf("%d %f %f\n", h->GetDetID()/100, li.X(), li.Y()  );
	
      gd->SetHitEnergy(ngdata, h->GetData(1)*1e6 ); // Gives eV
      gd->SetParticleType(ngdata, (UInt_t)h->GetData(13) );//  Track type (1 primary, >1 secondary) 
      gd->SetTrackID(ngdata, (UInt_t) h->GetData(17) );// track ID
      gd->SetParticleID(ngdata, h->GetData(18) );//  PID 
      //cout << "ngdata = " << ngdata << " TSBSGeant4File::GetGEMData: pid " <<h->GetData(18)<<" : "<<gd->GetParticleID(ngdata)<<endl;
      
      // gd->SetHitChamber(ngdata, h->GetData(23)*fManager->GetNChamber()+h->GetData(0));
      gd->SetHitChamber(ngdata, fManager->GetGEMID(h->GetData(0), h->GetData(19)));
      //cout<<(h->GetData(23)*fManager->GetNChamber()+h->GetData(0))<<" : "<<(h->GetData(0)*3+h->GetData(19))<<endl;
      
      gd->SetHitPlane(ngdata, h->GetData(0));
      gd->SetHitModule(ngdata,h->GetData(19));
      
      //cout << "ngdata = " << ngdata << "TSBSGeant4File::GetGEMData: igem " << gd->GetHitChamber(ngdata) << ", plane " << gd->GetHitPlane(ngdata) << ", module" <<  gd->GetHitModule(ngdata) << endl;
      
      ngdata++;
    }
  }
  gd->SetNHit(ngdata);
}


///////////////////////////////////////////////////////////////
// hitdata classes 

g4sbshitdata::g4sbshitdata(int detid, unsigned int size ){
  fDetID = detid;
  fData  = new double[size];
  fSize  = size;
  fFillbits = 0;
  
  if( size > sizeof( long long int )*8 ){
    fprintf(stderr, "%s %s line %d:  Error:  Event size too long for bit pattern storage (requested %d, have %ld)\n",
	    __FILE__, __PRETTY_FUNCTION__, __LINE__, size, 
	    sizeof(long long int)*8);
    exit(1);
  }
  
  // There is no value indexed at 0, so we'll just set it to 0 for
  // sanity's sake and not deal with crazy offsets all over
  
  fFillbits |= 1;
  fData[0] = 3.1415927;
}

void g4sbshitdata::SetData(unsigned int idx, double data ){
  if( idx >= fSize ){
    fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
    return;

  }

  fFillbits |= (1<<idx);

  fData[idx] = data;
  return;
}

double g4sbshitdata::GetData(unsigned int idx) const {
  if( idx >= fSize ){
    fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
    return 1e9;
  }

  if( !(fFillbits & (1<<idx)) ){
    fprintf(stderr, "%s %s line %d:  Error:  Accessing unset data (idx %d) val: %f\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fData[idx] );
    return 1e9;
  }

  return fData[idx];
}

bool g4sbshitdata::IsFilled() const {
  if( fFillbits == ((1<<fSize) - 1) ){
	return true;
    }

    return false;
}

g4sbshitdata::~g4sbshitdata(){
    delete fData;
}

///////////////////////////////////////////////////////////////
// gendata classes

// Size is 1 bigger because we are also including the weight
// Set that default to 1
g4sbsgendata::g4sbsgendata():g4sbshitdata(-1, __GENERATED_SIZE+2){
    SetData(8,1.0);
}

#endif//__CINT__
