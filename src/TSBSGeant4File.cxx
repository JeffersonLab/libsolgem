#include "TSBSGeant4File.h"
//#include "g4sbs_types.h"
#include "gemc_types.h"
#include "fstream"

#ifndef __CINT__

// Set following variables to 1 (and recompile) t get some useful printouts
#define DEBUG 0
#define WARNING 1

TSBSGeant4File::TSBSGeant4File() : fFile(0), fSource(0) {
  fFilename[0] = '\0';
}

TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  //TSBSGeant4File::TSBSGeant4File(const char *f) : fFile(0), fSource(0) {
  SetFilename(f);
  fManager = TSolDBManager::GetInstance();
  //InitMiscParam(filedbpath);
  
  //string gasdatafilename = fManager->GetGasDataFilename();
  //cout << " Gas data file container: -> " << &gasdatafilename << " <- " << endl;
  //cout << " Gas data file name: -> " << gasdatafilename.c_str() << " <- " << endl;
  //ReadGasData(gasdatafilename.c_str());
  ReadGasData("gasErange.txt");
  
  //Filling the table that will be used to calculate the low energy electron range in the gas. 
  
  cout << "Initialization completed" << endl;
}

void TSBSGeant4File::ReadGasData(const char* filename){
  double D_gas;
  double T, p, R;
  
  ifstream in(filename);
  
  if(in.is_open()){
    in.ignore(100,';');
    in >> D_gas;
    in.ignore(50,';');
    
    p = -10.0;
    while(T<0.1){//Kinetic energy cut to 0.1 MeV
      in >> T >> R;
      if(!in.good())break;
      
      p = T*sqrt(1.0+2.0*0.511/T)*1.0e-3;// in GeV
      feMom.push_back(p);
      fgasErange.push_back(R/D_gas*1.0e-2);// in m...
    }
  }else{
#if WARNING>1
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
      feMom.push_back(eMom[i]);
      fgasErange.push_back(gasErange[i]);
    }
  }
}

/*
//-----------------------------------------------------------------------------
// Reading database for miscellaneous parameters.
// This is done without proper DB request, but there is a set of default parameters.
// If those have to be used, user will be warned by a warning message.
// Data should be sorted as in file db/db_g4sbsmiscdata.dat
// This is the user's responsibility to make sure his input file is read correctly.
void TSBSGeant4File::InitMiscParam(const char* dbpath) {
  ifstream in(dbpath);
  if(!in.is_open()){
    printf("TSBSGeant4File Warning: May not read database at %s\n", dbpath);
    printf(" => Using sbs default params\n");
    
    fZSpecOffset = 3.38551;
    strcpy( fgasdatafile, "gasErange.txt");
  }else{
    cout << "TSBSGeant4File Info: reading database at location " << dbpath << endl;
    cout <<" This file should be written the same way db/db_g4sbsmiscdata.dat "<< endl;
    cout << "(same structure, same order of parameters)" << endl;
    Float_t dummy;
    //string read_str;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;

    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    in.ignore(100,'=');
    in >> dummy;
    
    in.ignore(100,'=');
    in >> fZSpecOffset;
    in.ignore(100,'=');
    in >> fgasdatafile;
    cout << "Gas range data input file: " << fgasdatafile << endl;
  }
  
  
  // cout << fNSECTORS << " " << fNPLANES << " " << fNPROJ << " " << fCHAN_PER_SLOT << " "
  // 	 << fmodules_per_readout << "; " << endl;
  //   << fgZ0 << " " << fgDoCalo << " " << fgCaloZ << " " << fgCaloRes << endl;
  
  in.close();
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

    fTree = new g4sbs_gep_tree_with_spin(C1, false);
    // g4sbs_gep_tree_with_spin declare all variables, branches, etc... 
    // to read, event by event, the varaibles stored in the tree. 
    // See comments in g4sbs_gep_tree_with_spin for more details...

    fEvNum = -1;
 
    return 1;
}

Int_t TSBSGeant4File::Close(){
    // Return 0 on fail, 1 on success
    Int_t ret = 1;
    
    if( !fFile->IsOpen() ){ return 0; }
    
    fFile->Close();
    
    delete fFile; fFile = 0;
    return ret;
}

Int_t TSBSGeant4File::ReadNextEvent(){
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
    bool newtrk, dupli;// These variables help avoid store many times the same MC track info
    bool res = false;
    
    fEvNum++;
    
    res = fTree->GetEntry(fEvNum);
    //Test that the next entry exist
    if( !res ){
      // Don't need to print this out.  Not really an error
#if DEBUG>0
      fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
	      __FILE__, __FUNCTION__, __LINE__ );
#endif //DEBUG
      return 0;
    }
    
    double weight = fTree->ev_solang*fTree->ev_sigma; 
    
    int det_id;//0: FT, 1: FPPs
    
    int pid;
    int type;
    int plane;
    double edep;
    double tmin;
    double tmax;
    
    TVector3 Mom;
    double pz;

    TVector3 X_in;
    TVector3 X_out;
    TVector3 X_RO;
    
    TVector3 Vtx;
    
    double hit_data_temp[23];
    double gen_data_temp[8];
    
    //variables for the correction of hits given by very small momenta
    double eRangeSlope;
    double eRangeGas;
    double temp;
    
    //if(fManager->)
    //Loop on the Forward Tracker detector hits: detectors 10 to 15
    for(int i = 0; i<fTree->Harm_FT_hit_nhits; i++){
      det_id = 1;
      
      pid = fTree->Harm_FT_hit_pid->at(i);
      type = fTree->Harm_FT_hit_mid->at(i)+1;
      plane = fTree->Harm_FT_hit_plane->at(i);
      edep = fTree->Harm_FT_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FT_hit_tmin->at(i);
      tmax = fTree->Harm_FT_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FT_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FT_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FT_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FT_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FT_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FT_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_FT_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_txp->at(i), // in mm 
		       fTree->Harm_FT_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_typ->at(i), // in mm
		       (fTree->Harm_FT_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3+9.185*fTree->Harm_FT_hit_txp->at(i), // in mm 
		      fTree->Harm_FT_hit_ty->at(i)*1.0e3+9.185*fTree->Harm_FT_hit_typ->at(i), // in mm 
		      (fTree->Harm_FT_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+9.185);// in mm
      
      //cout << "FT: momentum: " << fTree->Harm_FT_hit_p->at(i) << " < ? " << feMom.back() << endl;
      
      //Calculation of very low momentum electrons range ingas.
      if(fabs(fTree->Harm_FT_hit_pid->at(i))==11 && fTree->Harm_FT_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FT_hit_txp->at(i), 2)+pow(fTree->Harm_FT_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FT_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(fTree->Harm_FT_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(fTree->Harm_FT_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_typ->at(i)*eRangeGas/eRangeSlope);

	  X_RO.SetX(fTree->Harm_FT_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_RO.SetY(fTree->Harm_FT_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  //cout << "Coucou ! FT" << endl;
       	}
      }
      
      //Correcting X_out x and y if out of the GEM plane...
      if(fabs(X_out.X())>=749.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " << i 
	     << ": X_out.X " << X_out.X() << " outside FT plane " << 10+plane;
#endif //WARNING
	temp = fabs(X_out.X());
	X_out[0]*=749.99/temp;
#if WARNING>0
	cout  << "; set at limit: " << X_out.X() << " mm " << endl;
#endif //WARNING
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=199.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " << i 
	     << ": X_out.Y " << X_out.Y() << " outside FT plane " << 10+plane;
#endif //WARNING
	temp = fabs(X_out.Y());
	X_out[1]*=199.99/temp;
#if WARNING>0
	cout  << "; set at limit: " << X_out.Y() << " mm " << endl;	
#endif //WARNING
	X_RO.SetY(X_out.Y());
      }
      
      Vtx = TVector3(fTree->Harm_FT_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vz->at(i)*1.0e3);// in mm

      //Filling hit_data temporary array...
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
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
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;
      
      // ... to copy it in the actual g4sbsGenData structure.
      // only store new MC tracks
      if(n_gen==0){
	fg4sbsGenData.push_back(new g4sbsgendata());
	for(int j = 0; j<8; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }else{// this determines if the track is new or not
	newtrk = true; 
	for(int z = n_gen-1; z>=0; z--){
	  dupli = true;
	  if(fg4sbsGenData[z]->GetData(0)!=gen_data_temp[0]){
	    dupli=false;
	  }else{
	    for(int j = 4; j<8; j++){
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
	  for(int j = 0; j<8; j++){
	    fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	  }
	  n_gen++;
	}
      }
      
      // Print out block
#if DEBUG>0
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
#endif //DEBUG
    }
    
    //Loop on the Focal Plane Polarimeter 1 hits: detectors 0 to 4
    // This block is not well commented, 
    // as it is very similar to the previous block of instructions
    // where Forward Tracker data are unfolded.
    for(int i = 0; i<fTree->Harm_FPP1_hit_nhits; i++){
      det_id = 0;
      
      pid = fTree->Harm_FPP1_hit_pid->at(i);
      type = fTree->Harm_FPP1_hit_mid->at(i)+1;
      plane = fTree->Harm_FPP1_hit_plane->at(i);
      edep = fTree->Harm_FPP1_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FPP1_hit_tmin->at(i);
      tmax = fTree->Harm_FPP1_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FPP1_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP1_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP1_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FPP1_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP1_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP1_hit_ty->at(i)*1.0e3, // in mm
		      (fTree->Harm_FPP1_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_txp->at(i), 
		       fTree->Harm_FPP1_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_typ->at(i), 
		       (fTree->Harm_FPP1_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3+9.185*fTree->Harm_FPP1_hit_txp->at(i), 
		      fTree->Harm_FPP1_hit_ty->at(i)*1.0e3+9.185*fTree->Harm_FPP1_hit_typ->at(i), 
		      (fTree->Harm_FPP1_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+9.185);// in mm
      
      //cout << "FPP1: momentum: " << fTree->Harm_FPP1_hit_p->at(i) << " < ? " << feMom.back() << endl;
      if(fabs(fTree->Harm_FPP1_hit_pid->at(i))==11 && fTree->Harm_FPP1_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FPP1_hit_txp->at(i), 2)+pow(fTree->Harm_FPP1_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FPP1_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
       	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(fTree->Harm_FPP1_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_RO.SetY(fTree->Harm_FPP1_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  //cout << "Coucou ! FPP1 " << endl;
       	}
      }
         
      if(fabs(X_out.X())>=999.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " << fTree->Harm_FT_hit_nhits+i 
	     << ": X_out.X " << X_out.X() << " outside FPP1 plane " << plane;
#endif //WARNING
	temp = fabs(X_out.X());
	X_out[0]*=999.99/temp;
#if WARNING>0
	cout << "; set at limit: " << X_out.X() << " mm " << endl;
#endif //WARNING
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=299.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " << fTree->Harm_FT_hit_nhits+i 
	     << ": X_out.Y " << X_out.Y() << " outside FPP1 plane " << plane;
#endif //WARNING
	temp = fabs(X_out.Y());
	X_out[1]*=299.99/temp;
#if WARNING>0
	cout << "; set at limit: " << X_out.Y() << " mm " << endl;
#endif //WARNING
	X_RO.SetY(X_out.Y());
      }

      Vtx = TVector3(fTree->Harm_FPP1_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vz->at(i)*1.0e3);// in mm
      
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  data_size(__GEM_TAG)));

      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
      
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;
      
      newtrk = true; 
      for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(0)!=gen_data_temp[0]){
	  dupli=false;
	}else{
	  for(int j = 4; j<8; j++){
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
	for(int j = 0; j<8; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      
#if DEBUG>0
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

      cout << "detector ID: " << det_id << ", plane: " << plane << endl
	   << "particle ID: " << pid << ", type (1, primary, >1 secondary: " << type << endl
	   << "energy deposit (GeV): " << edep << endl;
      cout << "Momentum (GeV): ";
      Mom.Print();
      cout << endl;
      cout << "hit position at drift entrance (mm): ";
      X_in.Print();
      cout << " time : " << tmin << endl;
      cout << "hit position at drift exit (mm): ";
      X_out.Print();
      cout << " time : " << tmax << endl;
      cout << "hit position at Readout (mm): ";
      X_RO.Print();
      cout << endl;
      cout << "Vertex position (mm): ";
      Vtx.Print();
      cout << endl;
#endif //DEBUG      
    }
    
    //Loop on the Focal Plane Polarimeter 2 hits: detectors 5 to 9
    // This block is not well commented, 
    // as it is very similar to the previous block of instructions
    // where Forward Tracker data are unfolded.
    for(int i = 0; i<fTree->Harm_FPP2_hit_nhits; i++){
      det_id = 0;
      
      pid = fTree->Harm_FPP2_hit_pid->at(i);
      type = fTree->Harm_FPP2_hit_mid->at(i)+1;
      plane = 5+fTree->Harm_FPP2_hit_plane->at(i);
      edep = fTree->Harm_FPP2_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FPP2_hit_tmin->at(i);
      tmax = fTree->Harm_FPP2_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FPP2_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP2_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP2_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FPP2_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP2_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_ty->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_z->at(i)+fManager->Getg4sbsZSpecOffset()*1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_txp->at(i), // in mm 
		       fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_typ->at(i), // in mm
		       (fTree->Harm_FPP2_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+9.185*fTree->Harm_FPP2_hit_txp->at(i), // in mm
		      fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+9.185*fTree->Harm_FPP2_hit_typ->at(i), // in mm
		      (fTree->Harm_FPP2_hit_z->at(i)+fManager->Getg4sbsZSpecOffset())*1.0e3+9.185);// in mm
      
      //cout << "FPP2: momentum: " << fTree->Harm_FPP2_hit_p->at(i) << " < ? " << feMom.back() << endl;
      if(fabs(fTree->Harm_FPP2_hit_pid->at(i))==11 && fTree->Harm_FPP2_hit_p->at(i)<=feMom.back()){
	eRangeSlope = sqrt(pow(fTree->Harm_FPP2_hit_txp->at(i), 2)+pow(fTree->Harm_FPP2_hit_typ->at(i), 2))*3.0e-3;//m
	eRangeGas = FindGasRange(fTree->Harm_FPP2_hit_p->at(i));//m
	//cout << "range: " << eRangeGas << " < ? "  << eRangeSlope << endl;
       	if(eRangeSlope>eRangeGas){
       	  X_out.SetX(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_out.SetY(fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_typ->at(i)*eRangeGas/eRangeSlope);
	  
	  X_RO.SetX(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_txp->at(i)*eRangeGas/eRangeSlope);
	  X_RO.SetY(fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_typ->at(i)*eRangeGas/eRangeSlope);
       	}
      }
      
      if(fabs(X_out.X())>=999.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " 
	     << fTree->Harm_FPP1_hit_nhits+fTree->Harm_FT_hit_nhits+i 
	     << ": X_out.X " << X_out.X() << " outside FPP2 plane " << plane;
#endif //WARNING
	temp = fabs(X_out.X());
	X_out[0]*=999.99/temp;
#if WARNING>0
	cout  << "; set at limit: " << X_out.X() << " mm " << endl;
#endif //WARNING
	X_RO.SetX(X_out.X());
      }
      if(fabs(X_out.Y())>=299.99){
#if WARNING>0
	cout << "Warning: Evt " << fEvNum << ", hit " 
	     << fTree->Harm_FPP1_hit_nhits+fTree->Harm_FT_hit_nhits+i 
	     << ": X_out.Y " << X_out.Y() << " outside FPP2 plane " << plane;
#endif //WARNING
	temp = fabs(X_out.Y());
	X_out[1]*=299.99/temp;
#if WARNING>0
	cout  << "; set at limit: " << X_out.Y() << " mm " << endl;
#endif //WARNING
	X_RO.SetY(X_out.Y());
      }
      
      Vtx = TVector3(fTree->Harm_FPP2_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vz->at(i)*1.0e3);// in mm

      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id, data_size(__GEM_TAG)));
      
      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      n_hits++;
      
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;

      newtrk = true; 
      for(int z = n_gen-1; z>=0; z--){
	dupli = true;
	if(fg4sbsGenData[z]->GetData(0)!=gen_data_temp[0]){
	  dupli=false;
	}else{
	  for(int j = 4; j<8; j++){
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
	for(int j = 0; j<8; j++){
	  fg4sbsGenData[n_gen]->SetData(j, gen_data_temp[j]);
	}
	n_gen++;
      }
      
#if DEBUG>0
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
#endif //DEBUG          
    }
    
    return 1;
}


void TSBSGeant4File::Clear(){
    // Clear out hit and generated data

#if DEBUG>0
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

#if DEBUG>0
	fprintf(stderr, "%s %s line %d: Hits deleted\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif //DEBUG

    return;
}

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


TSolGEMData* TSBSGeant4File::GetGEMData()
{
  // Return TSolGEMData object filled with GEM data of present event.
  // The returned object pointer must be deleted by the caller!

  TSolGEMData* gd = new TSolGEMData();

  GetGEMData(gd);
  return gd;
}

void TSBSGeant4File::GetGEMData(TSolGEMData* gd)
{
  // Pack data into TSolGEMData
   
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
	// xy  where x is the det id (1 for FT, 0 for FPPs) y the plane number, 
	// labeled from 0 to 9 instead of 1 to 10, for convenience reasons:
	// FPPs: chambers 0-9, FT: chambers 10-15.
	
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
	gd->SetTrackID(ngdata, (UInt_t) h->GetData(13) );// track type...
	gd->SetParticleType(ngdata, h->GetData(18) );//  PID 
		
	gd->SetHitChamber(ngdata,  h->GetDetID()*10+h->GetData(0)-1);
	
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
g4sbsgendata::g4sbsgendata():g4sbshitdata(-1, __GENERATED_SIZE+1){
    SetData(7,1.0);
}

#endif//__CINT__
