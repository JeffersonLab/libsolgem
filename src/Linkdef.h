#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ defined_in "src/TSolAnalyzer.h";
#pragma link C++ defined_in "src/TSolEvData.h";
#pragma link C++ defined_in "src/TSolEVIOFile.h";
#pragma link C++ defined_in "src/TSolGEMChamber.h";
#pragma link C++ defined_in "src/TSolGEMCluster.h";
#pragma link C++ defined_in "src/TSolGEMData.h";
#pragma link C++ defined_in "src/TSolGEMPlane.h";
#pragma link C++ defined_in "src/TSolGEMVStrip.h";
#pragma link C++ defined_in "src/TSolSimAux.h";
#pragma link C++ defined_in "src/TSolSimGEMDigitization.h";
#pragma link C++ defined_in "src/TSolSpec.h";
#pragma link C++ defined_in "src/TSolWedge.h";

// Limited stuff in EVIO file.  We don't want to be
// able to call the evio functions in the interpreter

#endif
