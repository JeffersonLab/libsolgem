void DigDemo()
{
  gROOT->LoadMacro("libsolgem.so");

  TSolSpec* s = new TSolSpec ("spectrometer", "SOLiD spectrometer");
  s->Init();

  TSolGEMChamber* y = new TSolGEMChamber ("testchamber","Test chamber");
  y->SetName("testchamber");
  y->Init();

  s->AddGEM (*y);
  s->Print();

  TSolSimGEMDigitization* d = new TSolSimGEMDigitization (*s);
  d->Print();

  TSolGEMData* gd = new TSolGEMData (1); // 1 hit wonder
  gd->SetRun (1000);
  gd->SetEvent (1);
  gd->SetMomentum (0, TVector3 (0.01, 0.02, 3.00));
  gd->SetHitEntrance (0, TVector3 (-0.01, 0.003, 1.55));
  gd->SetHitExit (0, TVector3 (-0.02, 0.001, 1.58));
  gd->SetHitEntrance (0, TVector3 (-0.022, 0.0008, 1.59));
  gd->SetHitEnergy (0, 3.1);
  gd->SetHitChamber (0, 1);
  gd->SetEntryNumber (0, 2);
  gd->SetParticleID (0, 1);
  gd->SetParticleType (0, 0);
  gd->Print();
  gd->PrintHit (0);
}
