void PlaneDemo()
{
  gROOT->LoadMacro("libsolgem.so");
  TSolGEMChamber y("testchamber","Test chamber");
  y.SetName("testchamber");
  y.InitPlane(0,"testplanex","Test plane X");
  y.GetPlane(0).SetName("testplanex");
  y.GetPlane(0).Init();
  y.InitPlane(1,"testplaney","Test plane Y");
  y.GetPlane(1).SetName("testplaney");
  y.GetPlane(1).Init();

  y.GetPlane(0).Print();
  y.GetPlane(1).Print();
}
