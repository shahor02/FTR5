
// You have to load the class before ... ;-)
// .L DetectorK.cxx++

//void standardPlots() {

void testDetectorR5() {


  R5Detector its("ALICE","ITS");
  its.SetBField(1.);
  // new ideal Pixel properties?
  Double_t x0IB     = 0.001;
  Double_t x0OB     = 0.005;
  Double_t xRho     = 0.;
  Double_t resRPhiIB     = 0.0001;
  Double_t resZIB        = 0.0001;
  Double_t resRPhiOB     = 0.0005;
  Double_t resZOB        = 0.0005;
  Double_t eff           = 0.98;
  //
  // select Z span in such a way to have +-1 unit eta coverage for vertex at 2sigmaZ (~12cm) from nominal IP
  // i.e. Zmax >= 12 + R/tan( 2*atan(exp(-1.)) ) = 12 + R/0.851
  
  its.AddLayer((char*)"vertex",  0.0,  0.1, 0, 0); // dummy vertex for matrix calculation

  its.AddLayer((char*)"bpipe",   1.6,  200., 0.0022);
  its.AddLayer((char*)"ddd1",    1.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff); 
  its.AddLayer((char*)"ddd2",    2.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff); 
  its.AddLayer((char*)"ddd3",    3.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff);
  its.AddLayer((char*)"ddd3a",   8.0,  21.0, x0IB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"ddd4",   20.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"ddd5",   25.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff); 

  //its.AddLayer((char*)"ddd6",  35.0, 80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"ddd7",   40.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"ddd8",   55.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  
  //  its.AddLayer((char*)"dddZ",  90., 130., x0OB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"dddY",   80.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  its.AddLayer((char*)"dddX",  100.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff); 

  its.ProcessTrack(0.3, 0.3, 0.14, 1, 0., 4, 0, 4);

}
