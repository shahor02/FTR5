#include "TH1F.h"
#include "TRandom.h"
#include "TMath.h"

TH1F* dist = 0;
void tp(int ntst = 10000, double rho = 0.1, double smax = 0.4)
{
  dist = new TH1F("dist","dist",1000,0.,smax);
  // poisson mean
  double mu = smax*rho;
  
  for (int i=0;i<ntst;i++) {
    int n = gRandom->Poisson(mu);
    if (n<1) continue;
    double dmin = smax;
    for (int j=n;j--;) {
      double d = gRandom->Rndm()*smax; 
      if (d<dmin) dmin = d;
    }
    dist->Fill(dmin);
  }
}


Bool_t genFake(float rho, float sigA, float sigB, float sinT, float cosT, float nsig ,float &xw, float &yw)
{
  // generate fake whithin ellipsoid with axes sigA, sigB
  // assuming uniform background distribution
  //
  //  sigA = sigB = ss = -1;
  if (rho<1e-6) return kFALSE;
  // 1) probability to have a fake at all is Poisson distributed
  // with mu proportional to area of the ellipsoid
  float sTot = 0.5*sigA*sigB*nsig*nsig;
  double muPoisson = rho*sTot;
  int nCand = gRandom->Poisson(muPoisson);
  if (nCand<1) return kFALSE;
  // 2) generate randomly candidates and chose closest
  const int kBuff = 100;
  static float rndmArr[200];
  static int buffCnt = kBuff;
  float sw = 1.,aw=1.,bw=1.;
  while(nCand--) {
    if (buffCnt>=kBuff) {
      buffCnt = 0;
      gRandom->RndmArray(2*kBuff,rndmArr);
    }
    int ind = 2*buffCnt;
    float a = 2*(rndmArr[ind]-0.5), b = 2*(rndmArr[ind+1]-0.5);
    buffCnt++;
    float s = a*a+b*b;
    if (s<sw) {
      sw = s;
      aw = a;
      bw = b;
    }
  }
  // estimate axes of ellipse
  aw *= sigA*nsig;
  bw *= sigB*nsig;
  // rotate back to lab frame using provided sin and cos of rotated elipsoid
  xw = cosT*aw + sinT*bw;
  yw =-sinT*aw + cosT*bw;
  
  return kTRUE;
}
