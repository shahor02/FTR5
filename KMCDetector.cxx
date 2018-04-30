//Bool_t KMCDetector::SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, Double_t phi,
//KMCProbe* KMCDetector::PrepareKalmanTrack(double pt, double lambda, double mass, int charge, double phi, double x,double y, double z)

#include "KMCDetector.h"
#include <TMath.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFormula.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TText.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TH2.h>
#include <TH1.h>
#include <TArrayI.h>
#include <AliLog.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

***********************************************************/


#define RIDICULOUS 999999 // A ridiculously large resolution (cm) to flag a dead detector

#define Luminosity    1.e27       // Luminosity of the beam (LHC HI == 1.e27, RHIC II == 8.e27 )
#define SigmaD        6.0         // Size of the interaction diamond (cm) (LHC = 6.0 cm)
#define dNdEtaMinB    1//950//660//950           // Multiplicity per unit Eta  (AuAu MinBias = 170, Central = 700)
// #define dNdEtaCent    2300//15000 //1600//2300        // Multiplicity per unit Eta  (LHC at 5.5 TeV not known)

#define CrossSectionMinB         8    // minB Cross section for event under study (PbPb MinBias ~ 8 Barns)
#define AcceptanceOfTpcAndSi     1 //1//0.60 //0.35  // Assumed geometric acceptance (efficiency) of the TPC and Si detectors
#define UPCBackgroundMultiplier  1.0   // Increase multiplicity in detector (0.0 to 1.0 * UPCRate ) (eg 1.0)
#define OtherBackground          0.0   // Increase multiplicity in detector (0.0 to 1.0 * minBias)  (eg 0.0)
#define EfficiencySearchFlag     2     // Define search method:
                                       // -> ChiSquarePlusConfLevel = 2, ChiSquare = 1, Simple = 0.  

#define PionMass                 0.139  // Mass of the Pion
#define KaonMass                 0.498  // Mass of the Kaon
#define D0Mass                   1.865  // Mass of the D0

//TMatrixD *probKomb; // table for efficiency kombinatorics

ClassImp(KMCProbe)

Int_t    KMCProbe::fgNLayers = 0;
Double_t KMCProbe::fgMissingHitPenalty = 2.;


//__________________________________________________________________________
KMCProbe::KMCProbe() :
  fMass(0.14),
  fChi2(0),
  fHits(0),
  fFakes(0),
  fNHits(0),
  fNHitsFake(0),
  fInnerChecked(-1),
  fOuterChecked(-1)
{
}

//__________________________________________________________________________
KMCProbe& KMCProbe::operator=(const KMCProbe& src) 
{
  if (this!=&src) {
    AliExternalTrackParam::operator=(src);
    fMass = src.fMass;
    fChi2 = src.fChi2;
    fHits = src.fHits;
    fFakes = src.fFakes;
    fNHits = src.fNHits;
    fNHitsFake = src.fNHitsFake;
    fInnerChecked     = src.fInnerChecked;
    fOuterChecked     = src.fOuterChecked;
  }
  return *this;
}

//__________________________________________________________________________
KMCProbe::KMCProbe(KMCProbe& src) 
  : AliExternalTrackParam(src),
    fMass(src.fMass),
    fChi2(src.fChi2),
    fHits(src.fHits),
    fFakes(src.fFakes),
    fNHits(src.fNHits),    
    fNHitsFake(src.fNHitsFake),
    fInnerChecked(src.fInnerChecked),
    fOuterChecked(src.fOuterChecked)
{
}

//__________________________________________________________________________
void KMCProbe::Reset()
{
  fMass = 0.14;
  fChi2 = 0;
  fHits = fFakes = 0;
  fNHits = fNHitsFake = 0;
  fInnerChecked = fOuterChecked = -1;
  AliExternalTrackParam::Reset();
}


//__________________________________________________________________________
void KMCProbe::ResetCovMat()
{
  // reset errors
  double *trCov  = (double*)GetCovariance();
  double *trPars = (double*)GetParameter();
  const double kLargeErr2Coord = 50*50;
  const double kLargeErr2Dir = 0.6*0.6;
  const double kLargeErr2PtI = 0.5*0.5;
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  //
}

//__________________________________________________________________________
void KMCProbe::Print(Option_t* option) const
{
  printf("M=%.3f Chi2=%7.2f (Norm:%6.2f) Checked: %d-%d, Hits: Total:%d Fakes:%d | Y:%+8.4f Z: %+8.4f |", 
	 fMass,fChi2,GetNormChi2(kTRUE),fInnerChecked,fOuterChecked, fNHits,fNHitsFake, GetY(),GetZ());
  for (int i=0;i<fgNLayers;i++) {
    if (!IsHit(i)) printf(".");
    else printf("%c",IsHitFake(i) ? '-':'+');
  }
  printf("|%s\n",IsKilled() ? " KILLED":"");
  TString opt = option;
  if (!opt.IsNull()) AliExternalTrackParam::Print(option);
}

//__________________________________________________________________________
Bool_t KMCProbe::GetXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir) const
{
  // Get local X of the track position estimated at the radius lab radius r. 
  // The track curvature is accounted exactly
  //
  // The flag "dir" can be used to remove the ambiguity of which intersection to take (out of 2 possible)
  // 0  - take the intersection closest to the current track position
  // >0 - go along the track (increasing fX)
  // <0 - go backward (decreasing fX)
  //
  // special case of R=0
  if (r<kAlmost0) {x=0; return kTRUE;}

  const double* pars = GetParameter();
  const Double_t &fy=pars[0], &sn = pars[2];
  //
  double fx = GetX();
  double crv = GetC(bz);
  if (TMath::Abs(crv)<=kAlmost0) { // this is a straight track
    if (TMath::Abs(sn)>=kAlmost1) { // || to Y axis
      double det = (r-fx)*(r+fx);
      if (det<0) return kFALSE;     // does not reach raduis r
      x = fx;
      if (dir==0) return kTRUE;
      det = TMath::Sqrt(det);
      if (dir>0) {                       // along the track direction
	if (sn>0) {if (fy>det)  return kFALSE;} // track is along Y axis and above the circle
	else      {if (fy<-det) return kFALSE;} // track is against Y axis amd belo the circle
      }
      else if(dir>0) {                                    // agains track direction
	if (sn>0) {if (fy<-det) return kFALSE;} // track is along Y axis
        else if (fy>det)  return kFALSE;        // track is against Y axis
      }
    }
    else if (TMath::Abs(sn)<=kAlmost0) { // || to X axis
      double det = (r-fy)*(r+fy);
      if (det<0) return kFALSE;     // does not reach raduis r
      det = TMath::Sqrt(det);
      if (!dir) {
	x = fx>0  ? det : -det;    // choose the solution requiring the smalest step
	return kTRUE;
      }
      else if (dir>0) {                    // along the track direction
	if      (fx > det) return kFALSE;  // current point is in on the right from the circle
	else if (fx <-det) x = -det;       // on the left
	else               x =  det;       // within the circle
      }
      else {                               // against the track direction
	if      (fx <-det) return kFALSE;  
	else if (fx > det) x =  det;
	else               x = -det;
      }
    }
    else {                                 // general case of straight line
      double cs = TMath::Sqrt((1-sn)*(1+sn));
      double xsyc = fx*sn-fy*cs;
      double det = (r-xsyc)*(r+xsyc);
      if (det<0) return kFALSE;    // does not reach raduis r
      det = TMath::Sqrt(det);
      double xcys = fx*cs+fy*sn;
      double t = -xcys;
      if (dir==0) t += t>0 ? -det:det;  // chose the solution requiring the smalest step
      else if (dir>0) {                 // go in increasing fX direction. ( t+-det > 0)
	if (t>=-det) t += -det;         // take minimal step giving t>0
	else return kFALSE;             // both solutions have negative t
      }
      else {                            // go in increasing fx direction. (t+-det < 0)
	if (t<det) t -= det;            // take minimal step giving t<0
	else return kFALSE;             // both solutions have positive t
      }
      x = fx + cs*t;
    }
  }
  else {                                 // helix
    // get center of the track circle
    double tR = 1./crv;   // track radius (for the moment signed)
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double r0 = TMath::Sqrt(x0*x0+y0*y0);
    //    printf("Xc:%+e Yc:%+e\n",x0,y0);
    //
    if (r0<=kAlmost0) {
      AliDebug(2,Form("r0 = %f",r0));
      return kFALSE;
    }            // the track is concentric to circle
    tR = TMath::Abs(tR);
    double tR2r0 = tR/r0;
    double g = 0.5*(r*r/(r0*tR) - tR2r0 - 1./tR2r0);
    double det = (1.-g)*(1.+g);
    if (det<0) {
      AliDebug(2,Form("g=%f tR=%f r0=%f\n",g,tR, r0));
      return kFALSE;
    }         // does not reach raduis r
    det = TMath::Sqrt(det);
    //
    // the intersection happens in 2 points: {x0+tR*C,y0+tR*S} 
    // with C=f*c0+-|s0|*det and S=f*s0-+c0 sign(s0)*det
    // where s0 and c0 make direction for the circle center (=x0/r0 and y0/r0)
    //
    double tmp = 1.+g*tR2r0;
    x = x0*tmp; 
    double y = y0*tmp;
    if (TMath::Abs(y0)>kAlmost0) { // when y0==0 the x,y is unique
      double dfx = tR2r0*TMath::Abs(y0)*det;
      double dfy = tR2r0*x0*TMath::Sign(det,y0);
      if (dir==0) {                    // chose the one which corresponds to smallest step 
	double delta = (x-fx)*dfx-(y-fy)*dfy; // the choice of + in C will lead to smaller step if delta<0
	if (delta<0) x += dfx;
	else         x -= dfx;
      }
      else if (dir>0) {  // along track direction: x must be > fx
	x -= dfx; // try the smallest step (dfx is positive)
	if (x<fx && (x+=dfx+dfx)<fx) return kFALSE;
      }
      else { // backward: x must be < fx
	x += dfx; // try the smallest step (dfx is positive)
	if (x>fx && (x-=dfx+dfx)>fx) return kFALSE;
      }
    }
    else { // special case: track touching the circle just in 1 point
      if ( (dir>0&&x<fx) || (dir<0&&x>fx) ) return kFALSE; 
    }
  }
  //
  return kTRUE;
}

//____________________________________
Bool_t KMCProbe::PropagateToR(double r, double b, int dir) 
{
  // go to radius R
  //
  double xR = 0;
  double rr = r*r;
  int iter = 0;
  const double kTiny = 1e-4;
  while(1) {
    if (!GetXatLabR(r ,xR, b, dir)) {
      //      printf("Track with pt=%f cannot reach radius %f\n",Pt(),r);
      //      Print("l");
      return kFALSE;
    }
    
    if (!PropagateTo(xR, b)) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to propagate to X=%f for R=%f\n",xR,r); 
	Print("l"); 
      }
      return kFALSE;
    }
    double rcurr2 = xR*xR + GetY()*GetY();
    if (TMath::Abs(rcurr2-rr)<kTiny || rr<kAlmost0) return kTRUE;
    //
    // two radii correspond to this X...
    double pos[3]; GetXYZ(pos);
    double phi = TMath::ATan2(pos[1],pos[0]);
    if (!Rotate(phi)) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to rotate to %f to propagate to R=%f\n",phi,r); 
	Print("l"); 
      }
      return kFALSE;
    }
    if (++iter>8) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to propagate to R=%f after %d steps\n",r,iter); 
	Print("l"); 
      }
      return kFALSE;
    }
  } 
  return kTRUE;
}


//__________________________________________________________________________
Bool_t KMCProbe::CorrectForMeanMaterial(const KMCLayer* lr, Bool_t inward)
{
  //  printf("before at r=%.1f p=%.4f\n",lr->fR, P());
  if (AliExternalTrackParam::CorrectForMeanMaterial(lr->fx2X0, inward ? lr->fXRho : -lr->fXRho, GetMass() , kTRUE)) {
    //  printf("after  at r=%.1f p=%.4f\n",lr->fR, P());
    return kTRUE;
  }
  AliDebug(2,Form("Failed to apply material correction, X/X0=%.4f", lr->fx2X0));
  if (AliLog::GetGlobalDebugLevel()>1) Print();
  return kFALSE;
}

/////////////////////////////////////////////////////////////////////////////
ClassImp(KMCCluster)

//_________________________________________________________________________
KMCCluster::KMCCluster(KMCCluster &src) 
: TObject(src),
  fY(src.fY),fZ(src.fZ),fX(src.fX),fPhi(src.fPhi)
{}

//__________________________________________________________________________
KMCCluster& KMCCluster::operator=(const KMCCluster& src) 
{
  if (this!=&src) {
    TObject::operator=(src);
    fY = src.fY;
    fZ = src.fZ;
    fX = src.fX;
    fPhi = src.fPhi;
  }
  return *this;
}

//_________________________________________________________________________
void KMCCluster::Print(Option_t *) const 
{
  printf(" Local YZ = (%3.4lf,%3.4lf) | X=%3.4lf  phi: %+.3f %s\n",fY,fZ,fX,fPhi,IsKilled()?"Killed":""); 
}

/////////////////////////////////////////////////////////////////////////////
ClassImp(KMCLayer)

//__________________________________________________________________________
KMCLayer::KMCLayer(char *name) : 
  TNamed(name,name),fR(0),fx2X0(0),fPhiRes(0),fZRes(0),fEff(0),fIsDead(kFALSE),fActiveID(-1),fClMC()
{
  Reset();
}

//__________________________________________________________________________
void KMCLayer::Print(Option_t *opt) const
{
  printf("Lr%3d(A%3d) %10s R=%5.1f X2X0=%.3f XRho=%.3f SigY=%.4f SigZ=%.4f Eff:%4.2f\n",
	 GetUniqueID(),fActiveID,GetName(), fR, fx2X0,fXRho,fPhiRes,fZRes,fEff);
  TString opts = opt; opts.ToLower();
  if (opts.Contains("c")) {
    printf("Cluster: MC: %+7.4f:%+7.4f\n",fClMC.fY,fClMC.fZ);
  }
}

/////////////////////////////////////////////////////////////////////////////
Double_t KMCDetector::fgVtxConstraint[2]={-1,-1};

ClassImp(KMCDetector)

KMCDetector::KMCDetector() :
TNamed("test_detector","detector"),
  fNLayers(0),
  fNActiveLayers(0),
  fFirstActiveLayer(-1),
  fLastActiveLayer(-1),
  fFirstActiveLayerTracked(-1),
  fLastActiveLayerTracked(-1),
  fBFieldG(5.),
  fIntegrationTime(0.02), // in ms
  fdNdEtaCent(2000),       // Multiplicity
  fDensFactorEta(1),
  fMaxChi2Cl(30.),
  fMaxNormChi2NDF(7.),
  fMinHits(4),
  fMaxSnp(0.8)
{
  //
  // default constructor
  //
}

KMCDetector::KMCDetector(char *name, char *title)
  : TNamed(name,title),
    fNLayers(0),
    fNActiveLayers(0),
    fFirstActiveLayer(-1),
    fLastActiveLayer(-1),
    fFirstActiveLayerTracked(-1),
    fLastActiveLayerTracked(-1),
    fBFieldG(5.),
    fIntegrationTime(0.02), // in ms
    fdNdEtaCent(2000),       // Multiplicity
    fDensFactorEta(1.),
    fMaxChi2Cl(30.),
    fMaxNormChi2NDF(7.),
    fMinHits(4),
    fMaxSnp(0.8)

{
  //
  // default constructor, that set the name and title
  //
  //  fLayers = new TObjArray();
}

KMCDetector::~KMCDetector() { // 
  // virtual destructor
  //
  //  delete fLayers;
}

void KMCDetector::AddLayer(char *name, Float_t radius, Float_t x2X0, Float_t xrho, Float_t phiRes, Float_t zRes, Float_t eff) {
  //
  // Add additional layer to the list of layers (ordered by radius)
  // 

  KMCLayer *newLayer = (KMCLayer*) fLayers.FindObject(name);

  if (!newLayer) {
    newLayer = new KMCLayer(name);
    newLayer->fR = radius;
    newLayer->fx2X0 = x2X0;
    newLayer->fXRho  = xrho;
    newLayer->fPhiRes = phiRes;
    newLayer->fZRes = zRes;
    eff = TMath::Min(1.f,eff);
    newLayer->fEff = eff;
    newLayer->fActiveID = -2;
    TString lname = name;
    if (lname.Contains("vertex")) newLayer->SetBit(KMCLayer::kBitVertex);
    //
    newLayer->fIsDead =  (newLayer->fPhiRes<0 && newLayer->fZRes<0) || newLayer->fEff<=0.;
    //
    if (fLayers.GetEntries()==0) 
      fLayers.Add(newLayer);
    else {
      //
      for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
	KMCLayer *l = (KMCLayer*)fLayers.At(i);
	if (radius<l->fR) { fLayers.AddBefore(l,newLayer); break; }
	  if (radius>l->fR && (i+1)==fLayers.GetEntries() ) fLayers.Add(newLayer); // even bigger then last one
      }
      //
    }
    //
    ClassifyLayers();
    fNLayers++;
    //
  } else {
    printf("Layer with the name %s does already exist\n",name);
  }
}

//____________________________________________________________
void KMCDetector::ClassifyLayers()
{
  // assign active Id's, etc
  fFirstActiveLayer = fLastActiveLayer = -1;

  fNActiveLayers = 0;
  //
  int nl = GetNLayers();
  for (int il=0;il<nl;il++) {
    KMCLayer* lr = GetLayer(il);
    lr->SetUniqueID(il);
    if (!lr->IsDead()) {
      fLastActiveLayer = il; 
      if (fFirstActiveLayer<0) fFirstActiveLayer = il;
      lr->SetActiveID(fNActiveLayers++);
    }
  }
  //
  KMCProbe::SetNLayers(fNActiveLayers);
}


//________________________________________________________________________________
Int_t KMCDetector::GetLayerID(int actID) const
{
  // find physical layer id from active id
  if (actID<0 || actID>fNActiveLayers) return -1;
  for (int i=fLastActiveLayer; i--;) {
    if (GetLayer(i)->GetActiveID()==actID) return i;   
  }
  return -1;
}

void KMCDetector::PrintLayout() {
  //
  // Prints the detector layout
  //
  printf("Detector %s: \"%s\"\n",GetName(),GetTitle());
  
  if (fLayers.GetEntries()>0)  printf("  Name \t\t r [cm] \t  X0 \t  phi & z res [um]\n");

  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    KMCLayer* tmp = (KMCLayer*)fLayers.At(i);
    
    printf("%d. %s \t %03.2f   \t%1.4f\t  ",i, tmp->GetName(), tmp->GetRadius(), tmp->GetRadL() );
    if (tmp->IsDead()) printf("  -  ");
    else               printf("%3.0f   ",tmp->GetPhiRes()*10000);
    if (tmp->IsDead()) printf("  -\n");
    else               printf("%3.0f\n",tmp->GetZRes()*10000);
  }
}

Double_t KMCDetector::ThetaMCS ( Double_t mass, Double_t x2X0, Double_t momentum ) const
{
  //
  // returns the Multiple Couloumb scattering angle (compare PDG boolet, 2010, equ. 27.14)
  //

  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( x2X0 > 0 ) theta  =  0.0136 * TMath::Sqrt(x2X0) / ( beta * momentum ) * (1+0.038*TMath::Log(x2X0)) ;
  return (theta) ;
}

Double_t KMCDetector::ProbGoodHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Howard Wieman: http://rnc.lbl.gov/~wieman/GhostTracks.htm 
  // and http://rnc.lbl.gov/~wieman/HitFinding2D.htm
  // This is the probability of getting a good hit using 2D Gaussian distribution function and infinite search radius
  Double_t sx, sy, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusRPhi * HitDensity(radius) ;
  sy = 2 * TMath::Pi() *  searchRadiusZ    * searchRadiusZ    * HitDensity(radius) ;
  goodHit =  TMath::Sqrt(1./((1+sx)*(1+sy)))  ;
  return ( goodHit ) ;
}


Double_t KMCDetector::ProbGoodChiSqHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Victor Perevoztchikov and Howard Wieman: http://rnc.lbl.gov/~wieman/HitFinding2DXsq.htm
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  Double_t sx, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusZ * HitDensity(radius) ;
  goodHit =  1./(1+sx) ;
  return ( goodHit ) ;  
}

Double_t KMCDetector::ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff,
						 Double_t searchRadiusRPhi, Double_t searchRadiusZ, double confLevel) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  // Plus, in addition, taking a "confidence level" and the "layer efficiency" into account 
  // Following is correct for 2 DOF

  Double_t c = -2 *TMath::Log(confLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t goodHit = leff/(2*alpha) * (1 - TMath::Exp(-alpha*c));
  return ( goodHit ) ;  
}

Double_t KMCDetector::ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff,
						 Double_t searchRadiusRPhi, Double_t searchRadiusZ, double confLevel) 
{
  // Based on work by Ruben Shahoyan 
  // This is the probability to not have any match to the track (see also :ProbGoodChiSqPlusConfHit:)

  Double_t c = -2 *TMath::Log(confLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t nullHit = (1-leff+confLevel*leff)*TMath::Exp(-c*(alpha-1./2));
  return ( nullHit ) ;  
}


Double_t KMCDetector::HitDensity ( Double_t radius ) 
{
  // Background (0-1) is included via 'OtherBackground' which multiplies the minBias rate by a scale factor.
  // UPC electrons is a temporary kludge that is based on Kai Schweda's summary of Kai Hainken's MC results
  // See K. Hencken et al. PRC 69, 054902 (2004) and PPT slides by Kai Schweda.
  // Note that this function assumes we are working in CM and CM**2 [not meters].
  // Based on work by Yan Lu 12/20/2006, all radii and densities in centimeters or cm**2.

  //  Double_t MaxRadiusSlowDet = 0.1; //?   // Maximum radius for slow detectors.  Fast detectors 
  if (radius<0.01) return 0;
  //
  double arealDensity  = OneEventHitDensity(fdNdEtaCent,radius)  + UpcHitDensity(radius) ;
  return ( arealDensity ) ;  
}

double KMCDetector::OneEventHitDensity( Double_t multiplicity, Double_t radius ) const
{
  // This is for one event at the vertex.  No smearing.
  double den   = multiplicity / (2.*TMath::Pi()*radius*radius) * fDensFactorEta ; // 2 eta ?
  // note: surface of sphere is  '4*pi*r^2'
  //       surface of cylinder is '2*pi*r* h' 
  return den ;
} 

double KMCDetector::UpcHitDensity(Double_t radius)
{ 
  // QED electrons ...

  Double_t mUPCelectrons = 0;
  /*
  //  mUPCelectrons =  fLhcUPCscale * (1.23 - radius/6.5)      ;  // Fit to Kai Schweda summary tables at RHIC * 'scale' for LHC
  mUPCelectrons = fLhcUPCscale*5456/(radius*radius)/dNdEtaMinB;      // Fit to 'Rossegger,Sadovsky'-Alice simulation
  if ( mUPCelectrons < 0 ) mUPCelectrons =  0.0             ;  // UPC electrons fall off quickly and don't go to large R
  mUPCelectrons *= IntegratedHitDensity(dNdEtaMinB,radius) ;  // UPCs increase Mulitiplicty ~ proportional to MinBias rate
  mUPCelectrons *= UPCBackgroundMultiplier                 ;  // Allow for an external multiplier (eg 0-1) to turn off UPC
  */
  return mUPCelectrons ;
}

void KMCDetector::CalcDensFactorEta(double eta)
{
  if (TMath::Abs(eta)<1e-3) fDensFactorEta = 1.;
  else {
    fDensFactorEta = TMath::Tan( 2.*TMath::ATan(TMath::Exp(-TMath::Abs(eta))) );
    fDensFactorEta = 1./TMath::Sqrt( 1. + 1./fDensFactorEta/fDensFactorEta);
  }
}

void KMCDetector::ApplyMS(KMCProbe* trc, double x2X0) const
{
  // simulate random modification of track params due to the MS
  if (x2X0<=0) return;
  double alpha = trc->GetAlpha(); // store original alpha
  double mass = trc->GetMass();
  //
  double snp = trc->GetSnp();
  double dip = trc->GetTgl();
  Double_t angle=TMath::Sqrt((1.+ dip*dip)/((1-snp)*(1.+snp)));
  x2X0 *= angle;
  //
  static double covCorr[15],covDum[21]={0};
  static double mom[3],pos[3];
  double *cov = (double*) trc->GetCovariance();
  memcpy(covCorr,cov,15*sizeof(double));
  trc->GetXYZ(pos);
  trc->GetPxPyPz(mom);
  double pt2 = mom[0]*mom[0]+mom[1]*mom[1];
  double pt = TMath::Sqrt(pt2);
  double ptot2 = pt2 + mom[2]*mom[2];
  double ptot  = TMath::Sqrt(ptot2);
  double beta = ptot/TMath::Sqrt(ptot2 + mass*mass);
  double sigth = TMath::Sqrt(x2X0)*0.014/(ptot*beta);
  //
  // a la geant
  double phiSC = gRandom->Rndm()*TMath::Pi();
  double thtSC = gRandom->Gaus(0,1.4142*sigth);
  //  printf("MS phi: %+.5f tht: %+.5f\n",phiSC,thtSC);
  double sn = TMath::Sin(thtSC);
  double dx = sn*TMath::Sin(phiSC);
  double dy = sn*TMath::Cos(phiSC);  
  double dz = TMath::Cos(thtSC);
  double v[3];
  //  printf("Before: %+.3e %+.3e %+.3e | MS: %+.3e %+.3e\n",mom[0],mom[1],mom[2],thtSC,phiSC);
  for (int i=3;i--;) mom[i] /= ptot;
  double vmm = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  if (!IsZero(pt)) {
    double pd1 = mom[0]/vmm;
    double pd2 = mom[1]/vmm;
    v[0] = pd1*mom[2]*dx - pd2*dy + mom[0]*dz;
    v[1] = pd2*mom[2]*dx + pd1*dy + mom[1]*dz;
    v[2] = -vmm*dx                + mom[2]*dz;
  }
  else {
    v[0] = dx;
    v[1] = dy;
    v[2] = dz*TMath::Sign(1.,mom[2]);
  }
  double nrm = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //  printf("before :%+e %+e %+e  || %+e %+e %+e %+e\n",mom[0],mom[1],mom[2],  sigth, x2X0, pt, beta);
  //  trc->Print();
  // direction cosines -> p
  for (int i=3;i--;) mom[i] = ptot*v[i]/nrm;
  //  printf("After : %+.3e %+.3e %+.3e\n",mom[0],mom[1],mom[2]);
  trc->Set(pos,mom,covDum,trc->Charge());
  //
  trc->Rotate(alpha);
  memcpy(cov,covCorr,15*sizeof(double));
  //
}

KMCProbe* KMCDetector::PrepareKalmanTrack(double pt, double eta, double mass, int charge, double phi, double x,double y, double z)
{
  // Prepare trackable Kalman track at the farthest position
  //
  // Set track parameters
  // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
  double lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-eta)); 
  fProbe.Reset();
  fProbe.SetMass(mass);
  KMCProbe* probe = new KMCProbe(fProbe);
  double *trPars = (double*)probe->GetParameter();
  double *trCov  = (double*)probe->GetCovariance();
  double xyz[3] = {x,y,z};
  probe->Global2LocalPosition(xyz,phi);
  probe->Set(xyz[0],phi,trPars,trCov);
  trPars[KMCProbe::kY] = xyz[1];
  trPars[KMCProbe::kZ] = xyz[2];
  trPars[KMCProbe::kSnp] = 0;                       //            track along X axis at the vertex
  trPars[KMCProbe::kTgl] = TMath::Tan(lambda);                // dip
  trPars[KMCProbe::kPtI] = charge/pt;               //            q/pt      
  //
  // put tiny errors to propagate to the outer-most radius
  trCov[KMCProbe::kY2] = trCov[KMCProbe::kZ2] = trCov[KMCProbe::kSnp2] = trCov[KMCProbe::kTgl2] = trCov[KMCProbe::kPtI2] = 1e-20;
  fProbe = *probe;  // store original track
  //
  Bool_t res = TransportKalmanTrackWithMS(probe);
  probe->ResetCovMat();// reset cov.matrix
  //
  return probe;
}


//________________________________________________________________________________
int KMCDetector::TransportKalmanTrackWithMS(KMCProbe *probTr, Bool_t applyMatCorr)
{
  // Transport track till layer maxLr, applying random MS
  //
  int nActLrOK = 0;
  double r = TMath::Sqrt(probTr->GetX()*probTr->GetX()+probTr->GetY()*probTr->GetY());
  fFirstActiveLayerTracked = -1;
  fLastActiveLayerTracked = 0;
  for (Int_t j=0; j<fNLayers; j++) {
    if (j>fLastActiveLayer) continue;
    KMCLayer* lr = (KMCLayer*)fLayers.At(j);
    if (lr->GetRadius() <= r) continue;
    if (!PropagateToLayer(probTr,lr,1)) break;
    if (TMath::Abs(probTr->GetSnp())>fMaxSnp) break;
    if (lr->GetRadL()>0 && applyMatCorr) {
      ApplyMS(probTr,lr->GetRadL()); // apply MS
      if (!probTr->CorrectForMeanMaterial(lr,kFALSE)) break; 
    }
    //
    if (lr->IsDead()) continue;
    if (fFirstActiveLayerTracked<0) fFirstActiveLayerTracked = lr->GetActiveID();
    fLastActiveLayerTracked = lr->GetActiveID();
    
    // store randomized cluster local coordinates and phi
    double rz,ry;
    gRandom->Rannor(rz,ry);
    lr->GetMCCluster()->Set(probTr->GetY()+ry*lr->GetPhiRes(),probTr->GetZ()+rz*lr->GetZRes(), 
			    probTr->GetX(), probTr->GetAlpha() );
    nActLrOK++;
    //
  }
  //
  return nActLrOK;
}

//____________________________________________________________________________
Bool_t KMCDetector::PropagateToLayer(KMCProbe* trc, KMCLayer* lr, int dir) const
{
  // bring the track to layer and rotat to frame normal to its surface
  if (!trc->PropagateToR(lr->fR,fBFieldG, dir)) return kFALSE;
  //
  // rotate to frame with X axis normal to the surface (defined by ideal track)
  if (!lr->IsVertex()) {
    double phi = trc->PhiPos();
    if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;
    if (!trc->Rotate(phi)) {
      return kFALSE;
    }
  }
  //
  return kTRUE;
}

//____________________________________________________________________________
Bool_t KMCDetector::UpdateTrack(KMCProbe* trc, KMCLayer* lr, KMCCluster* cl) const
{
  // update track with measured cluster
  // propagate to cluster
  double meas[2] = {cl->GetY(),cl->GetZ()}; // ideal cluster coordinate
  double measErr2[3] = {lr->fPhiRes*lr->fPhiRes,0,lr->fZRes*lr->fZRes};
  //
  if (!trc->PropagateToCluster(cl,fBFieldG)) return kFALSE; // track was not propagated to cluster frame
  //
  double chi2 = trc->GetPredictedChi2(meas,measErr2);
  //  if (chi2>fMaxChi2Cl) return kTRUE; // chi2 is too large
  //  
  if (!trc->Update(meas,measErr2)) {
    AliDebug(2,Form("layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		    lr->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
    if (AliLog::GetGlobalDebugLevel()>1) trc->Print("l");
    return kFALSE;
  }
  trc->AddHit(lr->GetActiveID(), chi2);
  //
  return kTRUE;
}



///LAST

/*
//________________________________________________________________________________
Bool_t KMCDetector::SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, Double_t phi,
				     Double_t xv, Double_t yv, Double_t zv,
				     int charge)
{
  // analytic and fullMC of track with given kinematics.
  //
  // prepare kalman track
  KMCProbe* probe = PrepareKalmanTrack(pt,eta,mass,charge,phi,xv,yv,zv);

  
  if (!SolveSingleTrackViaKalman(mass,pt,eta)) return kFALSE;
  //
  // Store non-updated track errors of inward propagated seed >>>>>>>>
  int maxLr = fLastActiveITSLayer + offset;
  if (maxLr >= fLastActiveLayerTracked-1) maxLr = fLastActiveLayerTracked;
  KMCProbe probeTmp = fProbe; // original probe at vertex
  KMCLayer* lr = 0;
  for (Int_t j=1; j<=maxLr; j++) {
    lr = GetLayer(j);
    //    printf("Here0: %d\n",j);
    if (!PropagateToLayer(&probeTmp,lr,1)) return 0;
    if (j!=maxLr) if (!probeTmp.CorrectForMeanMaterial(lr, kFALSE)) return 0;
    //    printf("Prelim. Err at lr:%8s | %7.3f %7.3f\n",lr->GetName(),TMath::Sqrt(probeTmp.GetSigmaY2()),TMath::Sqrt(probeTmp.GetSigmaZ2()));
  }
  for (Int_t j=maxLr; j>0; j--) {
    lr = GetLayer(j);
    //    printf("Here1: %d\n",j);
    if (j!=maxLr) if (!PropagateToLayer(&probeTmp,lr,-1)) return 0;
    lr->fSig2EstD = probeTmp.GetSigmaY2();
    lr->fSig2EstZ = probeTmp.GetSigmaZ2();
    //    probeTmp.Print("l");
    printf("Natural Err at lr:%8s | %7.3f %7.3f\n",lr->GetName(),TMath::Sqrt(lr->fSig2EstD),TMath::Sqrt(lr->fSig2EstZ));
    if (!probeTmp.CorrectForMeanMaterial(lr, kTRUE)) return 0;
  }
  // Store non-updated track errors of inward propagated seed <<<<<<<<
  //
  int nsm = sumArr ? sumArr->GetEntriesFast() : 0;
  KMCLayer* vtx = GetLayer(0);
  //
  for (int i=0;i<nsm;i++) {
    KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(i);
    if (!tsm) continue;
    tsm->SetRefProbe( GetProbeTrack() ); // attach reference track (generated)
    tsm->SetAnProbe( vtx->GetAnProbe() ); // attach analitycal solution
  }
  //
  TStopwatch sw;
  sw.Start();
  for (int it=0;it<nMC;it++) {
    printf("ev: %d\n",it);
    SolveSingleTrackViaKalmanMC(offset);
    KMCProbe* trc = vtx->GetWinnerMCTrack();
    vtx->GetMCTracks()->Print();
    if (progressP==1 || (progressP>0 &&  (it%progressP)==0)) {
      printf("%d%% done |",it*100/nMC); 
      sw.Stop(); sw.Print(); sw.Start(kFALSE);
    }
    for (int ism=nsm;ism--;) { // account the track in each of summaries
      KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(ism);
      if (!tsm) continue;
      tsm->AddUpdCalls(GetUpdCalls());
      tsm->AddTrack(trc); 
    }
  }
  //
  sw.Stop();
  printf("Total time: "); sw.Print();
  return kTRUE;
}

//________________________________________________________________________________
KMCProbe* KMCDetector::KalmanSmooth(int actLr, int actMin,int actMax) const
{
  // estimate kalman smoothed track params at given active lr 
  // from fit at layers actMin:actMax (excluding actLr)
  // SolveSingleTrackViaKalman must have been called before
  //
  if (actMin>actMax) swap(actMin,actMax);
  if (actMax>=fNActiveLayers) actMax = fNActiveLayers-1;
  int nlrfit = actMax-actMin;
  if (actLr>=actMin && actLr<=actMax) nlrfit-=1;
  if (nlrfit<2) {AliInfo("Need a least 2 active layers in the fit"); return 0;}
  static KMCProbe iwd,owd;
  //
  // find phisical layer id's
  int pLr  = GetLayerID(actLr);
  int pMin = GetLayerID(actMin);
  int pMax = GetLayerID(actMax);
  //
  //  printf(">>> %d %d %d\n",pLr, pMin,pMax);
  Bool_t useIwd=kFALSE, useOwd=kFALSE;
  if (pLr<pMax) { // need inward piece
    iwd = GetLayer(pMax)->fTrCorr;
    iwd.ResetCovMat();
    iwd.GetHitsPatt() = 0;
    for (int i=pMax;i>=pLr;i--) {
      KMCLayer* lr = GetLayer(i);
      //      printf("IWD %d\n",i);
      if (!lr->IsDead() && i!=pLr && i>=pMin) if (!UpdateTrack(&iwd,lr,&lr->fClCorr))  return 0;
      if (i!=pLr) {
	if (!iwd.CorrectForMeanMaterial(lr,kTRUE)) return 0; // correct for materials of this layer
	if (!PropagateToLayer(&iwd,GetLayer(i-1),-1)) return 0;      // propagate to next layer
      }
      //  printf("IWD%d:  ",i); iwd.Print("l");
    }
    useIwd = kTRUE;
  }
  if (pLr>pMin) { // need outward piece
    owd = GetLayer(pMin)->fTrCorr;
    owd.ResetCovMat();
    owd.GetHitsPatt() = 0;
    for (int i=pMin;i<=pLr;i++) {
      KMCLayer* lr = GetLayer(i);
      //      printf("OWD %d\n",i);
      if (!lr->IsDead() && i!=pLr && i<=pMax) if (!UpdateTrack(&owd,lr,&lr->fClCorr))  return 0;
      if (i!=pLr) {
	if (!owd.CorrectForMeanMaterial(lr,0)) return 0; // correct for materials of this layer
	if (!PropagateToLayer(&owd,GetLayer(i+1), 1)) return 0;      // propagate to next layer
      }
      //      printf("OWD%d:  ",i); owd.Print("l");
    }
    useOwd = kTRUE;
  }
  //
  // was this extrapolation outside the fit range?
  if (!useIwd) return (KMCProbe*)&owd; 
  if (!useOwd) return (KMCProbe*)&iwd;
  //
  // weight both tracks
  if (!iwd.Propagate(owd.GetAlpha(),owd.GetX(),fBFieldG)) return 0;
  double meas[2] = {owd.GetY(),owd.GetZ()};
  double measErr2[3] = {owd.GetSigmaY2(), owd.GetSigmaZY(), owd.GetSigmaZ2()};
  //  printf("Weighting\n");
  //  owd.Print("l");
  //  iwd.Print("l");
  if (!iwd.Update(meas,measErr2)) return 0;
  iwd.GetHitsPatt() |= owd.GetHitsPatt();

  //  printf("->\n");
  //  iwd.Print("l");

  return (KMCProbe*)&iwd;
  //
}

//________________________________________________________________________________
KMCProbe* KMCDetector::KalmanSmoothFull(int actLr, int actMin,int actMax) const
{
  // estimate kalman smoothed track params at given active lr 
  // from fit at layers actMin:actMax (excluding actLr)
  // SolveSingleTrackViaKalman must have been called before
  //
  static TClonesArray prediction("KMCProbe",10);
  static TClonesArray update("KMCProbe",10);
  static KMCProbe res;
  //
  if (actMin>actMax) swap(actMin,actMax);
  int nlrfit = actMax-actMin;
  if (actLr>=actMin && actLr<=actMax) nlrfit-=1;
  if (nlrfit<2) {AliInfo("Need a least 2 active layers in the fit"); return 0;}
  //
  // find phisical layer id's
  int pLr  = GetLayerID(actLr);
  int pMin = GetLayerID(actMin);
  int pMax = GetLayerID(actMax);
  //
  int dir=0,dirInt=0;
  if      (pLr<=pMin) dir=-1; // inward extrapolation
  else if (pLr>=pMax) dir= 1; // outward extrapolation
  else if (actMax-actLr >= actLr-actMin) dirInt = -1; // inward  interpolation (the test point is closer to inner layer)
  else    dirInt = 1;                                 // outward interpolation (the test point is closer to outer layer)
  //
  if (dir!=0) { // no sens to do smoothing: simple Kalman filtering extrapolation
    int start = dir<0 ? pMax : pMin;
    res = GetLayer(start)->fTrCorr;
    res.ResetCovMat();
    KMCLayer* lr = 0;
    for (int i=(dir<0?pMax:pMin); i!=pLr; i+=dir) { // track till nearest layer to pLr
      lr = GetLayer(i);
      if (!lr->IsDead() && !(i<pMin ||i>pMax)) if (!UpdateTrack(&res,lr,&lr->fClCorr))  return 0; // update only with layers in fit range
      if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this layer
      if (!PropagateToLayer(&res,GetLayer(i+dir),dir))            return 0; // propagate to next layer     
    }
    if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this nearest layer
    if (!PropagateToLayer(&res,GetLayer(pLr), dir)) return 0; // propagate to test layer
    return (KMCProbe*)&res;
  }
  //
  // too bad, need to do real filtering
  //
  int start = dirInt<0 ? pMax : pMin;
  int stop  = dirInt<0 ? pMin-1 : pMax+1;
  res = GetLayer(start)->fTrCorr;
  res.ResetCovMat();
  KMCLayer* lr = 0;
  int count = 0;
  for (int i=start; i!=stop; i+=dirInt) { // track in full range, storing updates and predictions
    new(prediction[count]) KMCProbe(res);
    lr = GetLayer(i);
    if (!lr->IsDead() && i!=pLr) if (!UpdateTrack(&res,lr,&lr->fClCorr))  return 0; // update only with layers in fit range
    new(update[count]) KMCProbe(res);
    if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this layer
    if (!PropagateToLayer(&res,GetLayer(i+dir),dir))            return 0; // propagate to next layer     
    count++;
  }
  return (KMCProbe*)&res;
  //
}

//________________________________________________________________________________
Bool_t KMCDetector::SolveSingleTrackViaKalman(Double_t mass, Double_t pt, Double_t eta)
{
  // analytical estimate of tracking resolutions
  //  fProbe.SetUseLogTermMS(kTRUE);
  //
  if (fMinITSHits>fNActiveITSLayers) {fMinITSHits = fNActiveITSLayers; printf("Redefined request of min N ITS hits to %d\n",fMinITSHits);}
  if (TMath::Abs(eta)<1e-3) fDensFactorEta = 1.;
  else {
    fDensFactorEta = TMath::Tan( 2.*TMath::ATan(TMath::Exp(-TMath::Abs(eta))) );
    fDensFactorEta = 1./TMath::Sqrt( 1. + 1./fDensFactorEta/fDensFactorEta);
  }
  double lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-eta)); 
  KMCProbe* probe = PrepareKalmanTrack(pt,lambda,mass,-1);
  if (!probe) return kFALSE;
  //
  KMCLayer *lr = 0;
  //
  //
  // Start the track fitting --------------------------------------------------------
  //
  // Back-propagate the covariance matrix along the track. 
  // Kalman loop over the layers
  //
  KMCProbe* currTr = 0;
  lr = (KMCLayer*)fLayers.At(fLastActiveLayerTracked);
  lr->fTrCorr = *probe;
  delete probe; // rethink...
  //
  for (Int_t j=fLastActiveLayerTracked; j--; ) {  // Layer loop
    //
    KMCLayer *lrP = lr;
    lr = (KMCLayer*)fLayers.At(j);
    //
    lr->fTrCorr = lrP->fTrCorr;
    currTr = &lr->fTrCorr;
    currTr->ResetHit(lrP->GetActiveID());
    //
    // if there was a measurement on prev layer, update the track
    if (!lrP->IsDead()) { // include measurement
      KMCCluster cl(currTr->GetY(),currTr->GetZ(), currTr->GetX(), currTr->GetAlpha());
      if (!UpdateTrack(currTr,lrP,&cl))  return kFALSE;
    }
    if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) return kFALSE; // correct for materials of this layer
    if (!PropagateToLayer(currTr,lr,-1)) return kFALSE;      // propagate to current layer
    //
  } // end loop over layers
  //
  return kTRUE;
}

//____________________________________________________________
Bool_t KMCDetector::SolveSingleTrackViaKalmanMC(int offset)
{
  // MC estimate of tracking resolutions/effiencies. Requires that the SolveSingleTrackViaKalman
  // was called before, since it uses data filled by this method
  //
  // The MC tracking will be done starting from fLastActiveITSLayer + offset (before analytical estimate will be used)
  //
  // At this point, the fProbe contains the track params generated at vertex.
  // Clone it and propagate to target layer to generate hit positions affected by MS
  //
  fUpdCalls = 0.;
  KMCProbe *currTrP=0,*currTr=0;
  int maxLr = fLastActiveITSLayer + offset;
  if (maxLr >= fLastActiveLayerTracked-1) maxLr = fLastActiveLayerTracked;
  ResetMCTracks(maxLr);
  KMCLayer* lr = (KMCLayer*)fLayers.At(maxLr);
  currTr = lr->AddMCTrack(&fProbe); // start with original track at vertex
  //
  if (!TransportKalmanTrackWithMS(currTr, maxLr)) return kFALSE; // transport it to outermost layer where full MC is done
  //
  if (fLastActiveITSLayer<fLastActiveLayerTracked) { // prolongation from TPC
    // start from correct track propagated from above till maxLr
    double *covMS = (double*)currTr->GetCovariance();
    const double *covIdeal =lr->fTrCorr.GetCovariance();
    for (int i=15;i--;) covMS[i] = covIdeal[i];
  }
  else { // ITS SA: randomize the starting point
    //    double *pars = (double*)currTr->GetParameter();
    //    pars[0] += gRandom->Gaus(0,TMath::Sqrt(currTr->GetSigmaY2()));
    //    pars[1] += gRandom->Gaus(0,TMath::Sqrt(currTr->GetSigmaZ2()));
    //
    currTr->ResetCovMat();
  }
  //
  for (Int_t j=maxLr; j--; ) {  // Layer loop
    //
    KMCLayer *lrP = lr;
    lr = (KMCLayer*)fLayers.At(j);
    int ntPrev = lrP->GetNMCTracks();
    //
    if (lrP->IsDead()) { // for passive layer just propagate the copy of all tracks of prev layer >>>
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill(); continue;} // correct for materials of prev. layer
	if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of dead layer <<<
    //
    if (lrP->IsTPC()) { // we don't consider bg hits in TPC, just update with MC cluster
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (!UpdateTrack(currTr, lrP, lrP->GetMCCluster(), kTRUE)) {currTr->Kill(); continue;} // update with correct MC cl.
	if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill(); continue;} // correct for materials of prev. layer
	if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of ideal (TPC?) layer <<<
    //
    // active layer under eff. study (ITS?): propagate copy of every track to MC cluster frame (to have them all in the same frame)
    // and calculate the limits of bg generation
    KMCCluster* clMC = lrP->GetMCCluster();
    if (lrP->GetLayerEff()<gRandom->Rndm()) clMC->Kill(); // simulate inefficiency
    ResetSearchLimits();
    int nseeds = 0;
    for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
      currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
      currTr = lr->AddMCTrack( currTrP );
      if (!currTr->PropagateToCluster(clMC,fBFieldG)) {currTr->Kill(); continue;} // go to MC cluster
      if ( !(currTr->GetNITSHits()>0 && currTr->GetNITSHits()==currTr->GetNFakeITSHits()) ) UpdateSearchLimits(currTr, lrP); // RS
      nseeds++;
    }
    //
    //    printf("%3d seeds\n",nseeds);
    if (fUseBackground && lrP->IsITS()) GenBgClusters(lrP); //  generate background hits
    //
    ntPrev = lr->GetNMCTracks();
    for (int itr=ntPrev;itr--;) { // loop over all tracks PROPAGATED from previous layer to clusters frame on previous layer
      currTrP = lr->GetMCTrack(itr); // this is a seed from prev layer. The new clusters are attached to its copies, the seed itself
                                     // will be propagated w/o cluster update if it does not violate requested "reconstructed" track settings
      if (currTrP->IsKilled()) continue;
      //printf("Check    %d %p %d\n",itr,currTrP,currTrP->GetUniqueID()); currTrP->Print();
      CheckTrackProlongations(currTrP, lr, lrP);
      if (NeedToKill(currTrP)) currTrP->Kill(); // kill track which was not updated at lrP
      //currTrP->Kill(); // kill track which was not updated at lrP
    }
    //  
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
    if (ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0) {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--)  lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    //
    for (int itr=ntTot;itr--;) {
      currTr = lr->GetMCTrack(itr);
      if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill();continue;} // correct for materials of prev. layer
      if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill();continue;} // propagate to current layer
    }
    AliDebug(1,Form("Got %d tracks on layer %s",ntTot,lr->GetName()));
    //    lr->GetMCTracks()->Print();
    //
  } // end loop over layers    
  //
  // do we use vertex constraint?
  KMCLayer *vtx = GetLayer(0);
  if (!vtx->IsDead() && vtx->IsITS()) {
    int ntr = vtx->GetNMCTracks();
    for (int itr=0;itr<ntr;itr++) {
      currTr = vtx->GetMCTrack(itr);
      if (currTr->IsKilled()) continue;
      KMCCluster* clv = vtx->GetMCCluster();
      double meas[2] = {clv->GetY(),clv->GetZ()};
      double measErr2[3] = {vtx->fPhiRes*vtx->fPhiRes,0,vtx->fZRes*vtx->fZRes};
      double chi2v = currTr->GetPredictedChi2(meas,measErr2);
      currTr->AddHit(vtx->GetActiveID(), chi2v, -1);
      currTr->SetInnerLrChecked(vtx->GetActiveID());
      if (NeedToKill(currTr)) currTr->Kill();
      // if (vtx->IsITS()) {if (!UpdateTrack(currTr, vtx, vtx->GetMCCluster(), kFALSE)) {currTr->Kill();continue;}}
    }
  }
  EliminateUnrelated();
  
  return kTRUE;
}



//____________________________________________________________________________
Int_t KMCDetector::GenBgClusters(KMCLayer* lr)
{
  // Generate fake clusters in precalculated RPhi,Z range
  if (fNBgLimits<1) return 0; // limits were not set - no track was prolongated
  //
  // Fix search limits to avoid seeds which will anyway point very far from the vertex
  double tolY = TMath::Sqrt(lr->fSig2EstD)*fMaxChi2ClSQ;
  double tolZ = TMath::Sqrt(lr->fSig2EstZ)*fMaxChi2ClSQ;
  
  //  printf("Before: Y: %+6.3f : %+6.3f tolY: %6.3f || Z: %+6.3f : %+6.3f tolZ: %6.3f\n",fBgYMin,fBgYMax,tolY, fBgZMin,fBgZMax,tolZ);
  if (fBgYMin < lr->fClCorr.fY-tolY) fBgYMin = lr->fClCorr.fY-tolY;
  if (fBgYMax > lr->fClCorr.fY+tolY) fBgYMax = lr->fClCorr.fY+tolY;
  if (fBgZMin < lr->fClCorr.fZ-tolZ) fBgZMin = lr->fClCorr.fZ-tolZ;
  if (fBgZMax > lr->fClCorr.fZ+tolZ) fBgZMax = lr->fClCorr.fZ+tolZ;
  //printf("After: Y: %+6.3f : %+6.3f tolY: %6.3f || Z: %+6.3f : %+6.3f tolZ: %6.3f\n",fBgYMin,fBgYMax,tolY, fBgZMin,fBgZMax,tolZ);
  //
  double dy = fBgYMax - fBgYMin;
  double dz = fBgZMax - fBgZMin;
  double surf = dy*dz;               // surface of generation
  if (surf<0) return 0;
  double poissProb = surf*HitDensity(lr->fR)*lr->GetLayerEff();
  AliDebug(2,Form("Bg for Lr %s (r=%.2f) : Density %.2f on surface %.2e [%+.4f : %+.4f][%+.4f %+.4f]",
		  lr->GetName(),lr->fR,HitDensity(lr->fR),surf,fBgYMin,fBgYMax,fBgZMin,fBgZMax));
  int nFakesGen = gRandom->Poisson( poissProb ); // preliminary number of extra clusters to test
  KMCCluster *refCl = lr->GetMCCluster();
  double sig2y = lr->GetPhiRes()*lr->GetPhiRes();
  double sig2z = lr->GetZRes()*lr->GetZRes();
  for (int ic=nFakesGen;ic--;) {
    double y = fBgYMin+dy*gRandom->Rndm();
    double z = fBgZMin+dz*gRandom->Rndm();
    double dfy = y-refCl->GetY();
    double dfz = z-refCl->GetZ();
    double dist = (dfy*dfy)/sig2y + (dfz*dfz)/sig2z;
    if (dist<4) continue; // avoid overlap with MC cluster
    lr->AddBgCluster(y, z, refCl->GetX(), refCl->GetPhi());
  }
  AliDebug(2,Form("Added %6d noise clusters on lr %s (poisson Prob=%8.2f for surface %.2e) DY:%7.4f DZ: %7.4f",
		   lr->GetNBgClusters(),lr->GetName(),poissProb,surf,dy,dz));
  return nFakesGen;
  //
}

//____________________________________________________________________________
void KMCDetector::UpdateSearchLimits(KMCProbe* probe, KMCLayer* lr)
{
  // define the search window for track on layer (where the bg hist will be generated)
  static double *currYMin = fBgYMinTr.GetArray();
  static double *currYMax = fBgYMaxTr.GetArray();
  static double *currZMin = fBgZMinTr.GetArray();
  static double *currZMax = fBgZMaxTr.GetArray();
  //
  double sizeY = probe->GetSigmaY2(), sizeZ = probe->GetSigmaZ2();
  //
  //  if (sizeY>2) sizeY=2;
  //  if (sizeZ>2) sizeZ=2;
  //  printf("Sizes at %s: %.5f %.5f\n",lr->GetName(), sizeY,sizeZ);
  //
  if (fNBgLimits>=fBgYMinTr.GetSize()) { // expand arrays, update pointers
    fBgYMinTr.Set(2*(fNBgLimits+1));
    fBgYMaxTr.Set(2*(fNBgLimits+1));
    fBgZMinTr.Set(2*(fNBgLimits+1));
    fBgZMaxTr.Set(2*(fNBgLimits+1));
    currYMin = fBgYMinTr.GetArray();
    currYMax = fBgYMaxTr.GetArray();
    currZMin = fBgZMinTr.GetArray();
    currZMax = fBgZMaxTr.GetArray();
  }
  if (fBgYMin > (currYMin[fNBgLimits]=probe->GetY()-sizeY) ) fBgYMin = currYMin[fNBgLimits];
  if (fBgYMax < (currYMax[fNBgLimits]=probe->GetY()+sizeY) ) fBgYMax = currYMax[fNBgLimits];
  if (fBgZMin > (currZMin[fNBgLimits]=probe->GetZ()-sizeZ) ) fBgZMin = currZMin[fNBgLimits];
  if (fBgZMax < (currZMax[fNBgLimits]=probe->GetZ()+sizeZ) ) fBgZMax = currZMax[fNBgLimits];
  if (AliLog::GetGlobalDebugLevel()>=2) {
    probe->Print("l");
    AliInfo(Form("Seed%3d Lr %s limits for y:%+8.4f z:%+8.4f [%+.4f : %+.4f][%+.4f %+.4f]",fNBgLimits,lr->GetName(),probe->GetY(),probe->GetZ(),currYMin[fNBgLimits],currYMax[fNBgLimits],currZMin[fNBgLimits],currZMax[fNBgLimits]));
    AliInfo(Form("Global Limits Lr %s                            [%+.4f : %+.4f][%+.4f %+.4f]",lr->GetName(),fBgYMin,fBgYMax,fBgZMin,fBgZMax));
    AliInfo(Form("MC Cluster: %+.4f : %+.4f",lr->fClMC.fY, lr->fClMC.fZ));
  }
  probe->SetUniqueID(fNBgLimits++);
  //
  if (lr->IsITS() && probe->GetNFakeITSHits()==0) {
    if (fHMCLrResidRPhi) fHMCLrResidRPhi->Fill(probe->GetY() - lr->GetMCCluster()->GetY(), lr->GetActiveID());
    if (fHMCLrResidZ)    fHMCLrResidZ->Fill(probe->GetZ() - lr->GetMCCluster()->GetZ(),lr->GetActiveID());
  }
  //
}

//____________________________________________________________________________
void KMCDetector::CheckTrackProlongations(KMCProbe *probe, KMCLayer* lr, KMCLayer* lrP)
{
  // explore prolongation of probe from lrP to lr with all possible clusters of lrP
  // the probe is already brought to clusters frame
  int nCl = lrP->GetNBgClusters();
  double measErr2[3] = {lrP->fPhiRes*lrP->fPhiRes,0,lrP->fZRes*lrP->fZRes};
  double meas[2] = {0,0};
  UInt_t tmpID = probe->GetUniqueID();
  double yMin = fBgYMinTr[tmpID];
  double yMax = fBgYMaxTr[tmpID];
  double zMin = fBgZMinTr[tmpID];
  double zMax = fBgZMaxTr[tmpID];
  //
  probe->SetInnerLrChecked(lrP->GetActiveID());
  for (int icl=-1;icl<nCl;icl++) {
    KMCCluster* cl = icl<0 ? lrP->GetMCCluster() : lrP->GetBgCluster(icl);  // -1 is for true MC cluster
    if (cl->IsKilled()) {
      if (AliLog::GetGlobalDebugLevel()>1) {printf("Skip cluster %d ",icl); cl->Print();}
      continue;
    }
    double y = cl->GetY();
    double z = cl->GetZ();
    AliDebug(2,Form("Check seed%d against cl#%d out of %d at layer %s | y:%+8.4f z:%+8.4f [%+.4f:%+.4f]  [%+.4f:%+.4f]",tmpID,icl,nCl,lrP->GetName(),y,z,yMin,yMax,zMin,zMax));
    if (AliLog::GetGlobalDebugLevel()>0) {
      if (icl==-1 && probe->GetNFakeITSHits()==0) {
	meas[0] = y; meas[1] = z;
	double chi2a = probe->GetPredictedChi2(meas,measErr2);
	if (chi2a>fMaxChi2Cl || (y<yMin || y>yMax) || (z<zMin || z>zMax)) {
	  probe->Print();
	  printf("Loosing good point (y:%+8.4f z:%+8.4f) on lr %s: chi2: %.2f  | dy:%+8.4f dz:%+8.4f [%+.4f:%+.4f]  [%+.4f:%+.4f] |x: %.2f %.2f | phi: %.2f %.2f\n",
		 y,z,lrP->GetName(),chi2a,y-probe->GetY(),z-probe->GetZ(),yMin,yMax,zMin,zMax, probe->GetX(), cl->GetX(), probe->GetAlpha(), cl->GetPhi());
	}
      }
    }
    if (y<yMin || y>yMax) continue; // preliminary check on Y
    if (z<zMin || z>zMax) continue; // preliminary check on Z
    meas[0] = y; meas[1] = z;
    double chi2 = probe->GetPredictedChi2(meas,measErr2);
    if (fHMCLrChi2 && probe->GetNFakeITSHits()==0 && icl==-1) fHMCLrChi2->Fill(chi2,lrP->GetActiveID());
    AliDebug(2,Form("Seed-to-cluster chi2 = Chi2=%.2f",chi2));
    if (chi2>fMaxChi2Cl) continue;
    // 
    // update track copy
    KMCProbe* newTr = lr->AddMCTrack( probe );
    fUpdCalls++;
    if (!newTr->Update(meas,measErr2)) {
      AliDebug(2,Form("Layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		      lrP->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
      if (AliLog::GetGlobalDebugLevel()>1) newTr->Print("l");
      newTr->Kill();
      continue;
    }
    newTr->AddHit(lrP->GetActiveID(), chi2, icl);
    if (AliLog::GetGlobalDebugLevel()>1) {
      AliInfo("Cloned updated track is:");
      newTr->Print();
    }
    if (NeedToKill(newTr)) newTr->Kill();
  }
  //
}

//____________________________________________________________________________
Bool_t KMCDetector::NeedToKill(KMCProbe* probe) const
{
  // check if the seed at given layer (last one where update was tried) 
  // still has chances to be reconstructed
  const Bool_t kModeKillMiss = kFALSE;
  //
  Bool_t kill = kFALSE;
  while (1) {
    int il = probe->GetInnerLayerChecked();
    int nITS = probe->GetNITSHits();
    int nITSMax = nITS + il; // maximum it can have
    if (nITSMax<fMinITSHits) {
      kill = kTRUE; 
      break;
    }    // has no chance to collect enough ITS hits
    //
    int ngr = fPattITS.GetSize();
    if (ngr>0) { // check pattern
      UInt_t patt = probe->GetHitsPatt();
      // complete the layers not checked yet
      for (int i=il;i--;) patt |= (0x1<<i);
      for (int ig=ngr;ig--;) 
	if (!(((UInt_t)fPattITS[ig]) & patt)) {
	  kill = kTRUE; 
	  break;
	}
      //
    }
    //
    if (nITS>2) {  // check if smallest possible norm chi2/ndf is acceptable
      double chi2min = probe->GetChi2();
      if (kModeKillMiss) {
	int nMiss = fNActiveITSLayers - probe->GetInnerLayerChecked() - nITS; // layers already missed
	chi2min = nMiss*probe->GetMissingHitPenalty();
      }
      chi2min /= ((nITSMax<<1)-KMCProbe::kNDOF);
      if (chi2min>fMaxNormChi2NDF) {
	kill = kTRUE; 
	break;
      }
    }
    //
    // loose vertex constraint
    double dst;
    if (nITS>=2) {
      probe->GetZAt(0,fBFieldG,dst);
      //printf("Zd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    if (nITS>=3) {
      probe->GetYAt(0,fBFieldG,dst);
      //printf("Dd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    //
    break;
  }
  if (kill && AliLog::GetGlobalDebugLevel()>1 && probe->GetNFakeITSHits()==0) {
    printf("Killing good seed, last upd layer was %d\n",probe->GetInnerLayerChecked());
    probe->Print("l");
  }
  return kill;
}

//_____________________________________________________________________
void KMCDetector::EliminateUnrelated()
{
  // kill useless tracks
  KMCLayer* lr = GetLayer(0);
  int ntr = lr->GetNMCTracks();
  int nval = 0;
  for (int itr=0;itr<ntr;itr++) {
    KMCProbe* probe = lr->GetMCTrack(itr);
    if (probe->IsKilled()) continue;
    if (probe->GetNITSHits()-probe->GetNFakeITSHits()<1) {probe->Kill(); continue;}
    nval++;
  }
  lr->GetMCTracks()->Sort();
  const int kDump = 0;
  if (kDump>0) {
    printf("Valid %d out of %d\n",nval, ntr);
    ntr = ntr>kDump ? kDump:0;
    for (int itr=0;itr<ntr;itr++) {
      lr->GetMCTrack(itr)->Print();
    }
  }
}

//_____________________________________________________________________
void KMCDetector::RequirePattern(UInt_t *patt, int groups)
{
  if (groups<1) {fPattITS.Set(0); return;}
  fPattITS.Set(groups);
  for (int i=0;i<groups;i++) fPattITS[i] = patt[i];
}

//_____________________________________________________________________
void KMCDetector::CalcHardSearchLimits(double dzv)
{
  //
  TArrayD zlims;
  zlims.Set(fNActiveITSLayers);
  for (int il0=0;il0<fNActiveITSLayers;il0++) {
    KMCLayer* lr0 = GetActiveLayer(il0);
    double angZTol = dzv/lr0->GetRadius();
    for (int il1=0;il1<fNActiveITSLayers;il1++) {
      if (il1==il0) continue;
      KMCLayer* lr1 = GetActiveLayer(il1);
      double ztol = angZTol*TMath::Abs(lr0->GetRadius() - lr1->GetRadius());
      if (ztol>zlims[il1]) zlims[il1] = ztol;
    }
  }
  //
  for (int il=0;il<fNActiveITSLayers;il++) printf("ZTol%d: %8.4f\n",il,zlims[il]);
}

//_______________________________________________________
double KMCDetector::PropagateBack(KMCProbe* trc) 
{
  static KMCProbe bwd;
  bwd = *trc;
  bwd.ResetCovMat();
  static double measErr2[3] = {0,0,0};
  static double meas[2] = {0,0};
  int icl = 0;
  double chi2Tot = 0;
  for (int il=1;il<=fLastActiveITSLayer;il++) {
    KMCLayer* lr = GetLayer(il);
    if (!PropagateToLayer(&bwd,lr,1)) return -1;
    int aID = lr->GetActiveID();
    if (aID>-1 && (icl=bwd.fClID[aID])>=-1) {
      KMCCluster* clMC =  icl<0 ? lr->GetMCCluster() : lr->GetBgCluster(icl);
      if (!bwd.PropagateToCluster(clMC,fBFieldG)) return -1;
      meas[0] = clMC->GetY(); meas[1] = clMC->GetZ();
      measErr2[0] = lr->fPhiRes*lr->fPhiRes;
      measErr2[2] = lr->fZRes*lr->fZRes;
      double chi2a = bwd.GetPredictedChi2(meas,measErr2);
      chi2Tot += chi2a;
      printf("Chis %d (cl%+3d): t2c: %6.3f tot: %6.3f\n",aID,icl,chi2a, chi2Tot);
      bwd.Update(meas,measErr2);
      bwd.AddHit(aID, chi2a, icl);	
    }
    if (!bwd.CorrectForMeanMaterial(lr,kFALSE)) return -1;
  }
  return chi2Tot;
}
*/
