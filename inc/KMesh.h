// @(#)root/html:$Id: KMesh.h 27910 2012-10-22 17:26:55Z Krambi $
// Author: Gregor Kramberger   18/10/12

#ifndef _KMesh
#define _KMesh

#include "TMath.h"
#include "Rtypes.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KMesh                                                               //
//                                                                      //
// Calculation mesh generator                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef MSVC
class  __declspec(dllexport) KMesh
#else
class  KMesh
#endif
{

private:
public:
    Int_t N;
    Float_t Max;
    Float_t Min;
    //_______________________________________________________________________________
    KMesh(Float_t x0, Float_t x1 = 0)
    {
        Max = x0;
        Min = x1;
        N = 0;
    };
    ~KMesh(){};
    Int_t GetBins(Int_t Num, Float_t SS, Float_t ES, Float_t *X);
    Int_t GetBins(Int_t size, Float_t *Pos, Float_t *Step, Float_t *Bins);
    //  Double_t fdv() { return((Double_t) Neff1*1e-6*e_0*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
    //  Double_t fdv(Double_t Neff) {return((Double_t) Neff*1e-6*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
    ClassDef(KMesh,1) 
};

#endif

#ifndef _KElec
#define _KElec

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Electronics Class                                                    //
//                                                                      //
// Electronic hadling of current source                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//#include <iostream.h>
#ifdef MSVC
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

#include <stdio.h>
#include <stdlib.h>
#include "TObject.h"
#include "TH1.h"
#include "TMath.h"
#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "TGraph.h"
#include "TF1.h"



class EXPORT  KElec {

private:
 Int_t Method;
 Double_t Cp;
 Double_t Rp;
 Double_t Crc;
 Double_t R1rc;
 Double_t R2rc;
 Double_t Ccr;
 Double_t R1cr;
 Double_t R2cr;
 Double_t PeakTime;
 Double_t IntTime;
public:
  KElec(Double_t=5e-12,Double_t=50,Double_t=25e-9,Double_t=1e-12,Double_t=200,Double_t=200,Double_t=1e-12,Double_t=200,Double_t=200,Double_t=25e-9,Int_t=0);
  virtual ~KElec();
  Double_t Trapez(TH1F *,Int_t,Double_t);
  Double_t  Simpson(TH1F *,Int_t,Double_t);
  void preamp(TH1F *his) {preamp(Cp,Rp,his,IntTime,Method);};
  void preamp(Double_t, Double_t,TH1F *, Double_t=-1111, Int_t = 0);
  void Revpreamp(Double_t, Double_t R,TH1F *,Double_t=1);
  void Revpreamp(TH1F *his,Double_t unit=1) {Revpreamp(Cp,Rp,his,unit);};
  void RCshape(Double_t, Double_t, Double_t,TH1F *, Int_t = 0);
  void RCshape(TH1F *his) {RCshape(Crc,R1rc,R2rc,his,Method);};
  void CRshape(Double_t, Double_t, Double_t,TH1F *, Int_t = 0);
  void CRshape(TH1F *his) {CRshape(Ccr,R1cr,R2cr,his,Method);};
  void SetCp(Double_t x) {Cp=x;};  
  void SetCrc(Double_t x) {Crc=x;};
  void SetR1rc(Double_t x) {R1rc=x;};  
  void SetR2rc(Double_t x) {R2rc=x;};
  void SetRp(Double_t x) {Rp=x;};  
  void SetCcr(Double_t x) {Ccr=x;};
  void SetR1cr(Double_t x) {R1cr=x;};  
  void SetR2cr(Double_t x) {R2cr=x;};
  void SetIntTime(Double_t x) {IntTime=x;};
  Double_t GetPeakTime() {return(PeakTime);};
  void PrintPar() {printf("KElec Parameters: Cp=%e, Rp=%e, IntTime=%e Crc=%e, R1rc=%e, R2rc=%e,Ccr=%e, R1cr=%e, R2cr=%e, Method=%d\n",Cp,Rp,IntTime,Crc,R1rc,R2rc,Ccr,R1cr,R2cr,Method);}
  void SetMethod(Int_t x) {Method=x;};
  ClassDef(KElec,1) 
};
#endif

