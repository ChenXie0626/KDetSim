#include "KMesh.h"
#include "TMath.h"

ClassImp(KMesh)
    //////////////////////////////////////////////////////////////////////////
    //                                                                      //
    // KMesh                                                                //
    // A mesh generator for the non-equividistant bins                      //
    //                                                                      //
    //////////////////////////////////////////////////////////////////////////

    Int_t KMesh::GetBins(Int_t Num, Float_t SS, Float_t ES, Float_t *Bins)
{

    Float_t dX = (Max - Min) / Num;
    Float_t dS = (SS - ES) / (Num - 1);
    Float_t dSi, SumS = 0;
    Int_t Ni = 0, k = 0;
    Bins[0] = 0;
    N = 0;
    printf("Steps=%d, SumS=%f \n", k, SumS);
    k++;
    for (Int_t i = 0; i < Num; i++) {
        dSi = SS - i * dS;
        Ni = TMath::Nint(dX / dSi);
        printf("New region : %d %d %f\n", i, Ni, dSi);
        for (Int_t j = 0; j < Ni; j++) {
            SumS += dSi;
            Bins[k] = SumS;
            if (Bins[k] > Max) {
                Bins[k] = Max;
                printf("END\n");
                break;
            }
            printf("Steps=%d, SumS=%f \n", k, Bins[k]);
            k++;
        }
    }
    if (Bins[k - 1] < Max) {
        Bins[k] = Max;
        printf("Steps=%d, SumS=%f \n", k, Bins[k]);
    };

    return k;
}

Int_t KMesh::GetBins(Int_t size, Float_t *Pos, Float_t *Step, Float_t *Bins)
{
    Float_t N[100];
    Int_t num, NN = 0, i, j, k;
    for (i = 0; i < size; i++) {
        //      if(i>0) N[i]=(Pos[i]-Pos[i-1])/Step[i]; else
        N[i] = Pos[i] / Step[i];
        if (N[i] != TMath::Nint(N[i])) {
            printf("Step size not integer ...\n");
            return -1;
        } else
            NN += (Int_t)N[i];
    }
    printf("%d \n", NN);
    //  Bins=new Float_t [NN+1];

    Bins[0] = 0;
    k = 1;
    for (i = 0; i < size; i++) {
        for (j = 0; j < N[i]; j++) {

            Bins[k] = Bins[k - 1] + Step[i];
            printf("%d %f\n", k, Bins[k]);
            k++;
        }
    }

    return NN;
}


ClassImp(KElec)

KElec::KElec(Double_t x1,Double_t x2 ,Double_t x3,Double_t x4,Double_t x5,Double_t x6,Double_t x7,Double_t x8,Double_t x9,Double_t x10,Int_t met)
{
Cp=x1; Rp=x2; IntTime=x3; Crc=x4; R1rc=x5; R2rc=x6; Ccr=x7; R1cr=x8; R2cr=x9; PeakTime=x10;  Method=met;
}
KElec::~KElec()
{
  //Clear();
}

Double_t KElec::Trapez(TH1F *histo,Int_t i,Double_t tau)
{
Double_t f1=0,f2=0,h=0,t1=0,t2=0;
 if(i==1)  return(histo->GetBinContent(i)*histo->GetBinWidth(i)/2); else
   {     
 t1=histo->GetBinCenter(i-1); f1=histo->GetBinContent(i-1)*TMath::Exp(t1/tau); 
 t2=histo->GetBinCenter(i);   f2=histo->GetBinContent(i)*TMath::Exp(t2/tau);	
 h=histo->GetBinWidth(i); //printf("tutut %d,t1=%e,t2=%e,f1=%e,f2=%e\n",i,t1,t2,f1,f2);
 return(h*0.5*(f1+f2)); 

   }
}

Double_t KElec::Simpson(TH1F *histo,Int_t i,Double_t tau)
{
Double_t f1=0,f2=0,f3=0,h=0,t1=0,t2=0,t3=0;
 if(i<3 || i%2==0) return(Trapez(histo,i,tau)); else {
 h=histo->GetBinWidth(i);
 t1=histo->GetBinCenter(i-2); f1=histo->GetBinContent(i-2)*TMath::Exp(t1/tau); 
 t2=histo->GetBinCenter(i-1); f2=histo->GetBinContent(i-1)*TMath::Exp(t2/tau);
 t3=histo->GetBinCenter(i);   f3=histo->GetBinContent(i)*TMath::Exp(t3/tau);
 return(h * (0.33333333333*f1 +1.333333333*f2 + 0.3333333333333*f3)-Trapez(histo,i-1,tau));
  }
}

void KElec::Revpreamp(Double_t C, Double_t R, TH1F *histo,Double_t unit)
{
Int_t i;
Double_t tau=R*C*unit;
Float_t prev,nextprev;;
Int_t start;
//start=histo->GetXaxis()->FindBin(0);
prev=histo->GetBinContent(1);
histo->SetBinContent(1,(Float_t) tau*((histo->GetBinContent(2)-histo->GetBinContent(1))/(histo->GetBinWidth(1)) + histo->GetBinContent(1)/tau));
start=2;
//printf("%e %e\n", histo->GetBinWidth(1), tau);
for(i=start;i<histo->GetNbinsX()-1;i++)
  // histo->SetBinContent(i,(Float_t) tau*((histo->GetBinContent(i+1)-histo->GetBinContent(i-1))/(histo->GetBinWidth(i)) + histo->GetBinContent(i)/tau));
{
 nextprev=histo->GetBinContent(i);
 histo->SetBinContent(i,(Float_t) tau*((histo->GetBinContent(i+1)-prev)/(2*histo->GetBinWidth(i)) + histo->GetBinContent(i)/tau));
 prev=nextprev;
}
  
}

void KElec::preamp(Double_t C, Double_t R, TH1F *histo,Double_t cut, Int_t method)
{
static double e_0    = 1.60217733e-19; 
Double_t suma=0;
Double_t tau=R*C;
Double_t t;
Int_t i;
TH1F *whisto=new TH1F();
histo->Copy(*whisto);

if(cut==-1111) cut=histo->GetNbinsX()*histo->GetBinWidth(1);

for(i=1;i<histo->GetNbinsX()-1;i++)
  {
t=whisto->GetBinCenter(i);
  if(t<=cut)
         { 
	   if(method==0) suma+=Trapez(whisto,i,tau);
	   if(method==1) suma+=Simpson(whisto,i,tau);
	 }
       	     histo->SetBinContent(i,(Float_t) (1/C*suma*TMath::Exp(-t/tau)));
	     // printf("%d,t=%e,suma=%e,tau=%e,tapez=%e\n",i,t,suma,tau,Trapez(histo,i,tau));
  }
delete whisto;
}


void KElec::RCshape(Double_t C, Double_t R1, Double_t R2,TH1F *histo, Int_t method)
{
Double_t suma=0;
Double_t tau=(R1*R2)*C/(R1+R2);
Double_t tau1=R1*C;
Double_t t;
Int_t i;
TH1F *whisto=new TH1F();
histo->Copy(*whisto);

for(i=1;i<histo->GetNbinsX()-1;i++)
  {
    t=whisto->GetBinCenter(i);	
    if(method==0)  suma+=Trapez(whisto,i,tau);
    if(method==1)  suma+=Simpson(whisto,i,tau);
    histo->SetBinContent(i,(Float_t) (1/tau1*suma*TMath::Exp(-t/tau)));
  }
delete whisto;
}

void KElec::CRshape(Double_t C, Double_t R1, Double_t R2,TH1F *histo, Int_t method)
{
Double_t suma=0;
Double_t tau=(R1*R2)*C/(R1+R2);
Double_t tau1=R1*C;
Double_t t,f;
Int_t i;
TH1F *whisto=new TH1F();
histo->Copy(*whisto);

for(i=1;i<histo->GetNbinsX()-1;i++)
  {
    t=whisto->GetBinCenter(i);
    f=whisto->GetBinContent(i);	
    if(method==0)  suma+=Trapez(whisto,i,tau);
    if(method==1)  suma+=Simpson(whisto,i,tau);
    histo->SetBinContent(i,f-(1/tau-1/tau1)*suma*TMath::Exp(-t/tau));
  }
delete whisto;
}


