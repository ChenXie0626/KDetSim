{
  //Impact ionization model 0
  KMaterial::ImpactIonization=0;

  // Doping level function
  TF1 *neff=new TF1("neff","x<[0]?[2]:([3]*(x-[0])/([1]-[0]))+[2]",0,51);
  neff->SetParameter(0,47); //border of the implant region linear increase of the doping profile
  neff->SetParameter(1,50); //end of border region implant
  neff->SetParameter(2,-100);  //set bulk doping 1e14 -> Vfd_bulk~20 V
  neff->SetParameter(3,-8000); //doping of the implant region

  TCanvas cdop; cdop.cd(); neff->Draw(); //draw doping profile

  TGraph *Ef[10]; // Array of electric fields
  Float_t vol[10],max[10],maxdrift[10];

  KPad det(50,50); //create a simple pad detectors
  det.SetDriftHisto(5e-9,200);
  det.MTresh=1.1; //treshold for impact ionization (gain in step not more that 5%) - electrons
  det.BDTresh=1.1; //break down treshold - hole multiplication treshold
  det.Temperature=293; //temperature of operation 
  det.diff=0; //diffusion must be off in multiplication simulations - it doesn't play a role due to short drift times


  //Calculation of the field and simulations of the drift
  det.SetEntryPoint(25,1.5,0.5);
  det.SetExitPoint(25,49.5,0.5);
  det.SStep=0.05;         //stepping size in calculation 25 nm  - crucial for precision 

  for(int i=0;i<10;i++) //loop over all the volates
    { 
     det.Voltage=-270+i*(-40); //set voltage
     vol[i]=TMath::Abs(det.Voltage); 
     det.Neff=neff; //set the doping profile
     det.SetUpVolume(0.05);  //this parameter should be 0.05 to get to the right accuracy - 50 nm mesh size for simulation
     det.SetUpElectrodes();  //set up detector

     Ef[i]=det.DrawPad("f");   //draw electric field
     Ef[i]->SetLineColor(i+1);  //set line color of the field lines

     det.MipIR(10);  // calulation from drift
     maxdrift[i]=-det.sum->Integral()*400; //conversion into electrons
     printf("Max. multiplication = %f (drift simulation)\n",maxdrift[i]); //plot 
     
    }
  
  Ef[0]->Draw("AL"); 

  //calculation of the gain from ionization integral for electron injection (at the p+p contact) 
  for(int i=0;i<10;i++) 
    {
     Ef[i]->Draw("L");  //Plot lines 
     TH1D *EfH=GraphToHist(Ef[i]); //convert to histogram so that ionization integral can be calculated
     max[i]=KM(EfH,1)*3450; //Calculation of ionization integral
     printf("Max. multiplication = %f (ionization integral)\n",max[i]);
    }

// CALCULATION OF THE TRACK MULTIPLICATION
/*  double tot=0; */
/*      for(int i=1;i<=EfH->GetNbinsX();i++) */
/*        { */
/* 	 tot+=KM(EfH,EfH->GetBinCenter(i)); */
/*        } */
/*      tot/=EfH->GetNbinsX(); */

  TGraph *MGD=new TGraph(10,vol,max); 
  TGraph *MGC=new TGraph(10,vol,maxdrift);
  SetStyle(MGD,1,1,1,20);
  SetStyle(MGC,1,1,2,21);
  TCanvas c2;

  MGD->Draw("ALP");
  MGC->Draw("LP");
  for(int i=0;i<10;i++) printf("%f %f %f\n",vol[i],max[i],maxdrift[i]); 
} 
