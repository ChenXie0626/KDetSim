{
  // gROOT->ProcessLine(".x root5logon.C");
  Float_t Iw=0.18; //width to the Gaussion distribution of the gain layer [um]
  Float_t Ip=2.2;  //position of the implant inside the sensors in [um]
  Float_t Id=-20e3; //peak doping concentration in gain layer [1e12 cm-3]
  Float_t SubCh=-40; //doping of the detector bulk [1e12 cm-3]
  Float_t BoxG;
  
  KMaterial::ImpactIonization = 1;  // impact ionization model - in root5 only =0 works

  //////////////////////////////////////////////////////////////////////////
  ///   DEFINITION OF THE DOPING FUNCTION FOR THE WHOLE DETECTOR          //
  //////////////////////////////////////////////////////////////////////////
  TF1 *neffG=new TF1("neff","TMath::Gaus(x,[2],[3])*[1]+[0]",0,50);
  neffG->SetParameter(0,SubCh);
  neffG->SetParameter(1,Id);
  neffG->SetParameter(2,Ip);
  neffG->SetParameter(3,Iw);
  BoxG=neffG->Integral(Ip-4*Iw,Ip+4*Iw); //Get the space charge under Gaussian
  printf("Integral of space charge Gaus=%f [cm-2] - Required Minimum Implantation Dose\n",BoxG*1e-4);
  
  //////////////////////////////////////////////////////////////////////////
  ///   DEFINITION OF THE DETECTOR                                        //
  //////////////////////////////////////////////////////////////////////////
  KPad *detG=new KPad(50,50);  // pad detector of 50 um size (not important) and 50 thick (important)
  detG->SetDriftHisto(20e-9);  // limit the induced current to 20 ns
  detG->Neff=neffG;            // set the doping function (defined above)
  detG->Voltage=+300;          // set the voltage
  detG->SetUpVolume(0.1);      // set the mesing of the volume -  100 nm mesh (1D problems only!)
  detG->SetUpElectrodes();     // set up both electrodes
  detG->SStep=0.1;             // set the drift step of e-h pairs
  detG->MTresh=1.2;            // set up the treshold for multiplication
  detG->BDTresh=1.2;           // set up the treshold for the detector breakdown 
  detG->Temperature=243;       // set the operation temperature
  detG->diff=1;                // switch on diffusion
  detG->SetEntryPoint(25,0,0.5); //set entry point of the track
  detG->SetExitPoint(25,50,0.5); //set exit point of the track

  //////////////////////////////////////////////////////////////////////////
  //      DEFINE GRAPHS                                                   //
  TGraph *ElField;      // electric field
  TGraph *ElPotential;  // electric potential
  TGraph *DopingProf;   // doping profile

  //////////////////////////////////////////////////////////////////////////

  TCanvas c2("Plots","Plots",1400,1000);   //open canvas
  c2.Divide(2,3);                          //divide canvas


  c2.cd(1);                                // plot in the first pad
  ElField=detG->DrawPad("f");              // draw electric field
  ElField->SetTitle("Electric field");
  ElField->GetXaxis()->SetTitle("depth y[#mum] (front electrode is at y=0)");
  ElField->GetYaxis()->SetTitle("E [V/#mum]");

  c2.cd(2);                                // plot in the second pad
  ElPotential=detG->DrawPad("p");          // draw electric potential
  ElPotential->SetTitle("Electric Potential");
  ElPotential->GetXaxis()->SetTitle("depth y[#mum] (front electrode at y=0)");
  ElPotential->GetYaxis()->SetTitle("U [V]");

  c2.cd(3);                                 // plot in the third pad
  TF1 *neffGc;                                  
  neffG->SetRange(0,5);
  neffGc=neffG->DrawCopy();
  neffGc->SetTitle("Doping profile (gain layer)");
  neffGc->GetXaxis()->SetTitle("depth y[#mum] (front electrode at 0)");
  neffGc->GetYaxis()->SetTitle("N_{eff} [10^{12} cm^{-3}]");
  neffGc->GetYaxis()->SetTitleOffset(1.5);

  
  c2.cd(4);                                 // plot in the fourth pad
  neffG->SetRange(0,50);
  neffG->SetTitle("Doping profile (full detector)");
  neffG->GetXaxis()->SetTitle("depth y[#mum] (front electrode at 0)");
  neffG->GetYaxis()->SetTitle("N_{eff} [10^{12} cm^{-3}]");
  neffG->GetYaxis()->SetTitleOffset(1.5);
  neffG->Draw();

  c2.cd(5);                                 // plot in the fifth pad
  detG->ShowMipIR(20);                      // plot the drift paths of the mip track 
                                            // subdivided into 20 buckets
  
  c2.cd(6);                                  // plot in the fifth pad
  Int_t divs=50;                             // simulated the track with 50 buckets - each
                                             // having ~100 e-h pairs bucket/um 

        detG->MipIR(divs);                   // simulate signal
  	detG->sum->SetLineColor(1); detG->sum->Draw("HIST");      //plot total induced current 
	detG->neg->SetLineColor(4); detG->neg->Draw("SAME HIST"); //plot electrons' induced current 
	detG->pos->SetLineColor(2); detG->pos->Draw("SAME HIST"); //plot holes' induced current

	//plot legends
	TLatex hole(10e-9,detG->sum->GetMaximum()*0.7,"hole contribution");
	hole.SetTextColor(2); hole.Draw("SAME");
        TLatex elec(10e-9,detG->sum->GetMaximum()*0.6,"hole contribution");
	elec.SetTextColor(4); elec.Draw("SAME");
        TLatex eph(10e-9,detG->sum->GetMaximum()*0.5,"holes and electrons");
	eph.SetTextColor(1); eph.Draw("SAME");
	// plot results
  	printf("Q=%f (%d e-h buckets created):: G=%f\n",detG->sum->Integral(),divs,detG->sum->Integral()/divs);
 

 //////////////////////////////////////////////////////////////////////////
 

 
  
 
}
