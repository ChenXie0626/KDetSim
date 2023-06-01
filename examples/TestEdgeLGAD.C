{ 
 gStyle->SetCanvasPreferGL(kTRUE); 
 gStyle->SetOptStat(0); 
 KDetector  det; 
 det.SetDriftHisto(5e-9); 
 
 
 // Detector binning 
 det.nx = 100; 
 det.ny = 50; 
 det.nz = 1; 
 
 
 // Detector dimensions 
 Float_t dimX = 100; 
 Float_t dimY = 50; 
 Float_t dimZ = 1; 
 
 
 // Pad sizes 
 Float_t pDimX = 100; 
 Float_t pDimY = 1; 
 Float_t pDimZ = 1; 
 
 
 //Pad position 
 Float_t pXpos = dimX/1.0; 
 Float_t pZpos = 1; 
 Float_t pYlow = 0; 
 Float_t pYhigh = dimY; 
 
 
// Material constants 
 Int_t diamond = 0; 
 Int_t aluminum = 100; //actually we use Cr/Au but that should not matter 
 
 
// Space charge 
Float_t space_charge = -18; //equivalent to a 2.5e11 space charge 
 
 
// Bias voltages 
 det.Voltage = -130; 
 
 
// Electrode setup 
 Int_t gndBit = 2; //defines ground electrode 
 Int_t padBit = 16385; // bit2 = 1 (HV electrode), bit14 = 1 (electrode for which Ramo field is calculated) 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 
 
//  Sensor geometry 
 det.EG = new TH3I("EG","EG",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ);  
 det.EG->GetXaxis()->SetTitle("x [#mum]");  
 det.EG->GetYaxis()->SetTitle("y [#mum]");  
 det.EG->GetZaxis()->SetTitle("z [#mum]"); 
 
 
 //  Detector material 
 det.DM = new TH3I("DM","DM",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ); 
 det.DM->GetXaxis()->SetTitle("x [#mum]");  
 det.DM->GetYaxis()->SetTitle("y [#mum]");  
 det.DM->GetZaxis()->SetTitle("z [#mum]"); 
 
 
 //  Space charge 
 det.NeffH = new TH3F("Neff","Neff",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ); 
 det.NeffH->GetXaxis()->SetTitle("x [#mum]"); 
 det.NeffH->GetYaxis()->SetTitle("y [#mum]"); 
 det.NeffH->GetZaxis()->SetTitle("z [#mum]"); 
 
 
//  Material and space charge setup: 
 Int_t i, j, k; 
 for(int k=1; k<=det.nz; k++)  
   for(int j=1; j<=det.ny; j++) 
     for(int i=1; i<=det.nx; i++)  
       { 
 	det.DM->SetBinContent(i, j, k, diamond); 
 	det.NeffH->SetBinContent(i, j, k, space_charge); 
       } 
 
 
 //  Electrode
 Float_t ElPos[3]={28, 50,0.5};  
 Float_t ElSiz[3]={28, 1, 0.5};  
 det.ElRectangle(ElPos, ElSiz, padBit , 100); 
 
 
 
 //  Top Pad 
 Float_t PadPos[3]={50, 1, 0.5}; 
 Float_t PadSiz[3]={50, 1, 0.5}; 
 det.ElRectangle(PadPos, PadSiz, 2, aluminum); 
 
 

 
 

 
 for(int j=50; j>=47; j--)
     for(int i=1; i<=50; i++)
       {
	        det.NeffH->SetBinContent(i, j, 1, -10000);
       }

 for(int j=50; j>=45; j--)
     for(int i=51; i<=56; i++)
       {
		det.NeffH->SetBinContent(i, j, 1, 1000);
       }

  det.SetBoundaryConditions(); 
  det.CalField(0); 
  det.CalField(1); 
det.SetExitPoint(100,5,0.5);
 det.SetEntryPoint(1,5,0.5);
 det.ShowMipIR(50);
}
/* } */
/*  TCanvas *c2 = new TCanvas("c2","c2");  */
/*  c2->cd();  */
/*  TH1 *projy = hEFxy->ProjectionY("",50, 50);  */
/*  projy.Draw();  */
 
 
/* entryPointX = 300;  */
/* entryPointY = 100;  */
/* entryPointZ = 0;  */

 
/* det.diff=1;  */
/* det.SetEntryPoint(entryPointX,entryPointY,entryPointZ);  */
/* det.SetExitPoint(entryPointX+20,entryPointY,entryPointZ+1);  */
/* det.MipIR(100);  */

 
/* TCanvas *c3 = new TCanvas("c3","c3", 1200,400);  */
/*  c3->cd();  */
/* det.ShowMipIR(50);  */

 
/* TCanvas *c4 = new TCanvas("c4","c4");  */
/* c4->cd();  */
/* TH1 *charge = det.sum; //draw current  */

 
/* charge->GetXaxis()->SetRangeUser(-5e-9,10e-9);  */
/* charge->Draw();  */

 
} 
