{
gSystem->Load("../lib/KDetSim.sl");
gStyle->SetCanvasPreferGL(kTRUE);

KDetector  det;

 det.SetDriftHisto(15e-9);
// Detector binning
det.nx = 100;
det.ny = 250;
det.nz = 25;

// Detector dimensions
Float_t dimX = 400;
Float_t dimY = 400;
Float_t dimZ = 100;

// Pad sizes
Float_t pDimX = 180;
Float_t pDimY = 180;
Float_t pDimZ = 1;

//Pad position
Float_t pXpos = dimX/2.0;
Float_t pYpos = dimY/2.0;
Float_t pZlow = 0.5;
Float_t pZhigh = dimZ-0.1;

// Material constants
Int_t diamond = 10;
Int_t aluminum = 100; //actually we use Cr/Au but that should not matter

// Space charge
Float_t space_charge = 0.1; //equivalent to a 1e11 space charge

// Bias voltages
det.Voltage = -400;

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

//  GND Pad
Float_t GndPos[3]={pXpos, pYpos, pZlow}; 
Float_t GndSiz[3]={pDimX, pDimY, pDimZ}; 
det.ElRectangle(GndPos, GndSiz, gndBit, aluminum);


//  Top Pad
Float_t PadPos[3]={pXpos, pYpos, pZhigh};
Float_t PadSiz[3]={pDimX, pDimY, pDimZ};
det.ElRectangle(PadPos, PadSiz, padBit, aluminum);

det.SetBoundaryConditions();
det.CalField(0);
det.CalField(1);
 
//define entry and exit point of mip
det.enp[0]=200;
det.enp[1]=200;
det.enp[2]=1;
det.exp[0]=100;
det.exp[1]=100;
det.exp[2]=50;

TCanvas c1; c1.cd();
det.ShowMipIR(20); //visualize track
TCanvas c2; c2.cd();

det.MipIR(200); //calculate drift 
det.sum->Draw(); //draw current
det.pos->Draw("SAME"); // hole current
det.neg->Draw("SAME");  // electron current

TCanvas c3; c3.cd();
det->Draw("EFxz")->Draw("COLZ");
}

