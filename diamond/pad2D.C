{
gSystem->Load("../lib/KDetSim.sl");
gStyle->SetCanvasPreferGL(kTRUE);
gStyle->SetOptStat(0);
KDetector  det;
det.SetDriftHisto(15e-9); //timerange=15e-9 (numbins=200 by default)

// Detector binning
// Inherited from Class KGeometry
det.nx = 450;
det.ny = 100;
det.nz = 1;

// Detector dimensions (Diamond: Heisenberg)
Float_t dimX = 4500; //micro
Float_t dimY = 510; //micro
Float_t dimZ = 1; //thickness: 1 micro

// Pad sizes
Float_t pDimX = 2000; // Side Length = 2.0*pDimX
Float_t pDimY = 0.1; //thickness: 0.1 micro
Float_t pDimZ = 0.5;

//Pad position
Float_t pXpos = dimX/2.0;
Float_t pZpos = dimZ/2.0;
Float_t pYlow = 0; //lower Pad
Float_t pYhigh = dimY; //Upper Pad

// Material constants
Int_t diamond = 10;
Int_t aluminum = 100; //actually we use Cr/Au but that should not matter

// Space charge
// Float_t space_charge = -0.25; //equivalent to a 2.5e11 space charge
Float_t space_charge = 0;

// Bias voltages
det.Voltage = 400;

// Electrode setup
Int_t gndBit = 2; //defines ground electrode
Int_t padBit = 16385; // bit2 = 1 (HV electrode), bit14 = 1 (electrode for which Ramo field is calculated)

// MIP properties
// the injection point of edge-TCT laser
Float_t entryPointX = 2000;
Float_t entryPointY = 100;
Float_t entryPointZ = 0;

//Float_t extX = 40;
Float_t extX = 0;
Float_t extY = 0;
Float_t extZ = 1;

// Diffusion
det.diff=1;

// Output directory
const char *currentProfile = "../Results/current_profile.dat";

// Other
Bool_t drawing = true;

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
	det.NeffH->SetBinContent(i, j, k, space_charge); // wait for change
	}

//  GND Pad
Float_t GndPos[3]={pXpos, pYlow, pZpos}; 
Float_t GndSiz[3]={pDimX, pDimY, pDimZ};
det.ElRectangle(GndPos, GndSiz, gndBit, aluminum);

//  Top Pad
Float_t PadPos[3]={pXpos, pYhigh, pZpos};
Float_t PadSiz[3]={pDimX, pDimY, pDimZ};
det.ElRectangle(PadPos, PadSiz, padBit, aluminum);

det.SetBoundaryConditions();
det.CalField(0); //calculate electric field
det.CalField(1); //calculate weighting field

// MIP Properties
det.SetEntryPoint(entryPointX,entryPointY,entryPointZ);
// 2000 100 0
det.SetExitPoint(entryPointX+extX,entryPointY+extY,entryPointZ+extZ);
// 2000 100 1

det.MipIR(100,0);

// Drawing
if(drawing){

	// Draw Electric Field 2D Map
	TCanvas *c1 = new TCanvas("c1","c1",1200,400);
	c1->cd();
	TH2 *hEFxy = det->Draw("EFxy");
	hEFxy->GetYaxis()->SetRangeUser(20,500);
	//hEFxy->GetYaxis()->SetRangeUser(0,510);
	//hEFxy->GetXaxis()->SetRangeUser(20,4480);
	hEFxy->GetZaxis()->SetTitle("Electric Field [V/#mum]");
	hEFxy->SetTitle("4.5x4.5mm Diamond Pad Detector");
	c1->SetRightMargin(0.15);
	hEFxy->Draw("COLZ");

	// Draw Electric Field Profile on Y
	TCanvas *c2 = new TCanvas("c2","c2");
	c2->cd();
	TH1 *projy = hEFxy->ProjectionY("",50, 50);
	projy.Draw();
    
	/*
	// Draw Electric Field Profile on X
	TCanvas *c3 = new TCanvas("c3", "c3");
	c3->cd();
	TH1 *projx = hEFxy->ProjectionX("",50, 50);
	projx.Draw();
	*/

	TCanvas *c3 = new TCanvas("c3","c3", 1200,400);
	c3->cd();
	det.ShowMipIR(100);
	//det.

	TCanvas *c4 = new TCanvas("c4","c4");
	c4->cd();
	TH1 *charge = det.sum; //draw current

	charge->GetXaxis()->SetRangeUser(-5e-9,10e-9);
	charge->Draw();

	//draw hole
	TH1 *hole = det.pos;
	hole->Draw("SAME");

	//draw electron
	TH1 *elec = det.neg;
	elec->Draw("SAME");

}



//write current profile to file for later analysis
std::ofstream ofs (currentProfile, std::ofstream::out);
ofs << "MIP infos" << std::endl;


TAxis *axis = charge->GetXaxis();
for(bin=0; bin<charge->GetNbinsX(); bin++){
	Double_t binCenter = axis->GetBinCenter(bin);
	Double_t binContent = charge->GetBinContent(bin);
	ofs << binCenter << ' ' << binContent << std::endl;
}
ofs.close();

}

