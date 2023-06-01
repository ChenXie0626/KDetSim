{
gStyle->SetCanvasPreferGL(kTRUE);
gStyle->SetOptStat(0);
KDetector  det;
det.SetDriftHisto(15e-9);

// Detector binning
Int_t nx = 460;
Int_t ny = 50;
Int_t nz = 1;

// Detector dimensions in micrometer
Int_t dimX = 4600;
Int_t dimY = 500;
Int_t dimZ = 1;

// Pad sizes
Float_t pDimX = 2000;
Float_t pDimY = 1;
Float_t pDimZ = 1;

//Pad position
Float_t pXpos = dimX/2.0;
Float_t pZpos = 0;
Float_t pYlow = 0;
Float_t pYhigh = dimY;

// Material constants
Int_t diamond = 10;
Int_t aluminum = 100; //=Al, actually we use Cr/Au but that should not matter

// Space charge
Float_t space_charge = 1; //equivalent to a number * 1e12 1/ccm space charge

// Bias voltages
det.Voltage = -400;

// Electrode setup
Int_t gndBit = 1; //defines ground electrode
Int_t padBit = 16386; // bit2 = 1 (HV electrode), bit14 = 1 (electrode for which Ramo field is calculated)

// MIP properties
Int_t entryPointX = 2000;
Int_t entryPointY = 50;
Int_t entryPointZ = 0;
Int_t extX = 0;
Int_t extY = 0;
Int_t extZ = 1;

// Diffusion
det.diff=1;

// Output directory
const char *currentProfile="../Results/prof_y=50.dat";

// Other
Bool_t drawing = true;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//  Sensor geometry
det.EG = new TH3I("EG","EG", nx,0,dimX, ny,0,dimY, nz,0,dimZ); 
det.EG->GetXaxis()->SetTitle("x [#mum]"); 
det.EG->GetYaxis()->SetTitle("y [#mum]"); 
det.EG->GetZaxis()->SetTitle("z [#mum]");

//  Detector material
det.DM = new TH3I("DM","DM", nx,0,dimX, ny,0,dimY, nz,0,dimZ); 
det.DM->GetXaxis()->SetTitle("x [#mum]"); 
det.DM->GetYaxis()->SetTitle("y [#mum]"); 
det.DM->GetZaxis()->SetTitle("z [#mum]");

//  Space charge
det.NeffH = new TH3F("Neff","Neff", nx,0,dimX, ny,0,dimY, nz,0,dimZ); 
det.NeffH->GetXaxis()->SetTitle("x [#mum]");
det.NeffH->GetYaxis()->SetTitle("y [#mum]");
det.NeffH->GetZaxis()->SetTitle("z [#mum]");



std::cout << "pad bin location: " << det.EG->GetXaxis()->FindBin(pXpos) << std::endl;

//  Material and space charge setup:
Int_t i, j, k;
Float_t neff = 5;
for(int k=1; k<=nz; k++) 
  for(int j=1; j<=ny; j++)
    for(int i=1; i<=nx; i++){
		det.DM->SetBinContent(i, j, k, diamond);

		//linear
		//det.NeffH->SetBinContent(i, j, k, 0.6*(-1.0*neff+j/neff));
		
		//quadratic
		det.NeffH->SetBinContent(i, j, k, 10.0/15625*pow((j-25),3));
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
det.CalField(0);
det.CalField(1);

// MIP Properties
det.SetEntryPoint(entryPointX,entryPointY,entryPointZ);
det.SetExitPoint(entryPointX+extX,entryPointY+extY,entryPointZ+extZ);
det.MipIR(100,0);

// Drawing
if(drawing){
	TCanvas *c1 = new TCanvas("c1","c1",1200,400);
	c1->cd();
	TH2 *hEFxy = det->Draw("EFxy");
	//hEFxy->GetYaxis()->SetRangeUser(10,480);
	//hEFxy->GetXaxis()->SetRangeUser(20,4580);
	hEFxy->GetZaxis()->SetTitle("Electric Field [V/#mum]");
	hEFxy->SetTitle("Electric Field");
	c1->SetRightMargin(0.15);
	hEFxy->Draw("COLZ");

	TCanvas *c2 = new TCanvas("c2","c2",1200,400);
	c2->cd();
	TH2 *hPxy = det->Draw("P");
	hPxy->GetZaxis()->SetTitle("Potential [V]");
	hPxy->SetTitle("Potential");
	c2->SetRightMargin(0.15);
	hPxy->Draw("COLZ");

	TCanvas *c3 = new TCanvas("c3","c3");
	c3->cd();
	TH1 *projy = hEFxy->ProjectionY("",50, 50);
	projy.Draw();

	//TCanvas *c4 = new TCanvas("c4","c4", 1200,400);
	//c4->cd();
	//det.ShowMipIR(100);


	//all gaussian beam stuff:
	TCanvas *c4 = new TCanvas("c4","c4", 1200,400);
	c4->cd();
	//det.GaussBeam(100, 0, 20, 2, 0, 0);
	//det.ShowGaussBeam(100, 0, 20, 2, 1,0);



	TCanvas *c5 = new TCanvas("c5","c5");
	c5->cd();
	TH1 *charge = det.sum; //draw current
	TH1 *el = det.neg;
	TH1 *holes = det.pos;

	charge->GetXaxis()->SetRangeUser(-5e-9,10e-9);
	charge->Draw();
	el->Draw("SAME"); 
	holes->Draw("SAME");
	charge->Draw("SAME");



}

//write current profile to file for later analysis
std::ofstream ofs (currentProfile, std::ofstream::out);
TAxis *axis = charge->GetXaxis();
for(bin=0; bin<charge->GetNbinsX(); bin++){
	Double_t binCenter = axis->GetBinCenter(bin);
	Double_t binContent = charge->GetBinContent(bin);
	ofs << binCenter << ' ' << binContent << std::endl;
}
ofs.close();

}

