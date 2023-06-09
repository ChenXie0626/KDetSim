{
    TF1 *neff = new TF1("neff", "[0]+x[0]*0", 0, 1000);
    neff->SetParameter(0, 1);

    KPad det(50, 300);
    det.Neff = neff;
    det.Voltage = -200;
    det.SetUpVolume(1);
    det.SetUpElectrodes();
    //    det.MaxDriftLen=50; // uncomment if you want to restrict current calculation to first 50 microns of drift
    

    TCanvas c1;
    det.SetEntryPoint(25, 299.9, 0.5);
    det.SetExitPoint(25, 1., 0.5);
    det.Temperature = 253;
    det.diff = 1;
    det.ShowMipIR(200);

    TCanvas c2;
    det.MipIR(200);
    det.sum->DrawCopy();
    det.pos->Draw("SAME");
    det.neg->Draw("SAME");

    TCanvas c3;
    KElec tct;
    tct.RCshape(40e-12,50,10000,det.sum);
    det.sum->Draw("HIST");
    
}
