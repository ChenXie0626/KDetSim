{
    TF1 *neff = new TF1("neff", "[0]+x[0]*0", 0, 1000);
    neff->SetParameter(0, 1);

    KPad det(50, 300);
    det.Neff = neff;
    det.SetDebug(0);
    det.Voltage = -200;
    det.SetUpVolume(1);
    det.SetUpElectrodes();
    //    det.MaxDriftLen=50; // uncomment if you want to restrict current calculation to first 50 microns of drift
    
    TF3 *TauE=new TF3("taue","[0]",0,300,0,300,0,300);
    TauE->SetParameter(0,5e-9);
    
    TF3 *TauH=new TF3("tauh","[0]",0,300,0,300,0,300);
    TauH->SetParameter(0,5e-9);
    
        det->TauE=TauE;
        det->TauH=TauH;
    
    TCanvas c1;
    det.SetEntryPoint(25, 299.9, 0.5);
    det.SetExitPoint(25, 1., 0.5);
    det.Temperature = 253;
    det.diff = 1;
    det.ShowMipIR(100);

    TCanvas c2;
    det.MipIR(100);
    //    det.sum->Reset();
    //    SimpleTrap(det.neg,1e-9,0,25e-9);
    //    SimpleTrap(det.pos,1e-9,0,25e-9);
    //    det.sum.Add(det.pos);    det.sum.Add(det.neg);

    det.sum->Draw();
    det.pos->Draw("SAME");
    det.neg->Draw("SAME");


}
