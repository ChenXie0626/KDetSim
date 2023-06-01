{
    TF1 *neff = new TF1("neff", "x>[0]?[1]:0", 0, 1000);
    neff->SetParameter(0, 200);
    neff->SetParameter(1, -1);
    KPad det(50, 300);
    det.Neff = neff;
    det.Voltage = 200;
    det.SetUpVolume(1);
    det.SetUpElectrodes();
}
