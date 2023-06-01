{
    TF3 *f2 = new TF3("f2", "[0]", 0, 3000, 0, 3000, 0, 3000);
    f2->SetParameter(0, -2);

    KStrip *det = new KStrip(80, 10, 2, 3, 300);
    det->Voltage = -400;
    det->SetUpVolume(2);

    det->SetUpElectrodes();
    det->SetBoundaryConditions();
    det->SetUpMaterial(0); //Silicon
    //det->SetUpMaterial(10); //diamond

    det->NeffF = f2;
    // det->Mat=10;
    det->SetDebug(0);
    det->CalField(0);
    det->CalField(1);

    det->SetEntryPoint(25, 1, 0.5);
    det->SetExitPoint(25, 300, 0.5);

    det->diff = 1;

  det.SetEntryPoint(70,1,0.5);
  det.SetExitPoint(70,299,0.5);
  det.ShowMipIR(200);
   det.MipIR(200);
   printf("integral %f\n",det.sum.Integral());
SimpleTrap(det.neg,1e-9,0,10e-9);
det.neg.Draw();
SimpleTrap(det.pos,2e-9,0,10e-9);
det.pos.Draw("SAME");
det.sum.Reset();
det.sum.Add(det.pos);
det.sum.Add(det.neg);
det.sum.Draw("SAME");
 printf("integral %f\n",det.sum.Integral());

}
