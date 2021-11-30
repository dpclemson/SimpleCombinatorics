void callHK(int,int);
void OHist(){
  callHK(2,2);
  callHK(2,3);
  callHK(2,4);
  callHK(2,5);
  callHK(2,6);
  callHK(2,7);
  callHK(2,8);
  callHK(3,2);
  callHK(3,3);
  callHK(3,4);
  callHK(3,5);
  callHK(3,6);
  callHK(3,7);
  callHK(3,8);
}
void callHK(int H,int K){
  TFile* file = TFile::Open(Form("CondorOutput/Combined_h%d_k%d.root",H,K));
  if ( file == NULL ){
    cout << "Warning Missing File!"<< endl;
    return;
  }
  TProfile* histogram = (TProfile*)file->Get(Form("hmult_recursion_0_%d",K-2));
  if ( histogram == NULL ){
    cout << "Warning Histogram Missing!" << endl;
    return;
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas c1("c1","",1000,500);
  c1.Divide(2,1);
  c1.cd(1);
  histogram->SetMarkerStyle(kFullCircle);
  histogram->SetMarkerColor(kBlack);
  histogram->Draw("ex0p");
  histogram->GetXaxis()->SetTitle("N_{particles}");
  histogram->GetYaxis()->SetTitle(Form("#LT%d#GT_{%d}",K,H));
  //c1.SetLogy();
  //c1_1.SetLogy();
  //c1_1.SetGrid(1);
  gPad->SetLogy();
  gPad->SetGrid(1);
  //c1.Print(Form("Figures/TestHist_h%d_k%d.png",H,K));
  TF1 fun("fun","[0]/pow(x,[1])",K,250);
  fun.SetParameter(0,1.0);
  if(K==4)fun.FixParameter(1,K-1);
  if(K==4)fun.SetParLimits(1,K-2,K);
  fun.SetParameter(1,K-1);
  histogram->Fit("fun","R");
  TLatex tex;
  tex.SetNDC();
  tex.DrawLatex(0.6,0.75,Form("p0=%f",fun.GetParameter(0)));
  tex.DrawLatex(0.6,0.7,Form("p1=%f",fun.GetParameter(1)));
  tex.DrawLatex(0.6,0.65,Form("k-1=%d",K-1));
  //c1.Print(Form("Figures/TestHistFit_h%d_k%d.png",H,K));
  const int nbins=histogram->GetNbinsX();
  double residual[nbins];
  double Ntrack[nbins];
  for (int i=0; i < nbins; ++i )
    {
      Ntrack[i]=histogram->GetBinCenter(i+1);
      residual[i]=1;
      double fun_val=fun.Eval(histogram->GetBinCenter(i+1));
      residual[i]=histogram->GetBinContent(i+1) - fun_val;
      residual[i] /= fun_val;
      //      cout<<Ntrack[i]<<" "<<fun_val<<" "<<residual[i]<<endl;
    }
  c1.cd(2);
  //c1_2.SetLogy(0);
  //c1_2.SetGrid(1);
  gPad->SetGrid(1);
  residual[0]=1;
  TGraphErrors tge_residual(nbins, Ntrack, residual, 0, 0);
  tge_residual.SetMarkerStyle(kOpenCircle);
  tge_residual.SetMarkerColor(kGreen);
  tge_residual.Draw("ap");
  tge_residual.SetMaximum(0.99);
  tge_residual.SetMinimum(-0.99);
  c1.Print(Form("Figures/TestHist_h%d_k%d_residuals.png",H,K));
  
}

