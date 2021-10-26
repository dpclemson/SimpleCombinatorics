void callHK(int,int);
void OHist(){
  callHK(2,2);
  callHK(2,3);
  callHK(2,4);
  callHK(2,5);
  callHK(2,6);
  callHK(2,7);
  callHK(2,8);
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
  TCanvas c1("c1","");
  histogram->SetMarkerStyle(kFullCircle);
  histogram->SetMarkerColor(kBlack);
  histogram->Draw("ex0p");
  histogram->GetXaxis()->SetTitle("N_{particles}");
  histogram->GetYaxis()->SetTitle(Form("#LT%d#GT_{%d}",K,H));
  c1.SetLogy();
  c1.SetGrid(1);
  c1.Print(Form("Figures/TestHist_h%d_k%d.png",H,K));
}
