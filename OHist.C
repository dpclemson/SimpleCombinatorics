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
  TFile* file = TFile::Open(Form("OutputFiles/OutFile_h%d_k%d.root",H,K));
  if ( file == NULL ){
    cout << "Warning Missing File!"<< endl;
    return;
  }
  TProfile* histogram = (TProfile*)file->Get(Form("hmult_recursion_0_%d",K-2));
  if ( histogram == NULL ){
    cout << "Warning Histogram Missing!" << endl;
    return;
  }
  TCanvas c1("c1","");
  histogram->SetMarkerStyle(kFullCircle);
  histogram->SetMarkerColor(kBlack);
  histogram->Draw("ex0p");
  c1.SetLogy();
  c1.Print(Form("Figures/TestHist_h%d_k%d.png",H,K));
}
