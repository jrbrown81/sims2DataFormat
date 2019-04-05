{
   TFile *f=new TFile("2billion_1.root","recreate");
   TTree* t=new TTree("tree","lab setup entangled sim");
   Long64_t nlines=t->ReadFile("2billion (1).csv",
   "annihil_N:theta1:phi1:E1_a:X1_a:Y1_a:Z1_a:E1_b:X1_b:Y1_b:Z1_b:theta2:phi2:E2_a:X2_a:Y2_a:Z2_a:E2_b:X2_b:Y2_b:Z2_b");
   cout << nlines << " lines found in file" << endl;
   t->Write();
   f->Close();
}
