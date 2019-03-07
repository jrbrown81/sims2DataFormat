#define sims2DataFormat_cxx
#include "sims2DataFormat.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void sims2DataFormat::Loop()
{
// To run it on a different file, e.g.:
//                      Root > TFile *_file0 = TFile::Open("ANotherFile.root")
//                      Root > TTree* newTree
//                      Root > _file0->GetObject("tree",newTree)
//                      Root > .L sims2DataFormat.C
//                      Root > sims2DataFormat t(newTree)
//                      Root > t.Loop()

//   In a ROOT session, you can do:
//      root> .L sims2DataFormat.C
//      root> sims2DataFormat t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout << nentries << " lines to be processed from file." << endl;

// Parameters to replicate experimental configuration
   Float_t triggerThresh=50; // keV
   Bool_t errorFlags=0;
   Float_t fwhm=4.0; //energy resolution of detectors (%)
   Float_t sigma=fwhm/235.5; // faction rather than %
   
   cout << "Trigger threshold: " << triggerThresh << " keV" << endl
      << "FWHM for energy smearing: " << fwhm << "%" << endl;
   
   if(errorFlags==kTRUE) cout << endl << "Error warnings inhibited" << endl;
   
// My temporary data definations
   Int_t eventCounter=0;
   Int_t amID[2]={-1,-1};
   Int_t pixel[2][2]={{-1,-1},{-1,-1}};
   Float_t energy[2][2]={{0,0},{0,0}};
   Int_t pos[2][2]={{0,0},{0,0}};
   Int_t nTrigPixels[2]={0,0};

// Definitions for new trees for output of exp. data format
//   Long_t event = -1;
   Int_t nAMs = 2;
   Long_t *event_e = new Long_t[nAMs];
   Int_t *nTrigPixels_e = new Int_t[nAMs];
   Int_t *AM_e = new Int_t[nAMs];
   Int_t *GM_e = new Int_t[nAMs];
   Int_t *nAMsInEvent_e = new Int_t[nAMs];
   Int_t *pixel_e = new Int_t[nAMs];
   Float_t *E = new Float_t[nAMs];
   Float_t *E_neg = new Float_t[nAMs];
   Int_t *time_e = new Int_t[nAMs];
   Int_t *triggerFlag_e = new Int_t[nAMs];
   Int_t *nloop_e = new Int_t[nAMs];
   Int_t *timeDetect_e = new Int_t[nAMs];
   Int_t *timeDetectPos_e = new Int_t[nAMs];
   Float_t *Temp_e = new Float_t[nAMs];
   Int_t *pos_e = new Int_t[nAMs];
   Long_t event_eev;
   Int_t *AM_flag_eev = new Int_t[nAMs];
   Int_t *cath_flag_eev = new Int_t[nAMs];
   Int_t nAMsInEvent_eev;
   Int_t *nTrigPixels_eev = new Int_t[nAMs];
   
   Int_t *Nvar2 = new Int_t[nAMs];

   TString rfileOut = "sortedData_Ecal_simulated.root";
   TFile *f2 = new TFile(rfileOut,"recreate");
   // Event tree
   TTree *tree2_ev = new TTree("tree2_events","Event tree of E calib data");
   TTree **tree2 = new TTree*[nAMs];
   tree2_ev->Branch("event",&event_eev,"event/L");
   tree2_ev->Branch("nAMsInEvent",&nAMsInEvent_eev,"nAMsInEvent/I");
   tree2_ev->Branch("nTrigPixels",nTrigPixels_eev,Form("nTrigPixels[%d]/I",nAMs));
   tree2_ev->Branch("AM_flag",AM_flag_eev,Form("AM_flag[%d]/I",nAMs));
   tree2_ev->Branch("cath_flag",cath_flag_eev,Form("cath_flag[%d]/I",nAMs));
   // AM trees
   for (Int_t im = 0; im < nAMs; im++)
   {
      tree2[im] = new TTree(Form("tree2_AM%d",im), Form("E calib data tree AM%d",im));
      tree2[im]->Branch("event",&event_e[im],"event/L");
      tree2[im]->Branch("nTrigPixels",&nTrigPixels_e[im],"nTrigPixels/I");
      tree2[im]->Branch("AM",&AM_e[im],"AM/I");
      tree2[im]->Branch("GM",&GM_e[im],"GM/I");
      tree2[im]->Branch("nAMsInEvent",&nAMsInEvent_e[im],"nAMsInEvent/I");
      tree2[im]->Branch("pixel",&pixel_e[im],"pixel/I");
      tree2[im]->Branch("E",&E[im],"E/F");
      tree2[im]->Branch("E_neg",&E_neg[im],"E_neg/F");
      tree2[im]->Branch("pos",&pos_e[im],"pos/I");
      tree2[im]->Branch("time",&time_e[im],"time/I");
      tree2[im]->Branch("triggerFlag",&triggerFlag_e[im],"triggerFlag/I");
      tree2[im]->Branch("nloop",&nloop_e[im],"nloop/I");
      tree2[im]->Branch("timeDetect",&timeDetect_e[im],"timeDetect/I");
      tree2[im]->Branch("timeDetectPos",&timeDetectPos_e[im],"timeDetectPos/I");
      tree2[im]->Branch("Temp",&Temp_e[im],"Temp/F");
      Nvar2[im] = tree2[im]->GetNbranches();
   }
// End of definitions for output file

// Main loop over input file
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      eventCounter++;
      
//      cout << event_eev << " " << nAMsInEvent_eev << " " << nTrigPixels_eev[0] << " " << nTrigPixels_eev[1] << " " << AM_flag_eev[0] << " " << AM_flag_eev[1] << " " << cath_flag_eev[0] << " " << cath_flag_eev[1]  << endl;

// Find which AMs and pixels hit
// Find AMs
      amID[0]=findAM(Z1_a);
      amID[1]=findAM(Z2_a);
      if(amID[0]==-1 || amID[1]==-1) {
         if(errorFlags==kTRUE) cout << "Error: Gamma ray not in a detector!" << endl;
         continue;
      }
      if(amID[0]==amID[1]) {
         if(errorFlags==kTRUE) cout << "Error: Same AM hit by both gammas!" << endl;
         continue;
      }
// Find pixels
      pixel[0][0]=findPixelID(amID[0],X1_a,Y1_a);
      pixel[0][1]=findPixelID(amID[0],X1_b,Y1_b);
      pixel[1][0]=findPixelID(amID[1],X2_a,Y2_a);
      pixel[1][1]=findPixelID(amID[1],X2_b,Y2_b);

//      if(pixel[0][0]==-1 || pixel[0][1]==-1) {
//         if(errorFlags==kTRUE) cout << "Error: Gamma ray not in a pixel! " << endl;
//         continue;
//      }

// Smear energies
      energy[0][0]=gRandom->Gaus(E1_a,sigma*E1_a);
      energy[0][1]=gRandom->Gaus(E1_b,sigma*E1_b);
      energy[1][0]=gRandom->Gaus(E2_a,sigma*E2_a);
      energy[1][1]=gRandom->Gaus(E2_b,sigma*E2_b);

// Find ordering and number of triggered pixels
      for(int im=0; im<2; im++) {
         if(energy[im][0]>energy[im][1] && pixel[im][0]!=-1) {
            pos[im][0]=0;
            pos[im][1]=1;
         } else {
            pos[im][0]=1;
            pos[im][1]=0;
         }
         if(pixel[im][0]==pixel[im][1]) {
            energy[im][pos[im][0]]+=energy[im][pos[im][1]];
            energy[im][pos[im][1]]=0;
         }
         if(energy[im][pos[im][1]]>triggerThresh && pixel[im][1]!=-1) {
            nTrigPixels[im]=2;
         }
         else if(energy[im][pos[im][0]]>triggerThresh && pixel[im][0]!=-1) nTrigPixels[im]=1;
         else nTrigPixels[im]=0;; // skip if no pixels triggered in either AM
      }
      if(nTrigPixels[0]==0 && nTrigPixels[1]==0) continue;
//      if(pixel[0][0]==pixel[0][1]) {
////         cout << "Warning: Same pixel hit twice! " << endl;
//         energy[0][0]+=energy[0][1];
//         pixel[0][1]=-1;
////         continue;
//      }
      
//      cout << amID[0] << " " << pixel[0][0] << " " << pixel[0][1] << " " << energy[0][0] << " " << energy[0][1] << " "
//         << amID[1] << " " << pixel[1][0] << " " << pixel[1][1] << " " << energy[1][0] << " " << energy[1][1] << endl;
      
// Fill data for AM trees
      for(Int_t im = 0; im < nAMs; im++) {
         event_e[amID[im]]=eventCounter-1;
         nTrigPixels_e[amID[im]]=nTrigPixels[im];
//         nTrigPixels_e[amID[im]]=nTrigPixels[im];
         AM_e[amID[im]]=amID[im];
         GM_e[amID[im]]=1;
         nAMsInEvent_e[amID[im]]=2;
         pixel_e[amID[im]]=pixel[im][pos[im][0]];
         E[amID[im]]=energy[im][pos[im][0]];
//         if(E[amID[im]]<triggerThresh) cout << "Error! Error! " << event_e[amID[im]] << " " << amID[im] << " " << nTrigPixels_e[amID[im]] << endl;
         E_neg[amID[im]]=0;
         time_e[amID[im]]=0;
         triggerFlag_e[amID[im]]=0;
         nloop_e[amID[im]]=0;
         timeDetect_e[amID[im]]=0;
         timeDetectPos_e[amID[im]]=0;
         Temp_e[amID[im]]=0;
         pos_e[amID[im]]=0;
         if(nTrigPixels_e[amID[im]]!=0) tree2[amID[im]]->Fill();
         // if two pixel hits, every thing same except pixel and energy
         if(nTrigPixels_e[amID[im]]==2) {
            pixel_e[amID[im]]=pixel[im][pos[im][1]];
            E[amID[im]]=energy[im][pos[im][1]];
            pos_e[amID[im]]=1;
            tree2[amID[im]]->Fill();
         }
      }
      
      // Fill data for event tree
      event_eev=eventCounter-1;
      nAMsInEvent_eev=2;
      for(Int_t im = 0; im < nAMs; im++) {
         nTrigPixels_eev[amID[im]]=nTrigPixels_e[amID[im]];
         AM_flag_eev[im]=1;
         cath_flag_eev[im]=0;
      }
      tree2_ev->Fill();
      
   }
// End of loop over input file
   
// Write trees and close files
   f2->cd();
   tree2_ev->Write();
   for (Int_t im = 0; im < nAMs; im++)
   {
      tree2[im]->Write();
   }
   f2->Close();
   cout << "Output written to file: " << rfileOut << endl;
   
} // End of sims2DataFormat::Loop()

Int_t findAM(Float_t z)
{
   Int_t am=-1;
   if(z>=80) am=0;
   else if(z<=-80) am=1;
   
   return am;
}

Int_t findPixelID(Int_t am, Float_t x, Float_t y)
{
   Int_t pixelID=-1;
   Int_t row=-1, col=-1;
   
   // AM 1 if rotated 180 degrees so invert x
   if(am==1) x=-x;
   
   // added 100um to extremes
   Float_t pixEdges[12]={-4.5,-3.6,-2.8,-2.0,-1.2,-0.4,0.4,1.2,2.0,2.8,3.6,4.5};
   for(Int_t i=0; i<12; i++) {
      if(x>pixEdges[i] && x<=pixEdges[i+1]) col=i+1;
      if(y>pixEdges[i] && y<=pixEdges[i+1]) row=11-i;
   }
//   if(row!=-1 && col!=-1) pixelID=col+(row-1)*11;
   if(row!=-1 && col!=-1) pixelID=col+(row-1)*44;
   
//   if(pixelID!=-1){
   if(am==0) return pixelID;
   else if(am==1 && pixelID!=-1) return pixelID+22;
   else return pixelID;
}

