{
    TFile *_file0 = TFile::Open("file2Bprocessed.root");
    TTree* newTree;
    _file0->GetObject("tree",newTree);
    gROOT->ProcessLine(".L sims2DataFormat.C");
    gROOT->ProcessLine("sims2DataFormat t(newTree)");
    t.Loop();
}
