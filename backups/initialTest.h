//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  3 17:00:11 2018 by ROOT version 6.14/04
// from TTree EventTree/
// found on file: Output.ZZ.Tree.root
//////////////////////////////////////////////////////////

#ifndef initialTest_h
#define initialTest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <vector>
#include <string>

//=======CUSTOM DEFINITIONS========
#define MUON_MASS 0.10565837 // with no significant error at this precision
#define ELECTRON_MASS 0.0005109989
#define Z_MASS 91.19

struct Z_like
{
   bool operator() ( const std::pair<Int_t,TLorentzVector>& lhs,
                     const std::pair<Int_t,TLorentzVector>& rhs) const
   {
      return fabs(lhs.second.M() - Z_MASS) < fabs(rhs.second.M() - Z_MASS);
   }
};
//=================================

class initialTest : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TH1F *hist_combined_mass_P, *hist_pair_mass_P;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Float_t> Event_Weight = {fReader, "Event_Weight"};
   TTreeReaderValue<Float_t> MET_Et = {fReader, "MET_Et"};
   TTreeReaderValue<Float_t> MET_Phi = {fReader, "MET_Phi"};
   TTreeReaderValue<Int_t> Jet_n = {fReader, "Jet_n"};
   TTreeReaderArray<float> Jet_Pt = {fReader, "Jet_Pt"};
   TTreeReaderArray<float> Jet_Eta = {fReader, "Jet_Eta"};
   TTreeReaderArray<float> Jet_Phi = {fReader, "Jet_Phi"};
   TTreeReaderArray<float> Jet_M = {fReader, "Jet_M"};
   TTreeReaderArray<int> Jet_Flav = {fReader, "Jet_Flav"};
   TTreeReaderValue<Int_t> Muon_n = {fReader, "Muon_n"};
   TTreeReaderArray<float> Muon_Pt = {fReader, "Muon_Pt"};
   TTreeReaderArray<float> Muon_Eta = {fReader, "Muon_Eta"};
   TTreeReaderArray<float> Muon_Phi = {fReader, "Muon_Phi"};
   TTreeReaderArray<float> Muon_Isol = {fReader, "Muon_Isol"};
   TTreeReaderArray<float> Muon_Charge = {fReader, "Muon_Charge"};
   TTreeReaderValue<Int_t> Electron_n = {fReader, "Electron_n"};
   TTreeReaderArray<float> Electron_Pt = {fReader, "Electron_Pt"};
   TTreeReaderArray<float> Electron_Eta = {fReader, "Electron_Eta"};
   TTreeReaderArray<float> Electron_Phi = {fReader, "Electron_Phi"};
   TTreeReaderArray<float> Electron_Isol = {fReader, "Electron_Isol"};
   TTreeReaderArray<float> Electron_Charge = {fReader, "Electron_Charge"};


   initialTest(TTree * /*tree*/ =0) { }
   virtual ~initialTest() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(initialTest,0);

};

#endif

#ifdef initialTest_cxx
void initialTest::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t initialTest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef initialTest_cxx
