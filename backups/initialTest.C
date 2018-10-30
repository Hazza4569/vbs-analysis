#define initialTest_cxx
// To use this file, try the following session on your Tree T:
//
// root> T->Process("initialTest.C")
// root> T->Process("initialTest.C","some options")
// root> T->Process("initialTest.C+")
//


#include "initialTest.h"
#include "utils/pairset.cc"
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>

void initialTest::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   hist_combined_mass_P = new TH1F("","di-Z combined",60,0,150);
   hist_pair_mass_P = new TH1F("","dilepton",60,0,150);
}

void initialTest::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t initialTest::Process(Long64_t entry)
{

   std::vector<std::string> leptons = {"Muon","Electron"};

   int baseline = 0; //used to shift id number of particles (i.e. their index in the array)
                     //s.t. id's are unique across all particle types.
                     
   utils::PairSet particle_pairs;

   //============READ=============
   std::map< std::string, Int_t> n_map;
   std::map< std::string, std::vector<Float_t> > pt_map;
   std::map< std::string, std::vector<Float_t> > eta_map;
   std::map< std::string, std::vector<Float_t> > phi_map;
   std::map< std::string, std::vector<Float_t> > charge_map;

   for (int iLepton = 0; iLepton < leptons.size(); iLepton++)
   {
      const std::string key = leptons[iLepton];
      TTreeReaderValue<Int_t> 
         n_read = { fReader, (leptons[iLepton] + std::string("_n")).c_str() };
      TTreeReaderArray<Float_t>
         pt_read = { fReader, (leptons[iLepton] + std::string("_Pt")).c_str() };
      TTreeReaderArray<Float_t>
         eta_read = { fReader, (leptons[iLepton] + std::string("_Eta")).c_str() };
      TTreeReaderArray<Float_t>
         phi_read = { fReader, (leptons[iLepton] + std::string("_Phi")).c_str() };
      TTreeReaderArray<Float_t>
         charge_read = { fReader, (leptons[iLepton] + std::string("_Charge")).c_str() };
      fReader.SetEntry(entry);

      Int_t n_val = *n_read;
      std::vector<Float_t> pt_vec, eta_vec, phi_vec, charge_vec;
      for (Int_t i = 0; i < n_val; i++)
      {
         pt_vec.push_back(pt_read[i]);
         eta_vec.push_back(eta_read[i]);
         phi_vec.push_back(phi_read[i]);
         charge_vec.push_back(charge_read[i]);
      }

      n_map.insert(
         std::pair< const std::string, Int_t > ( key, n_val )
      );
      pt_map.insert(
         std::pair< const std::string, std::vector<Float_t> > ( key, pt_vec )
      );
      eta_map.insert(
         std::pair< const std::string, std::vector<Float_t> > ( key, eta_vec )
      );
      phi_map.insert(
         std::pair< const std::string, std::vector<Float_t> > ( key, phi_vec )
      );
      charge_map.insert(
         std::pair< const std::string, std::vector<Float_t> > ( key, charge_vec )
      );        
   }
   //=============================
   
   //========================MAP PAIRS=============================
   //for each entry, check for pairs of muons/electrons
   for (int iLepton = 0; iLepton < leptons.size(); iLepton++)
   { 
      std::string l = leptons[iLepton];
      Int_t n = n_map[l];
      std::vector<Float_t> pt(pt_map[l]), eta(eta_map[l]), phi(phi_map[l]), charge(charge_map[l]);

      //============FILL=============
      //fill lists of particles and anti-particles. Lists contain the indices
      //of the objects they're referencing:
      std::list<Int_t> particles, antiparticles;
      for ( Int_t i = 0; i < n; i++ )
      {
         if ( fabs(charge[i]) != 1 )
         {
            std::cout << "ERROR: Expected lepton but read non-unit charge.\n";
            continue;
         }
         if ( charge[i] == -1 ) particles.push_back(i);
         else antiparticles.push_back(i);
      }
      //=============================

      //check there is at least one of each:
      if ( particles.empty() || antiparticles.empty() ) continue; 

      Float_t mass = (l == "Muon") ? MUON_MASS : ELECTRON_MASS;

      //iterate over particles and fill map of particles and candidate antiparticles:
      for ( std::list<Int_t>::iterator iP = particles.begin(); iP != particles.end(); iP++ )
      {
         TLorentzVector particle_fourmomentum;
         particle_fourmomentum.SetPtEtaPhiM( pt[*iP], eta[*iP], phi[*iP], mass );

         for ( std::list<Int_t>::iterator iAP = antiparticles.begin();
               iAP != antiparticles.end(); iAP++)
         {
            TLorentzVector antiparticle_fourmomentum; 
            antiparticle_fourmomentum.SetPtEtaPhiM( pt[*iAP], eta[*iAP], phi[*iAP], mass );

            TLorentzVector combined_fourmomentum = (particle_fourmomentum + antiparticle_fourmomentum);

            particle_pairs.AddPair( (*iP)+baseline, (*iAP)+baseline, combined_fourmomentum );
         }
      }
      baseline += n;
   }
   //==============================================================

   //=======================FIND OPTIMUM===========================
   utils::PairSet chosen_pairs = particle_pairs.GetBestNPairs(2);

   Int_t n_pairs = chosen_pairs.GetNPairs();
   if ( n_pairs == 0 ) return kTRUE;

   TLorentzVector total_combined_fourmomentum(0,0,0,0);
   for ( Int_t i = 0; i < n_pairs; i++ )
   {
      TLorentzVector pair_fourmomentum = chosen_pairs.Fourmomentum(i);
      hist_pair_mass_P->Fill( pair_fourmomentum.M() );
      total_combined_fourmomentum += pair_fourmomentum; 
   }

   hist_combined_mass_P->Fill( total_combined_fourmomentum.M() );

   return kTRUE;
}

void initialTest::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void initialTest::Terminate()
{
   TCanvas *c1 = new TCanvas("c1");
   hist_combined_mass_P->Draw();

   TCanvas *c2 = new TCanvas("c2");
   hist_pair_mass_P->Draw();
}
