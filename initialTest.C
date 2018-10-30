#define initialTest_cxx
// To use this file, try the following session on your Tree T:
//
// root> T->Process("initialTest.C")
// root> T->Process("initialTest.C","some options")
// root> T->Process("initialTest.C+")
//

#include "initialTest.h"
#include "utils/pairset.cc"
#include "utils/selector.cc"
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>

void initialTest::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   hist_combined_mass_P = new TH1F("","di-Z combined",60,80,500);
   hist_pair_mass_P = new TH1F("","dilepton",60,0,150);

   for ( std::string lepton : { "Muon", "Electron" } )
   { 
      n_map.insert(
            std::pair< const std::string, TTreeReaderValue<Int_t>* >
            ( lepton, new TTreeReaderValue<Int_t>
              ({ fReader, (lepton + std::string("_n")).c_str() })
            ));
      pt_map.insert(
            std::pair< const std::string, TTreeReaderArray<Float_t>* >
            ( lepton, new TTreeReaderArray<Float_t>
              ({ fReader, (lepton + std::string("_Pt")).c_str() })
            ));
      eta_map.insert(
            std::pair< const std::string, TTreeReaderArray<Float_t>* >
            ( lepton, new TTreeReaderArray<Float_t>
              ({ fReader, (lepton + std::string("_Eta")).c_str() })
            ));
      phi_map.insert(
            std::pair< const std::string, TTreeReaderArray<Float_t>* >
            ( lepton, new TTreeReaderArray<Float_t>
              ({ fReader, (lepton + std::string("_Phi")).c_str() })
            ));
      charge_map.insert(
            std::pair< const std::string, TTreeReaderArray<Float_t>* >
            ( lepton, new TTreeReaderArray<Float_t>
              ({ fReader, (lepton + std::string("_Charge")).c_str() })
            ));
      isol_map.insert(
            std::pair< const std::string, TTreeReaderArray<Float_t>* >
            ( lepton, new TTreeReaderArray<Float_t>
              ({ fReader, (lepton + std::string("_Isol")).c_str() })
            ));
   }
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
   fReader.SetEntry(entry);

   int baseline = 0; //used to shift id number of particles (i.e. their index in the array)
   //s.t. id's are unique across all particle types.
   utils::PairSet particle_pairs;
   //particle_pairs.SetBenchmark(Z_MASS);

   //========================MAP PAIRS=============================
   //for each entry, check for pairs of muons/electrons
   for ( std::string lepton : { "Muon", "Electron" } )
   { 
      Int_t n = *(*n_map[lepton]);
      std::vector<Float_t> pt, eta, phi, charge, isol;
      for (Int_t i = 0; i < n; i++)
      {
         pt.push_back(pt_map[lepton]->At(i));
         eta.push_back(eta_map[lepton]->At(i));
         phi.push_back(phi_map[lepton]->At(i));
         charge.push_back(charge_map[lepton]->At(i));
         isol.push_back(isol_map[lepton]->At(i));
      }
      
      //============FILL=============
      //fill lists of particles and anti-particles. Lists contain the indices
      //of the objects they're referencing.
      //Lepton selection is done during this process using the Selector class.
      std::list<Int_t> particles, antiparticles;
      utils::Selector select(lepton);
      for ( Int_t i = 0; i < n; i++ )
      {
         //LEPTON SELECTION:
         //if ( !select.Pass(pt[i],eta[i],phi[i],charge[i],isol[i]) ) continue;

         if ( charge[i] == -1 )
         {
            particles.push_back(i);
         }
         else
         {
            antiparticles.push_back(i);
         }
      }
      //=============================
      
      //check there is at least one of each:
      if ( particles.empty() || antiparticles.empty() ) continue; 
      
      //if ( !select.GoodLeading() ) continue;

      Float_t mass = (lepton == "Muon") ? MUON_MASS : ELECTRON_MASS;

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
   if ( n_pairs == 0 )
   {
      //zero_pairs++;
      //std::cout << "No pairs found... " << particles.size() << ","
      //            << antiparticles.size() << ".\n";
      return kTRUE;
   }

   TLorentzVector total_combined_fourmomentum(0,0,0,0);
   for ( Int_t i = 0; i < n_pairs; i++ )
   {
      TLorentzVector pair_fourmomentum = chosen_pairs.Fourmomentum(i);
      hist_pair_mass_P->Fill( pair_fourmomentum.M() );
      total_combined_fourmomentum += pair_fourmomentum; 
   }

   if ( ( chosen_pairs.M(1) - 91.19 ) < 10 )
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
