#define zInvestigation_cxx
// To use this file, try the following session on your Tree T:
//
// root> T->Process("zInvestigation.C")
// root> T->Process("zInvestigation.C","some options")
// root> T->Process("zInvestigation.C+")
//

#include "zInvestigation.h"
#include "utils/pairset.cc"
#include "utils/selector.cc"
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>

void zInvestigation::Begin(TTree * /*tree*/)
{
   TString option = GetOption();

   hist_eta_P = new TH1F("","Z pseudorapidity",30,-5,5);
   hist_mass_P = new TH1F("","Z mass",60,0,150);

   total_leptons_ = 0;
   isol_rejected_ = 0;
   pt_rejected_ = 0;
   empty_lists_ = 0;
   one_pair_onlys_ = 0;
   nSignalEvents = 0;

   //hist_ranks = new TH1F("ranks_count","ranks",10,0,10);
   //hist_options = new TH1F("options_count","options",10,0,10);
}

void zInvestigation::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t zInvestigation::Process(Long64_t entry)
{
   fReader.SetEntry(entry);

   int baseline = 0; //used to shift id number of particles (i.e. their index in the array)
   //s.t. id's are unique across all particle types.
   utils::PairSet particle_pairs;
   particle_pairs.SetMetric(this->optimisingMetric);
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
         if ( !select.Pass(pt[i],eta[i],phi[i],charge[i],isol[i]) ) continue;

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

      total_leptons_ += select.totalCount;
      isol_rejected_ += select.isolCut;
      pt_rejected_ += select.ptCut; 

      //check there is at least one of each:
      if ( particles.empty() || antiparticles.empty() )
      {
         empty_lists_++;
         continue; 
      }

      if ( !select.GoodLeading() ) continue;

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
   utils::PairSet chosen_pairs = particle_pairs.GetBestNPairs(2,false);

   Int_t n_pairs = chosen_pairs.GetNPairs();
   if ( n_pairs == 0 )
   {
      //zero_pairs++;
      //std::cout << "No pairs found... " << particles.size() << ","
      //            << antiparticles.size() << ".\n";
      one_pair_onlys_++;
      return kTRUE;
   }

   //hist_ranks->Add(particle_pairs.hist_ranks);
   //hist_options->Add(particle_pairs.hist_options);

   for ( Int_t i = 0; i < n_pairs; i++ )
   {
      TLorentzVector pair_fourmomentum = chosen_pairs.Fourmomentum(i);
      hist_mass_P->Fill( pair_fourmomentum.M() );
      hist_eta_P->Fill( pair_fourmomentum.Eta() );
   }

   return kTRUE;
}

void zInvestigation::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void zInvestigation::Terminate()
{
   printf("Total leptons: %d, Lost to isolation cuts: %d, Lost to pt cuts: %d, No single pairs: %d, \
         No two pairs: %d.\n",
         total_leptons_, isol_rejected_, pt_rejected_, empty_lists_, one_pair_onlys_);

   TCanvas *c1 = new TCanvas("c1");
   hist_mass_P->Draw();

   TCanvas *c2 = new TCanvas("c2");
   hist_eta_P->Draw();
   
   TF1* fit1 = new TF1("fit1","[0]*exp(-[2]*(x-[1])**[3])",-5,5);
   fit1->SetParameters(400, 0, 1.4, 2.1);
   hist_eta_P->Fit("fit1");
}
