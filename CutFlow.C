#define MUON_MASS 0.10565837
#define ELECTRON_MASS 0.0005109989
#define Z_MASS 91.19
#include "CutFlow.h"
#include "utils/selector.cc"
#include "utils/pairset.cc"

CutFlow::CutFlow(std::string signalFile, std::string bgFile)
{
   std::string path("/home/user108/y4p/root/"), rootFt(".root"), datFt(".dat");

   files_.insert( make_pair("signal", signalFile) );
   files_.insert( make_pair("background", bgFile) );

   for ( std::string sample : { "signal", "background" } )
   {
      //get cross section data and number of events:
      Int_t n; Double_t sigma;
      std::fstream fs;
      fs.open( (path+files_[sample]+datFt).c_str() );
      fs >> n;
      fs >> sigma; 
      fs.close();
      num_events_.insert( make_pair(sample, n) );
      cross_section_.insert( make_pair(sample, sigma) );

      //Get sig and bg trees:
      TFile* f = TFile::Open( (path + files_[sample] + rootFt).c_str() );
      trees_.insert( make_pair(sample, (TTree*)f->Get("EventTree")) );
      readers_.insert( make_pair(sample, new TTreeReader(trees_[sample])) );
   }

   //Put all variables from file into one big map with a string seperated by underscores;
   for ( std::string sample : { "signal", "background" } )
   {
      for ( std::string obj_type : { "Muon", "Electron", "Jet" } )
      {
         std::vector<std::string> vars = { "Pt", "Eta", "Phi" };

         if ( obj_type == "Muon" || obj_type == "Electron" )
            vars.insert( vars.end(), { "Isol", "Charge" } );

         else if ( obj_type == "Jet" )
            vars.insert( vars.end(), { "M" } );

         for ( std::string variable : vars )
         {
            event_vars_.insert( make_pair(
                     sample+delim_+obj_type+delim_+variable,
                     new TTreeReaderArray<Float_t>({
                        *readers_[sample],
                        (obj_type+delim_+variable).c_str()
                        })
                     ) );     
         }
         object_numbers_.insert( make_pair(
                     sample+delim_+obj_type,
                     new TTreeReaderValue<Int_t>({
                        *readers_[sample],
                        (obj_type+std::string("_n")).c_str()
                        })
                     ) );
      } 
      fileEventNums_.insert( std::make_pair(sample,
         new std::fstream()) );
      eventFileName_.insert( std::make_pair( sample,
         std::string(".event_nums_")+sample+std::string(".dat") ));
   } 
   FileInit();
}

CutFlow::~CutFlow()
{
   //for (auto& tree : trees_)
   //{
   //delete tree.second;
   //}
}

void
CutFlow::PreSelection( Int_t n_on_shell )
{
   for ( std::string sample : { "signal", "background" } )
   {
      fileEventNums_[sample]->open(eventFileName_[sample],ios::out);
      for ( Long64_t entry = 0; entry < num_events_[sample]; entry++ )
      {
         readers_[sample]->SetEntry(entry);

         int baseline = 0; //used to shift id number of particles (i.e. their index in the array)
         //s.t. id's are unique across all particle types.
         utils::PairSet particle_pairs;
         particle_pairs.SetMetric("smds");
         particle_pairs.SetBenchmark(Z_MASS);

         //lepton selection
         for ( std::string object : { "Muon", "Electron" } )
         {
            Int_t n = *(*object_numbers_[sample+delim_+object]);
            std::vector<Float_t> pt, eta, phi, charge, isol;
            for (Int_t i = 0; i < n; i++)
            {  
               pt.push_back(event_vars_[sample+delim_+object+delim_+std::string("Pt")]->At(i));
               eta.push_back(event_vars_[sample+delim_+object+delim_+std::string("Eta")]->At(i));
               phi.push_back(event_vars_[sample+delim_+object+delim_+std::string("Phi")]->At(i));
               charge.push_back(event_vars_[sample+delim_+object+delim_+std::string("Charge")]->At(i));
               isol.push_back(event_vars_[sample+delim_+object+delim_+std::string("Isol")]->At(i));

            }

            //============FILL=============
            //fill lists of particles and anti-particles. Lists contain the indices
            //of the objects they're referencing.
            //Lepton selection is done during this process using the Selector class.
            std::list<Int_t> particles, antiparticles;
            utils::Selector select(object);
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

            //check there is at least one of each:
            if ( particles.empty() || antiparticles.empty() )
            {
               continue;
            }

            //if ( !select.GoodLeading() ) continue;
            Float_t mass = (object == "Muon") ? MUON_MASS : ELECTRON_MASS;

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

                  TLorentzVector combined_fourmomentum
                     = (particle_fourmomentum + antiparticle_fourmomentum);

                  particle_pairs.AddPair( (*iP)+baseline, (*iAP)+baseline, combined_fourmomentum );
               }
            }
            baseline += n;
         }

         //Z selection
         utils::PairSet chosen_pairs = particle_pairs.GetBestNPairs(2);

         Int_t n_pairs = chosen_pairs.GetNPairs();
         if ( n_pairs == 0 )
         {
            continue;
         }

         if ( n_on_shell > 0 && n_on_shell < 3 )
         { 
            if ( fabs(chosen_pairs.M(n_on_shell-1) - Z_MASS) > 10 ) continue;
         }

         //Add to file
         (*fileEventNums_[sample]) << entry << '\t';
      }
      fileEventNums_[sample]->close();
   }
   filedEvents_ = true;   
}

   void
CutFlow::JetEta(bool tag)
{
   std::string obj("Jet"), var("Eta");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key = sample+delim_+obj+delim_+var;
      char strHist[200];
      Int_t nBins(50);
      Double_t etaMin(-5), etaMax(5);
      snprintf(strHist,200,"%s;#eta_{j};Events / %.1f",files_[sample].c_str(),(etaMax-etaMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,etaMin,etaMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);
         if (!tag)
         {
            for ( auto & jet_eta : *event_vars_[key] ) h->Fill(jet_eta);
         }
         else
         {
            std::vector<TLorentzVector> vecs = ReconstructVectors(sample,obj);

            if ( vecs.size() < 2 ) continue;

            std::sort(vecs.begin(), vecs.end(), pt_ordered());

            //            std::vector<Float_t> jet_etas;
            //            for ( auto & jet_eta : *event_vars_[key] ) jet_etas.push_back(jet_eta);
            //
            //            Int_t n = jet_etas.size();
            //            if ( n < 2 ) continue;
            //
            //            std::sort(jet_etas.begin(), jet_etas.end());
            for ( Int_t i = 0; i < 2; i++) h->Fill(vecs[i].Eta());
         }
      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s_%s.pdf",
            files_[sample].c_str(),obj.c_str(),var.c_str(),tag ? "tag":"all");
      gPad->SaveAs(strOutput);
   }
}

   void
CutFlow::JetEtaDiff()
{
   Double_t rec_cut_value = 2.5;

   std::string obj("Jet"), var("Eta");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key = sample+delim_+obj+delim_+var;
      char strHist[200];
      Int_t nBins(45);
      Double_t xMin(0), xMax(9);
      snprintf(strHist,200,"%s;#Delta#eta_{jj};Events / %.1f",
            files_[sample].c_str(),(xMax-xMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,xMin,xMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);

         std::vector<TLorentzVector> vecs = ReconstructVectors(sample,obj);

         if ( vecs.size() < 2 ) continue;

         std::sort(vecs.begin(), vecs.end(), pt_ordered());

         //         std::vector<Float_t> jet_etas;
         //         for ( auto & jet_eta : *event_vars_[key] ) jet_etas.push_back(jet_eta);
         //
         //         Int_t n = jet_etas.size();
         //         if ( n < 2 ) continue;
         //
         //         std::sort(jet_etas.begin(), jet_etas.end());
         h->Fill(fabs(vecs[0].Eta() - vecs[1].Eta()));
      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      gPad->Update(); 
      printf("%f %f %f %f\n",rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      TLine *l = new TLine(rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      l->SetLineColor(kRed);
      l->SetLineWidth(2);
      l->SetLineStyle(kDashed);
      l->Draw();

      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s.pdf",
            files_[sample].c_str(),obj.c_str(),"EtaDiff");
      gPad->SaveAs(strOutput);
   }
}

   void
CutFlow::JetPt(Int_t n)
{
   Double_t rec_cut_value = (n == 1 || n == 2) ? 40 : -999999;

   std::string obj("Jet"), var("Pt");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key = sample+delim_+obj+delim_+var;
      char strHist[200];
      Int_t nBins(100);
      Double_t xMin(0), xMax(500);
      snprintf(strHist,200,"%s;%s%sjet p_{T} [GeV];Events / %.1f GeV",
            files_[sample].c_str(), (n==2) ? "sub" : "",(n==1||n==2) ? "leading ":"",
            (xMax-xMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,xMin,xMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);
         if (n == 0)
         {
            for ( auto & jet_pt : *event_vars_[key] ) h->Fill(jet_pt);
         }
         else
         {
            std::vector<Float_t> jet_pts;
            for ( auto & jet_pt : *event_vars_[key] ) jet_pts.push_back(jet_pt);

            Int_t nJets = jet_pts.size();
            if ( nJets < n ) continue;

            std::sort(jet_pts.begin(), jet_pts.end(), greater<float>());
            h->Fill(jet_pts[n-1]);
         }

      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      gPad->Update(); 
      printf("%f %f %f %f\n",rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      TLine *l = new TLine(rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      l->SetLineColor(kRed);
      l->SetLineWidth(2);
      l->SetLineStyle(kDashed);
      l->Draw();

      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s_%d.pdf",
            files_[sample].c_str(),obj.c_str(),var.c_str(),n);
      gPad->SaveAs(strOutput);
   }
}

   void
CutFlow::DijetM()
{
   Double_t rec_cut_value = 500;

   std::string obj("Jet");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key[4];
      Int_t index = 0;
      for (std::string var : {"Pt","Eta","Phi","M"})
      {
         key[index] = sample+delim_+obj+delim_+var;
         index++;
      }
      char strHist[200];
      Int_t nBins(200);
      Double_t xMin(0), xMax(1000);
      snprintf(strHist,200,"%s;m_{jj} [GeV];Events / %.1f GeV",
            files_[sample].c_str(),(xMax-xMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,xMin,xMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);
         std::vector<TLorentzVector> vecs = ReconstructVectors(sample,obj);

         if ( vecs.size() < 2 ) continue;

         std::sort(vecs.begin(), vecs.end(), pt_ordered());
         h->Fill( (vecs[0] + vecs[1]).M() );
      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      gPad->Update(); 
      TLine *l = new TLine(rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      l->SetLineColor(kRed);
      l->SetLineWidth(2);
      l->SetLineStyle(kDashed);
      l->Draw();

      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s.pdf",
            files_[sample].c_str(),obj.c_str(),"DijetMass");
      gPad->SaveAs(strOutput);
   }
}

   void
CutFlow::ElectronPt(Int_t n)
{
   Double_t rec_cut_value = (n == 1) ? 20 : (n==2) ? 12 : 7;

   std::string obj("Electron"), var("Pt");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key = sample+delim_+obj+delim_+var;
      char strHist[200];
      Int_t nBins(100);
      Double_t xMin(0), xMax(500);
      snprintf(strHist,200,"%s;%s%selectron p_{T} [GeV];Events / %.1f GeV",
            files_[sample].c_str(), (n==2) ? "sub" : "",(n==1||n==2) ? "leading ":"",
            (xMax-xMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,xMin,xMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);
         if (n == 0)
         {
            for ( auto & pt : *event_vars_[key] ) h->Fill(pt);
         }
         else
         {
            std::vector<Float_t> pts;
            for ( auto & pt : *event_vars_[key] ) pts.push_back(pt);

            Int_t nParticles = pts.size();
            if ( nParticles < n ) continue;

            std::sort(pts.begin(), pts.end(), greater<float>());
            h->Fill(pts[n-1]);
         }

      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      gPad->Update(); 
      printf("%f %f %f %f\n",rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      TLine *l = new TLine(rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      l->SetLineColor(kRed);
      l->SetLineWidth(2);
      l->SetLineStyle(kDashed);
      l->Draw();

      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s_%d.pdf",
            files_[sample].c_str(),obj.c_str(),var.c_str(),n);
      gPad->SaveAs(strOutput);
   }
}

   void
CutFlow::MuonPt(Int_t n)
{
   Double_t rec_cut_value = (n == 1) ? 20 : (n==2) ? 10 : 5;

   std::string obj("Muon"), var("Pt");
   for ( std::string sample : { "signal", "background" } )
   {
      std::string key = sample+delim_+obj+delim_+var;
      char strHist[200];
      Int_t nBins(100);
      Double_t xMin(0), xMax(500);
      snprintf(strHist,200,"%s;%s%smuon p_{T} [GeV];Events / %.1f GeV",
            files_[sample].c_str(), (n==2) ? "sub" : "",(n==1||n==2) ? "leading ":"",
            (xMax-xMin)/nBins);
      TH1D *h = new TH1D("h",strHist,nBins,xMin,xMax);

      fileEventNums_[sample]->open(eventFileName_[sample],ios::in);
      for ( Long64_t entry; *fileEventNums_[sample] >> entry ;)
      {
         readers_[sample]->SetEntry(entry);
         if (n == 0)
         {
            for ( auto & pt : *event_vars_[key] ) h->Fill(pt);
         }
         else
         {
            std::vector<Float_t> pts;
            for ( auto & pt : *event_vars_[key] ) pts.push_back(pt);

            Int_t nParticles = pts.size();
            if ( nParticles < n ) continue;

            std::sort(pts.begin(), pts.end(), greater<float>());
            h->Fill(pts[n-1]);
         }

      }
      fileEventNums_[sample]->close();

      Double_t scale = target_luminosity*cross_section_[sample]/num_events_[sample];
      h->Scale(scale);

      new TCanvas();
      h->Draw("histo");
      gPad->Update(); 
      printf("%f %f %f %f\n",rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      TLine *l = new TLine(rec_cut_value,gPad->GetUymin(),rec_cut_value,gPad->GetUymax());
      l->SetLineColor(kRed);
      l->SetLineWidth(2);
      l->SetLineStyle(kDashed);
      l->Draw();

      char strOutput[400];
      snprintf(strOutput,400,"/home/user108/y4p/tex/selection/graphs/%s_%s_%s_%d.pdf",
            files_[sample].c_str(),obj.c_str(),var.c_str(),n);
      gPad->SaveAs(strOutput);
   }
}


   std::vector<TLorentzVector>
CutFlow::ReconstructVectors(std::string sample, std::string obj_type)
{
   std::string key = sample+delim_+obj_type;
   Int_t nObjs = **object_numbers_[key];
   std::vector<TLorentzVector> rtn;

   std::vector<Double_t> pt, eta, phi, m;
   for ( auto & pt_i : *event_vars_[key+std::string("_Pt")] ) pt.push_back(pt_i);
   for ( auto & eta_i : *event_vars_[key+std::string("_Eta")] ) eta.push_back(eta_i);
   for ( auto & phi_i : *event_vars_[key+std::string("_Phi")] ) phi.push_back(phi_i);
   if ( obj_type == "Jet")
   {
      for ( auto & m_i : *event_vars_[key+std::string("_M")] ) m.push_back(m_i);
   }   
   else
   {
      for ( Int_t i = 0; i < pt.size(); i++ )
         m.push_back( (obj_type == "Muon") ? MUON_MASS : ELECTRON_MASS );
   }

   int a(pt.size()), b(eta.size()), c(phi.size()), d(m.size());
   if (!( a==b && b==c && c==d ))
   {
      printf("size mismatch -- pt: %d, eta: %d, phi: %d, m: %d",a,b,c,d);
   }

   for ( Int_t i = 0; i < a; i++ )
   {
      TLorentzVector obj;
      obj.SetPtEtaPhiM(pt[i],eta[i],phi[i],m[i]);
      rtn.push_back(obj);
   }

   //   for (Int_t iObj; iObj < nObjs; iObj++)
   //   {
   //      TLorentzVector obj;
   //      Double_t pt =     (*event_vars_[key+std::string("_Pt")])[iObj];
   //      Double_t eta =  0;//  (*event_vars_[key+std::string("_Eta")])[iObj]; //these two lines
   //      Double_t phi =    (*event_vars_[key+std::string("_Phi")])[iObj];     //cause seg faults
   //      Double_t m =    0;//  (obj_type == "Jet") ? (*event_vars_[key+std::string("_M")])[iObj]:
   //                                          (obj_type == "Muon") ? MUON_MASS : ELECTRON_MASS;
   //
   //      obj.SetPtEtaPhiM(pt,eta,phi,m);
   //      rtn.push_back(obj);
   //   }

   return rtn;
}

bool
CutFlow::JetIsElectron( std::string sample, Int_t jet_n )
{
   //electron->jet fake removal:
   Double_t delta_r_min = 0.4;
   std::string strJet = "Jet";
   //for the jet at index jet_n in this event, check if it is within delta_r

   TVector3 j3; 
   {
      Float_t j_pt = event_vars_[sample+delim_+strJet+delim_+std::string("Pt")]->At(jet_n);
      Float_t j_eta = event_vars_[sample+delim_+strJet+delim_+std::string("Eta")]->At(jet_n);
      Float_t j_phi = event_vars_[sample+delim_+strJet+delim_+std::string("Phi")]->At(jet_n);
      j3.SetPtEtaPhi(j_pt,j_eta,j_phi);
   }

   std::vector<TLorentzVector> e4s = this->ReconstructVectors( sample, "Electron" );
   for ( TLorentzVector &e4 : e4s )
   {  
      if ( j3.DeltaR( e4.Vect() ) < delta_r_min ) return true;
   }

   return false;
}

   void
CutFlow::ProcessAll(int n_on_shell)
{
   if ( n_on_shell >= 0 )
   {
      PreSelection( n_on_shell );
   }
   sleep(1);
   JetEta(false);
   JetEta(true);
   JetEtaDiff();
   JetPt(1);
   JetPt(2);
   DijetM();
   ElectronPt(1);
   ElectronPt(2);
   MuonPt(1);
   MuonPt(2);
}

   void
CutFlow::FileInit()
{
   for ( std::string sample : {"signal", "background"} )
   {
      fileEventNums_[sample]->open(eventFileName_[sample],ios::out);
      for ( Long64_t entry = 0; entry < num_events_[sample]; entry++ )
         (*fileEventNums_[sample]) << entry << '\t';
      fileEventNums_[sample]->close();
   }  
}
