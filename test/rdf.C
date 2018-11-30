#include "../utils/constants.h"

using floats = ROOT::VecOps::RVec<float>;
void rdf(std::string file, bool _save = false)
{
   std::string base_path("/home/user108/y4p/"), ROOT_dir("root/"), output_dir("root_output/"), ROOT_ft(".root"), dat_ft(".dat");

   Int_t n_events; Double_t cross_section;
   std::fstream fs;
   fs.open( (base_path+ROOT_dir+file+dat_ft).c_str() );
   fs >> n_events;
   fs >> cross_section; 
   fs.close();

   ROOT::RDataFrame d("EventTree",base_path+ROOT_dir+file+ROOT_ft);

   auto get_pair = [](floats &charge, floats &isol, float isolreq) { return (int)min(charge[charge==-1 && isol<isolreq].size(), charge[charge==1 && isol<isolreq].size()); };
   auto reconstructi_pair = [](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
   {
      std::vector<TLorentzVector> p_z; 
      for (int i = 0; i < q[q==-1 && isol<isolreq].size(); i++)
      {
         TLorentzVector neg, pos;
         neg.SetPtEtaPhiM(pt[q==-1 && isol<isolreq].at(i),eta[q==-1 && isol<isolreq].at(i),phi[q==-1 && isol<isolreq].at(i),mass);
         for (int j = 0; j < q[q==1 && isol<isolreq].size(); j++)
         {
            pos.SetPtEtaPhiM(pt[q==1 && isol<isolreq].at(j),eta[q==1 && isol<isolreq].at(j),phi[q==1 && isol<isolreq].at(j),mass);
            p_z.emplace_back( pos + neg );
         }
      }
      //take only two best z candidates
      std::sort(p_z.begin(),p_z.end(), [](TLorentzVector lhs, TLorentzVector rhs){ return fabs(lhs.M()-Z_MASS) < fabs(rhs.M()-Z_MASS); } );
      if (p_z.size() > 2) p_z.erase(p_z.begin()+2,p_z.end());
      return p_z;
   };
   auto get_pair_pt = [&reconstruct_pair](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
      { floats rtn; for ( auto &p_i : reconstruct_pair(q,pt,eta,phi,isol,isolreq,mass) ) rtn.emplace_back(p_i.Pt()); return rtn; };
   auto get_pair_eta = [&reconstruct_pair](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
      { floats rtn; for ( auto &p_i : reconstruct_pair(q,pt,eta,phi,isol,isolreq,mass) ) rtn.emplace_back(p_i.Eta()); return rtn; };
   auto get_pair_phi = [&reconstruct_pair](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
      { floats rtn; for ( auto &p_i : reconstruct_pair(q,pt,eta,phi,isol,isolreq,mass) ) rtn.emplace_back(p_i.Phi()); return rtn; };
   auto get_pair_m = [&reconstruct_pair](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
      { floats rtn; for ( auto &p_i : reconstruct_pair(q,pt,eta,phi,isol,isolreq,mass) ) rtn.emplace_back(p_i.M()); return rtn; };
   auto combined_var = [] (floats &var_ee, floats &var_mumu, floats &m_ee, floats &m_mumu)
      {
         floats rtn;
         std::map<float,float> ordered_vars;
         for ( int i = 0; i < m_ee.size(); i++ ) ordered_vars.insert( make_pair(fabs(m_ee.at(i)-Z_MASS), var_ee.at(i)) );
         for ( int i = 0; i < m_mumu.size(); i++ ) ordered_vars.insert( make_pair(fabs(m_mumu.at(i)-Z_MASS), var_mumu.at(i)) );
         for ( auto &var_i : ordered_vars )
         {
            if (rtn.size() == 2) break;
            rtn.emplace_back(var_i.second);
         }
         return rtn;
      };

   char n_events_str[15],cross_section_str[15];
   snprintf(n_events_str,15,"%d",n_events);
   snprintf(cross_section_str,15,"%f",cross_section);
   
   //Z definitions
   auto d_new0 = d.Define("Event_Count",n_events_str).Define("Cross_Section",cross_section_str)
      .Define("Electron_Isol_Max","(float)0.05").Define("Muon_Isol_Max","(float)0.05")
      .Define("Electron_Mass","(float)ELECTRON_MASS").Define("Muon_Mass","(float)MUON_MASS")
      .Define("Electron_Pairs",get_pair,{"Electron_Charge","Electron_Isol","Electron_Isol_Max"})
      .Define("Muon_Pairs",get_pair,{"Muon_Charge","Muon_Isol","Muon_Isol_Max"})
      .Define("Lepton_Pairs","Electron_Pairs+Muon_Pairs")
      .Define("Dielectron_Pt",get_pair_pt,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_Eta",get_pair_eta,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_Phi",get_pair_phi,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_M",get_pair_m,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dimuon_Pt",get_pair_pt,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_Eta",get_pair_eta,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_Phi",get_pair_phi,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_M",get_pair_m,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dilepton_Pt",combined_var,{"Dielectron_Pt","Dimuon_Pt","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_Eta",combined_var,{"Dielectron_Eta","Dimuon_Eta","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_Phi",combined_var,{"Dielectron_Phi","Dimuon_Phi","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_M",combined_var,{"Dielectron_M","Dimuon_M","Dielectron_M","Dimuon_M"});

   //j Deinitions
   auto d_new = d_new0.Define("test","1.");
   if (_save) d_new.Snapshot("EventTree",base_path+output_dir+file+std::string("_o")+ROOT_ft);

   auto d_twopair = d_new.Filter("Lepton_Pairs > 1","two_pairs");
   auto d_twoshell = d_twopair.Filter("fabs(Dilepton_M.at(1)-Z_MASS) < 10","two_on_shell");

   d_twopair.Filter("Electron_Pairs > 1 && ( Muon_Pairs == 0 || fabs(Dielectron_M.at(1)-Z_MASS) < fabs(Dimuon_M.at(0)-Z_MASS) )","4e_optimal");
   d_twopair.Filter("Muon_Pairs > 1 && ( Electron_Pairs == 0 || fabs(Dimuon_M.at(1)-Z_MASS) < fabs(Dielectron_M.at(0)-Z_MASS) )","4mu_optimal");
   d_twopair.Filter("(Muon_Pairs == 1 &&  Electron_Pairs == 1) ||\
         ( Muon_Pairs > 1 && Electron_Pairs > 1 && fabs(Dimuon_M.at(1)-Z_MASS) > fabs(Dielectron_M.at(0)-Z_MASS) && fabs(Dielectron_M.at(1)-Z_MASS) > fabs(Dimuon_M.at(0)-Z_MASS) ) ||\
         ( Muon_Pairs > 1 && Electron_Pairs == 1 && fabs(Dimuon_M.at(1)-Z_MASS) > fabs(Dielectron_M.at(0))-Z_MASS ) ||\
         ( Electron_Pairs > 1 && Muon_Pairs == 1 && fabs(Dielectron_M.at(1)-Z_MASS) > fabs(Dimuon_M.at(0)-Z_MASS) )","2e2mu_optimal");

   d_twopair.Filter("Electron_Pairs == 2 && Muon_Pairs == 0","only_4e_present");
   d_twopair.Filter("Muon_Pairs == 2 && Electron_Pairs == 0","only_4mu_present");
   d_twopair.Filter("Electron_Pairs == 1 && Muon_Pairs == 1","only_2e2mu_present");

   d_twopair.Filter("Muon_n > 0","has_lead_mu").Filter("Muon_n > 1","has_sub_mu");
   d_twopair.Filter("Electron_n > 0","has_lead_e").Filter("Electron_n > 1","has_sub_e");

   d_twoshell.Filter("Electron_Pairs > 1 && ( Muon_Pairs == 0 || fabs(Dielectron_M.at(1)-Z_MASS) < fabs(Dimuon_M.at(0)-Z_MASS) )","4e_optimal");
   d_twoshell.Filter("Muon_Pairs > 1 && ( Electron_Pairs == 0 || fabs(Dimuon_M.at(1)-Z_MASS) < fabs(Dielectron_M.at(0)-Z_MASS) )","4mu_optimal");
   d_twoshell.Filter("(Muon_Pairs == 1 &&  Electron_Pairs == 1) ||\
         ( Muon_Pairs > 1 && Electron_Pairs > 1 && fabs(Dimuon_M.at(1)-Z_MASS) > fabs(Dielectron_M.at(0)) && fabs(Dielectron_M.at(1)-Z_MASS) > fabs(Dimuon_M.at(0)-Z_MASS) ) ||\
         ( Muon_Pairs > 1 && Electron_Pairs == 1 && fabs(Dimuon_M.at(1)-Z_MASS) > fabs(Dielectron_M.at(0)) ) ||\
         ( Electron_Pairs > 1 && Muon_Pairs == 1 && fabs(Dielectron_M.at(1)-Z_MASS) > fabs(Dimuon_M.at(0)) )","2e2mu_optimal");

   d_twoshell.Filter("Electron_Pairs == 2 && Muon_Pairs == 0","only_4e_present");
   d_twoshell.Filter("Muon_Pairs == 2 && Electron_Pairs == 0","only_4mu_present");
   d_twoshell.Filter("Electron_Pairs == 1 && Muon_Pairs == 1","only_2e2mu_present");

   d_twoshell.Filter("Muon_n > 0","has_lead_mu").Filter("Muon_n > 1","has_sub_mu");
   d_twoshell.Filter("Electron_n > 0","has_lead_e").Filter("Electron_n > 1","has_sub_e");

   FILE *of = fopen((base_path+output_dir+file+std::string("_cutflow.dat")).c_str(),"w");

   int i = 0;
   for (auto &&cut : d.Report())
   {
      if ( cut.GetName() == "4e_optimal" )
      {
         i++;
         if ( i == 1 ) fprintf( of, "\nRequiring any two pairs:\n" );
         else if ( i == 2 ) fprintf( of, "\nRequiring two pairs on the Z mass shell (10 GeV margin):\n" );
      }
      fprintf( of, "%-22spass=%-10llu\t\tall=%-10llu\t--\t%.3f %%\n", (cut.GetName()+std::string(":")).c_str(), cut.GetPass(), cut.GetAll(), cut.GetEff() );
   } 
}
