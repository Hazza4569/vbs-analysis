#include "../utils/constants.h"

using floats = ROOT::VecOps::RVec<float>;
using ints = ROOT::VecOps::RVec<int>;
void rdf(std::string file, bool _save = false)
{
   std::string base_path("/home/user108/y4p/"), ROOT_dir("root/"), output_dir("root_output/"), ROOT_ft(".root"), dat_ft(".dat");

   Int_t n_events; Double_t cross_section;
   std::fstream fs;
   fs.open( (base_path+ROOT_dir+file+dat_ft).c_str() );
   fs >> n_events; //deprecated, better to take from the histogram sizes, as overwritten below
   fs >> cross_section; 
   fs.close();

   ROOT::RDataFrame d("EventTree",base_path+ROOT_dir+file+ROOT_ft);
   
   n_events = d.Histo1D("Muon_n")->GetEntries();//*(d.Take<int>("Event_Count").begin());

   auto get_pair = [](floats &charge, floats &isol, float isolreq) { return (int)min( charge[charge<0 && isol<isolreq].size() , charge[charge>0 && isol<isolreq].size() ); };
   auto vec_4vecs = [](floats &pt, floats &eta, floats &phi, float &m){
      std::vector<TLorentzVector> rtn;
      for (int i = 0; i < pt.size(); i++)
      {
         TLorentzVector v; v.SetPtEtaPhiM( pt.at(i), eta.at(i), phi.at(i), m );
         rtn.emplace_back(v);
      }
      return rtn;
   };
   auto reconstruct_pair = [](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass) {
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
   auto get_pair_y = [&reconstruct_pair](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass)
   { floats rtn; for ( auto &p_i : reconstruct_pair(q,pt,eta,phi,isol,isolreq,mass) ) rtn.emplace_back(p_i.Rapidity()); return rtn; };
   auto combined_var_Zorder = [] (floats &var_ee, floats &var_mumu, floats &m_ee, floats &m_mumu) {
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
   auto combined_var = [] (floats &var_e, floats &var_mu, floats &pt_e, floats &pt_mu) {
      floats rtn;
      std::map<float,float> ordered_vars;
      for ( int i = 0; i < pt_e.size(); i++ ) ordered_vars.insert( make_pair(1e10-pt_e.at(i), var_e.at(i)) );
      for ( int i = 0; i < pt_mu.size(); i++ ) ordered_vars.insert( make_pair(1e10-pt_mu.at(i), var_mu.at(i)) );
      for ( auto &var_i : ordered_vars )
      {
         rtn.emplace_back(var_i.second);
      }
      return rtn;
   };
   auto combined_var_isolated = [] (floats &var_e, floats &var_mu, floats &pt_e, floats &pt_mu, floats &isol_e, floats &isol_mu, float isol_ereq, float isol_mureq) {
      floats rtn;
      std::multimap<float,float> ordered_vars;
      for ( int i = 0; i < pt_e[isol_e<isol_ereq].size(); i++ ) ordered_vars.insert( make_pair(1e10-pt_e[isol_e<isol_ereq].at(i), var_e[isol_e<isol_ereq].at(i)) );
      for ( int i = 0; i < pt_mu[isol_mu<isol_mureq].size(); i++ ) ordered_vars.insert( make_pair(1e10-pt_mu[isol_mu<isol_mureq].at(i), var_mu[isol_mu<isol_mureq].at(i)) );
      for ( auto &var_i : ordered_vars )
      {
         rtn.emplace_back(var_i.second);
      }
      return rtn;
   };
   auto reconstruct_quad = [] (floats &pt, floats &eta, floats &phi, floats &m) {
      if (m.size() < 2) return TLorentzVector(0,0,0,-100);
      TLorentzVector v1, v2;
      v1.SetPtEtaPhiM(pt.at(0), eta.at(0), phi.at(0), m.at(0));
      v2.SetPtEtaPhiM(pt.at(1), eta.at(1), phi.at(1), m.at(1));
      return v1+v2;
   };
   auto get_quad_pt = [&reconstruct_quad](floats &pt, floats &eta, floats &phi, floats &m)
   { return reconstruct_quad(pt,eta,phi,m).Pt(); };
   auto get_quad_eta = [&reconstruct_quad](floats &pt, floats &eta, floats &phi, floats &m)
   { return reconstruct_quad(pt,eta,phi,m).Eta(); };
   auto get_quad_phi = [&reconstruct_quad](floats &pt, floats &eta, floats &phi, floats &m)
   { return reconstruct_quad(pt,eta,phi,m).Phi(); };
   auto get_quad_m = [&reconstruct_quad](floats &pt, floats &eta, floats &phi, floats &m)
   { return reconstruct_quad(pt,eta,phi,m).M(); };
   auto get_quad_y = [&reconstruct_quad](floats &pt, floats &eta, floats &phi, floats &m)
   { return reconstruct_quad(pt,eta,phi,m).Rapidity(); };

   char n_events_str[15],cross_section_str[15];
   snprintf(n_events_str,15,"%d",n_events);
   snprintf(cross_section_str,15,"%f",cross_section);

   //Z definitions
   auto d_new0 = d.Define("Event_Count",n_events_str).Define("Cross_Section",cross_section_str)
      .Define("Electron_Isol_Max","(float)0.35").Define("Muon_Isol_Max","(float)0.35")
      .Define("Electron_Mass","(float)ELECTRON_MASS").Define("Muon_Mass","(float)MUON_MASS")
      .Define("Electron_M","ROOT::VecOps::RVec<float>(Electron_n,ELECTRON_MASS)")
      .Define("Muon_M","ROOT::VecOps::RVec<float>(Muon_n,MUON_MASS)")
      .Define("Electron_Pairs",get_pair,{"Electron_Charge","Electron_Isol","Electron_Isol_Max"})
      .Define("Muon_Pairs",get_pair,{"Muon_Charge","Muon_Isol","Muon_Isol_Max"})
      .Define("Lepton_Pairs","Electron_Pairs+Muon_Pairs")
      .Define("Lepton_n","Electron_n+Muon_n")
      .Define("Lepton_Pt",combined_var,{"Electron_Pt","Muon_Pt","Electron_Pt","Muon_Pt"})
      .Define("Lepton_Eta",combined_var,{"Electron_Eta","Muon_Eta","Electron_Pt","Muon_Pt"})
      .Define("Lepton_Phi",combined_var,{"Electron_Phi","Muon_Phi","Electron_Pt","Muon_Pt"})
      .Define("Lepton_Isol",combined_var,{"Electron_Isol","Muon_Isol","Electron_Pt","Muon_Pt"})
      .Define("Lepton_Charge",combined_var,{"Electron_Charge","Muon_Charge","Electron_Pt","Muon_Pt"})
      .Define("Lepton_M",combined_var,{"Electron_M","Muon_M","Electron_Pt","Muon_Pt"})
      .Define("Isolated_Lepton_Pt",combined_var_isolated,{"Electron_Pt","Muon_Pt","Electron_Pt","Muon_Pt","Electron_Isol","Muon_Isol","Electron_Isol_Max","Muon_Isol_Max"})
      .Define("Isolated_Lepton_Eta",combined_var_isolated,{"Electron_Eta","Muon_Eta","Electron_Pt","Muon_Pt","Electron_Isol","Muon_Isol","Electron_Isol_Max","Muon_Isol_Max"})
      .Define("Isolated_Lepton_Phi",combined_var_isolated,{"Electron_Phi","Muon_Phi","Electron_Pt","Muon_Pt","Electron_Isol","Muon_Isol","Electron_Isol_Max","Muon_Isol_Max"})
      .Define("Isolated_Lepton_Charge",combined_var_isolated,{"Electron_Charge","Muon_Charge","Electron_Pt","Muon_Pt","Electron_Isol","Muon_Isol","Electron_Isol_Max","Muon_Isol_Max"})
      .Define("Isolated_Lepton_M",combined_var_isolated,{"Electron_M","Muon_M","Electron_Pt","Muon_Pt","Electron_Isol","Muon_Isol","Electron_Isol_Max","Muon_Isol_Max"})
      .Define("Isolated_Lepton_n","(int)Isolated_Lepton_Pt.size()")
      .Define("Dielectron_Pt",get_pair_pt,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_Eta",get_pair_eta,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_Phi",get_pair_phi,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_M",get_pair_m,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dielectron_Rapidity",get_pair_y,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass"})
      .Define("Dimuon_Pt",get_pair_pt,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_Eta",get_pair_eta,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_Phi",get_pair_phi,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_M",get_pair_m,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dimuon_Rapidity",get_pair_y,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass"})
      .Define("Dilepton_Pt",combined_var_Zorder,{"Dielectron_Pt","Dimuon_Pt","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_Eta",combined_var_Zorder,{"Dielectron_Eta","Dimuon_Eta","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_Phi",combined_var_Zorder,{"Dielectron_Phi","Dimuon_Phi","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_M",combined_var_Zorder,{"Dielectron_M","Dimuon_M","Dielectron_M","Dimuon_M"})
      .Define("Dilepton_Rapidity",combined_var_Zorder,{"Dielectron_Rapidity","Dimuon_Rapidity","Dielectron_M","Dimuon_M"})
      .Define("Tetralepton_Pt",get_quad_pt,{"Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"})
      .Define("Tetralepton_Eta",get_quad_eta,{"Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"})
      .Define("Tetralepton_Phi",get_quad_phi,{"Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"})
      .Define("Tetralepton_M",get_quad_m,{"Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"})
      .Define("Tetralepton_Rapidity",get_quad_y,{"Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"});

//==================================================================TRUTH LEPTON VALUES============================================================================

   auto d_new00 = d_new0
      .Define("Truth_Electron_Pt","Truth_Lepton_Pt[abs(Truth_Lepton_ID)==11]")
      .Define("Truth_Electron_Eta","Truth_Lepton_Eta[abs(Truth_Lepton_ID)==11]")
      .Define("Truth_Electron_Phi","Truth_Lepton_Phi[abs(Truth_Lepton_ID)==11]")
      .Define("Truth_Electron_n","(int)Truth_Electron_Pt.size()")
      .Define("Truth_Electron_M","ROOT::VecOps::RVec<float>(Truth_Electron_n,ELECTRON_MASS)")
      .Define("Truth_Electron_Charge","auto id = Truth_Lepton_ID[abs(Truth_Lepton_ID)==11]; ROOT::VecOps::RVec<float> rtn;"
                                    "for (int i=0; i < Truth_Electron_n; i++) rtn.push_back(id.at(i)/abs(id.at(i))); return rtn;")
      .Define("Truth_Electron_Isol","ROOT::VecOps::RVec<float> rtn; for (int i = 0; i < Truth_Electron_n; i++) rtn.emplace_back(0); return rtn;")

      .Define("Truth_Muon_Pt","Truth_Lepton_Pt[abs(Truth_Lepton_ID)==13]")
      .Define("Truth_Muon_Eta","Truth_Lepton_Eta[abs(Truth_Lepton_ID)==13]")
      .Define("Truth_Muon_Phi","Truth_Lepton_Phi[abs(Truth_Lepton_ID)==13]")
      .Define("Truth_Muon_n","(int)Truth_Muon_Pt.size()")
      .Define("Truth_Muon_M","ROOT::VecOps::RVec<float>(Truth_Muon_n,MUON_MASS)")
      .Define("Truth_Muon_Charge","auto id = Truth_Lepton_ID[abs(Truth_Lepton_ID)==13]; ROOT::VecOps::RVec<float> rtn;"
                                    "for (int i=0; i < Truth_Muon_n; i++) rtn.push_back(id.at(i)/abs(id.at(i))); return rtn;")
      .Define("Truth_Muon_Isol","ROOT::VecOps::RVec<float> rtn; for (int i = 0; i < Truth_Muon_n; i++) rtn.emplace_back(0); return rtn;")

      .Define("Truth_Electron_Isol_Max","Electron_Isol_Max").Define("Truth_Muon_Isol_Max","Muon_Isol_Max")
      .Define("Truth_Electron_Mass","(float)ELECTRON_MASS").Define("Truth_Muon_Mass","(float)MUON_MASS")
      .Define("Truth_Electron_Pairs",get_pair,{"Truth_Electron_Charge","Truth_Electron_Isol","Truth_Electron_Isol_Max"})
      .Define("Truth_Muon_Pairs",get_pair,{"Truth_Muon_Charge","Truth_Muon_Isol","Truth_Muon_Isol_Max"})
      .Define("Truth_Lepton_Pairs","Truth_Electron_Pairs+Truth_Muon_Pairs")
      .Define("Truth_Lepton_Isol",combined_var,{"Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Pt","Truth_Muon_Pt"})
      .Define("Truth_Lepton_Charge",combined_var,{"Truth_Electron_Charge","Truth_Muon_Charge","Truth_Electron_Pt","Truth_Muon_Pt"})
      .Define("Truth_Lepton_M",combined_var,{"Truth_Electron_M","Truth_Muon_M","Truth_Electron_Pt","Truth_Muon_Pt"})
      .Define("Truth_Isolated_Lepton_Pt",combined_var_isolated,
         {"Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Isol_Max","Truth_Muon_Isol_Max"})
      .Define("Truth_Isolated_Lepton_Eta",combined_var_isolated,
         {"Truth_Electron_Eta","Truth_Muon_Eta","Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Isol_Max","Truth_Muon_Isol_Max"})
      .Define("Truth_Isolated_Lepton_Phi",combined_var_isolated,
         {"Truth_Electron_Phi","Truth_Muon_Phi","Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Isol_Max","Truth_Muon_Isol_Max"})
      .Define("Truth_Isolated_Lepton_Charge",combined_var_isolated,
         {"Truth_Electron_Charge","Truth_Muon_Charge","Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Isol_Max","Truth_Muon_Isol_Max"})
      .Define("Truth_Isolated_Lepton_M",combined_var_isolated,
         {"Truth_Electron_M","Truth_Muon_M","Truth_Electron_Pt","Truth_Muon_Pt","Truth_Electron_Isol","Truth_Muon_Isol","Truth_Electron_Isol_Max","Truth_Muon_Isol_Max"})
      .Define("Truth_Isolated_Lepton_n","(int)Truth_Isolated_Lepton_Pt.size()")
      .Define("Truth_Dielectron_Pt",get_pair_pt,{"Truth_Electron_Charge","Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Isol","Truth_Electron_Isol_Max","Truth_Electron_Mass"})
      .Define("Truth_Dielectron_Eta",get_pair_eta,{"Truth_Electron_Charge","Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Isol","Truth_Electron_Isol_Max","Truth_Electron_Mass"})
      .Define("Truth_Dielectron_Phi",get_pair_phi,{"Truth_Electron_Charge","Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Isol","Truth_Electron_Isol_Max","Truth_Electron_Mass"})
      .Define("Truth_Dielectron_M",get_pair_m,{"Truth_Electron_Charge","Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Isol","Truth_Electron_Isol_Max","Truth_Electron_Mass"})
      .Define("Truth_Dielectron_Rapidity",get_pair_y,{"Truth_Electron_Charge","Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Isol","Truth_Electron_Isol_Max","Truth_Electron_Mass"})
      .Define("Truth_Dimuon_Pt",get_pair_pt,{"Truth_Muon_Charge","Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Isol","Truth_Muon_Isol_Max","Truth_Muon_Mass"})
      .Define("Truth_Dimuon_Eta",get_pair_eta,{"Truth_Muon_Charge","Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Isol","Truth_Muon_Isol_Max","Truth_Muon_Mass"})
      .Define("Truth_Dimuon_Phi",get_pair_phi,{"Truth_Muon_Charge","Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Isol","Truth_Muon_Isol_Max","Truth_Muon_Mass"})
      .Define("Truth_Dimuon_M",get_pair_m,{"Truth_Muon_Charge","Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Isol","Truth_Muon_Isol_Max","Truth_Muon_Mass"})
      .Define("Truth_Dimuon_Rapidity",get_pair_y,{"Truth_Muon_Charge","Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Isol","Truth_Muon_Isol_Max","Truth_Muon_Mass"})
      .Define("Truth_Dilepton_Pt",combined_var_Zorder,{"Truth_Dielectron_Pt","Truth_Dimuon_Pt","Truth_Dielectron_M","Truth_Dimuon_M"})
      .Define("Truth_Dilepton_Eta",combined_var_Zorder,{"Truth_Dielectron_Eta","Truth_Dimuon_Eta","Truth_Dielectron_M","Truth_Dimuon_M"})
      .Define("Truth_Dilepton_Phi",combined_var_Zorder,{"Truth_Dielectron_Phi","Truth_Dimuon_Phi","Truth_Dielectron_M","Truth_Dimuon_M"})
      .Define("Truth_Dilepton_M",combined_var_Zorder,{"Truth_Dielectron_M","Truth_Dimuon_M","Truth_Dielectron_M","Truth_Dimuon_M"})
      .Define("Truth_Dilepton_Rapidity",combined_var_Zorder,{"Truth_Dielectron_Rapidity","Truth_Dimuon_Rapidity","Truth_Dielectron_M","Truth_Dimuon_M"})
      .Define("Truth_Tetralepton_Pt",get_quad_pt,{"Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
      .Define("Truth_Tetralepton_Eta",get_quad_eta,{"Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
      .Define("Truth_Tetralepton_Phi",get_quad_phi,{"Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
      .Define("Truth_Tetralepton_M",get_quad_m,{"Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
      .Define("Truth_Tetralepton_Rapidity",get_quad_y,{"Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"});

//=================================================================================================================================================================

   auto get_delta_r = [&vec_4vecs](floats &pt, floats &eta, floats &phi, floats &m, floats &e_pt, floats &e_eta, floats &e_phi, float &e_m){
      floats rtn;
      for (int i = 0; i < pt.size(); i++)
      {
         TLorentzVector j; j.SetPtEtaPhiM( pt.at(i), eta.at(i), phi.at(i), m.at(i) );
         double dr_min = INT_MAX;
         for ( auto &e_i : vec_4vecs(e_pt,e_eta,e_phi,e_m) )
         {
            double dr = j.DeltaR(e_i);
            dr_min = (dr < dr_min) ? dr : dr_min;
         }
         rtn.emplace_back(dr_min);
      }
      return rtn;
   };
   auto n_good_jets = [](ints &flav, int flav_max, floats &dR, float dR_min) { return (int) flav[dR>=dR_min && flav<flav_max].size(); };
   auto jet_y = [] (floats &pt, floats &eta, floats &phi, floats &m) {
      floats rtn;
      for (int i = 0; i < pt.size(); i++)
      {
         TLorentzVector v; v.SetPtEtaPhiM(pt.at(i),eta.at(i),phi.at(i),m.at(i));
         rtn.emplace_back(v.Rapidity());
      }
      return rtn;
   };
   auto reconstruct_dijet = [] (floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min) {
      if (m[dR>=dR_min && flav<flav_max].size() < 2) return TLorentzVector(0,0,0,-100);
      TLorentzVector v1, v2;
      v1.SetPtEtaPhiM(pt[dR>=dR_min && flav<flav_max].at(0), eta[dR>=dR_min && flav<flav_max].at(0), phi[dR>=dR_min && flav<flav_max].at(0), m[dR>=dR_min && flav<flav_max].at(0));
      v2.SetPtEtaPhiM(pt[dR>=dR_min && flav<flav_max].at(1), eta[dR>=dR_min && flav<flav_max].at(1), phi[dR>=dR_min && flav<flav_max].at(1), m[dR>=dR_min && flav<flav_max].at(1));
      return v1+v2;
   };
   auto get_dijet_pt = [&reconstruct_dijet](floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min)
   { return reconstruct_dijet(pt,eta,phi,m,flav,flav_max,dR,dR_min).Pt(); };
   auto get_dijet_eta = [&reconstruct_dijet](floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min)
   { return reconstruct_dijet(pt,eta,phi,m,flav,flav_max,dR,dR_min).Eta(); };
   auto get_dijet_phi = [&reconstruct_dijet](floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min)
   { return reconstruct_dijet(pt,eta,phi,m,flav,flav_max,dR,dR_min).Phi(); };
   auto get_dijet_m = [&reconstruct_dijet](floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min)
   { return reconstruct_dijet(pt,eta,phi,m,flav,flav_max,dR,dR_min).M(); };
   auto get_dijet_y = [&reconstruct_dijet](floats &pt, floats &eta, floats &phi, floats &m, ints &flav, int flav_max, floats &dR, float dR_min)
   { return reconstruct_dijet(pt,eta,phi,m,flav,flav_max,dR,dR_min).Rapidity(); };
   auto get_eta_diff = [](floats &eta, ints &flav, int flav_max, floats &dR, float dR_min)
   { return (eta[dR>=dR_min && flav<flav_max].size()<2) ? -1 : fabs(eta[dR>=dR_min && flav<flav_max].at(1)-eta[dR>=dR_min && flav<flav_max].at(0)); };

   //j Deinitions
   auto d_new = d_new00.Define("Jet_DeltaR_Min","(float)0.3").Define("Jet_Flav_Max","(int)100")
                      .Define("Jet_DeltaR",get_delta_r,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Mass"})
                      .Define("Jet_Good_n",n_good_jets,{"Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Jet_Rapidity",jet_y,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M"})
                      .Define("Dijet_Pt",get_dijet_pt,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Dijet_Eta",get_dijet_eta,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Dijet_Phi",get_dijet_phi,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Dijet_M",get_dijet_m,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Dijet_Rapidity",get_dijet_y,{"Jet_Pt","Jet_Eta","Jet_Phi","Jet_M","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Jet12_Eta_Diff",get_eta_diff,{"Jet_Eta","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"})
                      .Define("Jet12_Rapidity_Diff",get_eta_diff,{"Jet_Rapidity","Jet_Flav","Jet_Flav_Max","Jet_DeltaR","Jet_DeltaR_Min"});

   if (_save) d_new.Snapshot("EventTree",base_path+output_dir+file+std::string("_o")+ROOT_ft);
   printf("saved\n");

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
   
   auto d_twojet = d_new.Filter("Jet_Good_n > 1", "two_good_jets");
   d_new.Filter("Jet_n > 1", "two_jets");   

//   FILE *of = fopen((base_path+output_dir+file+std::string("_cutflow.dat")).c_str(),"w");
//
//   int i = 0;
//   for (auto &&cut : d.Report())
//   {
//      if ( cut.GetName() == "4e_optimal" )
//      {
//         i++;
//         if ( i == 1 ) fprintf( of, "\nRequiring any two pairs:\n" );
//         else if ( i == 2 ) fprintf( of, "\nRequiring two pairs on the Z mass shell (10 GeV margin):\n" );
//      }
//      fprintf( of, "%-22spass=%-10llu\t\tall=%-10llu\t--\t%.3f %%\n", (cut.GetName()+std::string(":")).c_str(), cut.GetPass(), cut.GetAll(), cut.GetEff() );
//   } 
}
