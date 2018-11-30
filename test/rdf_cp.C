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
   auto get_mz = [](floats &q, floats &pt, floats &eta, floats &phi, floats &isol, float isolreq, float mass, int n)
   {
      floats m_z; 
      for (int i = 0; i < q[q==-1 && isol<isolreq].size(); i++)
      {
         TLorentzVector neg, pos;
         neg.SetPtEtaPhiM(pt[q==-1 && isol<isolreq].at(i),eta[q==-1 && isol<isolreq].at(i),phi[q==-1 && isol<isolreq].at(i),mass);
         for (int j = 0; j < q[q==1 && isol<isolreq].size(); j++)
         {
            pos.SetPtEtaPhiM(pt[q==1 && isol<isolreq].at(j),eta[q==1 && isol<isolreq].at(j),phi[q==1 && isol<isolreq].at(j),mass);
            m_z.emplace_back( (float) (pos+neg).M() );
         }
      }
      //take only two best z candidates
      if (m_z.size() < n) printf("size error -- n_pairs: %d, individual vector size: %lu\n", n, m_z.size());
      std::sort(m_z.begin(),m_z.end(), [](float lhs, float rhs){ return fabs(lhs-Z_MASS) < fabs(rhs-Z_MASS); } );
      if (m_z.size() > 2) m_z.erase(m_z.begin()+2,m_z.end());
      return m_z;
   };
   auto combine_mll = [](floats &m_ll, floats &m_ll2, int n)
   {
      std::vector<float> tmp;
      int init_size = m_ll.size() + m_ll2.size();
      if ( init_size < n ) printf("size error -- n_pairs: %d, combined vector sizes: %d\n",n,init_size);
      tmp.reserve( init_size );
      tmp.insert( tmp.end(), m_ll.begin(), m_ll.end() );
      tmp.insert( tmp.end(), m_ll2.begin(), m_ll2.end() );
      floats m_z = tmp;
      std::sort(m_z.begin(),m_z.end(), [](float lhs, float rhs){ return fabs(lhs-Z_MASS) < fabs(rhs-Z_MASS); } );
      if (m_z.size() > 2) m_z.erase(m_z.begin()+2,m_z.end());
      if ( (init_size >= 2) && (m_z.size() < 2) ) printf("size error -- input size %d, output size %lu",init_size,m_z.size());
      return m_z;
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
      .Define("m_ee",get_mz,{"Electron_Charge","Electron_Pt","Electron_Eta","Electron_Phi","Electron_Isol","Electron_Isol_Max","Electron_Mass","Electron_Pairs"})
      .Define("m_mumu",get_mz,{"Muon_Charge","Muon_Pt","Muon_Eta","Muon_Phi","Muon_Isol","Muon_Isol_Max","Muon_Mass","Muon_Pairs"})
      .Define("Lepton_Pairs","Electron_Pairs+Muon_Pairs").Define("m_ll",combine_mll,{"m_ee","m_mumu","Lepton_Pairs"});
   //j Deinitions
   auto d_new = d_new0
   if (_save) d_new.Snapshot("EventTree",base_path+output_dir+file+std::string("_o")+ROOT_ft);

   auto d_twopair = d_new.Filter("Lepton_Pairs > 1","two_pairs");
   auto d_twoshell = d_twopair.Filter("fabs(m_ll.at(1)-Z_MASS) < 10","two_on_shell");

   d_twopair.Filter("Electron_Pairs > 1 && ( Muon_Pairs == 0 || fabs(m_ee.at(1)-Z_MASS) < fabs(m_mumu.at(0)-Z_MASS) )","4e_optimal");
   d_twopair.Filter("Muon_Pairs > 1 && ( Electron_Pairs == 0 || fabs(m_mumu.at(1)-Z_MASS) < fabs(m_ee.at(0)-Z_MASS) )","4mu_optimal");
   d_twopair.Filter("(Muon_Pairs == 1 &&  Electron_Pairs == 1) ||\
         ( Muon_Pairs > 1 && Electron_Pairs > 1 && fabs(m_mumu.at(1)-Z_MASS) > fabs(m_ee.at(0)) && fabs(m_ee.at(1)-Z_MASS) > fabs(m_mumu.at(0)-Z_MASS) ) ||\
         ( Muon_Pairs > 1 && Electron_Pairs == 1 && fabs(m_mumu.at(1)-Z_MASS) > fabs(m_ee.at(0)) ) ||\
         ( Electron_Pairs > 1 && Muon_Pairs == 1 && fabs(m_ee.at(1)-Z_MASS) > fabs(m_mumu.at(0)) )","2e2mu_optimal");

   d_twopair.Filter("Electron_Pairs == 2 && Muon_Pairs == 0","only_4e_present");
   d_twopair.Filter("Muon_Pairs == 2 && Electron_Pairs == 0","only_4mu_present");
   d_twopair.Filter("Electron_Pairs == 1 && Muon_Pairs == 1","only_2e2mu_present");

   d_twopair.Filter("Muon_n > 0","has_lead_mu").Filter("Muon_n > 1","has_sub_mu");
   d_twopair.Filter("Electron_n > 0","has_lead_e").Filter("Electron_n > 1","has_sub_e");

   d_twoshell.Filter("Electron_Pairs > 1 && ( Muon_Pairs == 0 || fabs(m_ee.at(1)-Z_MASS) < fabs(m_mumu.at(0)-Z_MASS) )","4e_optimal");
   d_twoshell.Filter("Muon_Pairs > 1 && ( Electron_Pairs == 0 || fabs(m_mumu.at(1)-Z_MASS) < fabs(m_ee.at(0)-Z_MASS) )","4mu_optimal");
   d_twoshell.Filter("(Muon_Pairs == 1 &&  Electron_Pairs == 1) ||\
         ( Muon_Pairs > 1 && Electron_Pairs > 1 && fabs(m_mumu.at(1)-Z_MASS) > fabs(m_ee.at(0)) && fabs(m_ee.at(1)-Z_MASS) > fabs(m_mumu.at(0)-Z_MASS) ) ||\
         ( Muon_Pairs > 1 && Electron_Pairs == 1 && fabs(m_mumu.at(1)-Z_MASS) > fabs(m_ee.at(0)) ) ||\
         ( Electron_Pairs > 1 && Muon_Pairs == 1 && fabs(m_ee.at(1)-Z_MASS) > fabs(m_mumu.at(0)) )","2e2mu_optimal");

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
