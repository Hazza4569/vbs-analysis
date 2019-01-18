#include "../utils/constants.h"
using floats = ROOT::VecOps::RVec<float>;
void selection_truth(string date)
{
   FILE *of = fopen((string("/home/user108/y4p/cutflows/truth_selection_cutflow_")+date+string(".dat")).c_str(),"w");
   new TCanvas();
   gStyle->SetOptStat(1111111);
   for (std::string file : {"ZZjj_ATLAS_500K","inclusive_ATLAS_500K"})
   {
      //setup:
      double target_luminosity = 35.9; //As CMS analysis
      std::string pre("/home/user108/y4p/root_output/"),post("_o.root");

      ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
      int n_events = *(d.Take<int>("Event_Count").begin());
      double cross_section = *(d.Take<double>("Cross_Section").begin());
      Double_t scale = target_luminosity*cross_section/n_events;

      int bins(1); double rmin(0), rmax(1); char histstr[100];
      auto histname = [&bins,&rmin,&rmax,&file](string xlabel, string units){
         char rtn[500]; bool nu=(units=="");
         snprintf(rtn,500,"%s (truth level);%s %s%s%s; Events / %g %s",file.c_str(),xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
               (rmax-rmin)/bins,units.c_str());
         return rtn;
      };
      int cutnum = 0;
      auto savefile = [&](string var){
         char rtn[500];
         snprintf(rtn,500,"/home/user108/y4p/graph_logs/%s/selection1_truth_%02d_%s_%s.pdf",date.c_str(),cutnum++,var.c_str(),file.c_str());
         return rtn;
      };
      TH1F *h;

      double muon_leadpt(20), electron_leadpt(20);
      double muon_subpt(10), electron_subpt(12);
      double z_shell_proximity(30);
      double jet_subpt(40);
      double jet_etamax(4.5);
      double dijet_mass_loose(100), dijet_mass_tight(400);
      double jet_eta_diff(2.4);
      double ATLAS_e_etamax(2.4), ATLAS_mu_etamax(2.4);

      auto d_new = d.Define("Truth_Lead_Lepton_Pt","Truth_Isolated_Lepton_Pt.size() < 1 ? -1 : Truth_Isolated_Lepton_Pt.at(0)")
         .Define("Truth_Sublead_Lepton_Pt","Truth_Isolated_Lepton_Pt.size() < 2 ? -1 : Truth_Isolated_Lepton_Pt.at(1)")
         .Define("Sublead_Jet_Pt","Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].size() < 2 ? -1 : Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1)")
         .Define("Truth_Leptons_Contained_ATLAS",
                 string("int rtn(0); for (int i=0; i < Truth_Lepton_n; i++) if ( abs(Truth_Lepton_Eta.at(i)) < ((abs(Truth_Lepton_ID.at(i)) == 13) ? ")
                 + to_string(ATLAS_mu_etamax) + string(" : ") + to_string(ATLAS_e_etamax) + string(") ) rtn++; return rtn;"))
         .Define("Truth_Leptons_Contained_CMS",
                 "int rtn(0); for (int i=0; i < Truth_Lepton_n; i++) if ( abs(Truth_Lepton_Eta.at(i)) < ((abs(Truth_Lepton_ID.at(i)) == 13) ? 2.4 : 5.0) ) rtn++; return rtn;");

//---------------------------------------hacky rdatafame filter workaround (needing RNodes in future ROOT versions-----------------------------------
//         [&muon_subpt,&electron_subpt](floats &mu_pt, floats &e_pt, floats &mu_isol, floats &e_isol, float mu_max, float e_max){
//            int mu_size( mu_pt[mu_isol<mu_max].size() ), e_size( e_pt[e_isol<e_max].size() );
//            double mu_lead( (mu_size > 0) ? mu_pt[mu_isol<mu_max].at(0) : 0 ), e_lead( (e_size > 0) ? e_pt[e_isol<e_max].at(0) : 0 );
//            double lead, H0_sub, H0_subreq, H1_subreq, H1_isol_max; floats *H1, *H1_isol;
//            //H0 is that the subleading lepton will be of opposite type to the leading.
//            if ( e_lead > mu_lead ) {
//               lead = e_lead; H0_sub = mu_lead; H0_subreq = muon_subpt; H1_subreq = electron_subpt;
//               H1 = &e_pt; H1_isol = &e_isol; H1_isol_max = e_max;
//            } else {
//               lead = mu_lead; H0_sub = e_lead; H0_subreq = electron_subpt; H1_subreq = muon_subpt;
//               H1 = &mu_pt; H1_isol = &mu_isol; H1_isol_max = mu_max;
//            }
//            for (auto &h1pt : (*H1)[*H1_isol<H1_isol_max]) if (h1pt > H0_sub && h1pt != lead) return (h1pt > H1_subreq);
//            return (H0_sub > H0_subreq);
//         },
//         {"Muon_Pt","Electron_Pt","Muon_Isol","Electron_Isol","Muon_Isol_Max","Electron_Isol_Max"},

      char subptfuncstring[2000];
      snprintf(subptfuncstring,2000,
            "double muon_subpt(%f), electron_subpt(%f);"
            "int mu_size( Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].size() ), e_size( Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].size() );"
            "double mu_lead( (mu_size > 0) ? Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].at(0) : 0 ), e_lead( (e_size > 0) ? Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].at(0) : 0 );"
            "double lead, H0_sub, H0_subreq, H1_subreq, H1_isol_max; ROOT::VecOps::RVec<float> *H1, *H1_isol;"
            //H0 is that the subleading lepton will be of opposite type to the leading.
            "if ( e_lead > mu_lead ) {"
            "   lead = e_lead; H0_sub = mu_lead; H0_subreq = muon_subpt; H1_subreq = electron_subpt;"
            "   H1 = &Truth_Electron_Pt; H1_isol = &Truth_Electron_Isol; H1_isol_max = Truth_Electron_Isol_Max;"
            "} else {"
            "   lead = mu_lead; H0_sub = e_lead; H0_subreq = electron_subpt; H1_subreq = muon_subpt;"
            "   H1 = &Truth_Muon_Pt; H1_isol = &Truth_Muon_Isol; H1_isol_max = Truth_Muon_Isol_Max;"
            "}"
            "for (auto &h1pt : (*H1)[*H1_isol<H1_isol_max]) if (h1pt > H0_sub && h1pt != lead) return (h1pt > H1_subreq);"
            "return (H0_sub > H0_subreq);",
         muon_subpt, electron_subpt); 
//----------------------------------------------------------------------------------------------------------------------------------------------------

      auto d_start = d_new.Filter("true");
      auto *d_curr = &d_start;
      cutnum = 0;

      //Basic preselection:
      auto dtrig = d_curr->Filter(
            //SINGLE LEPTONS:
            "(Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].size() > 0 &&\
            Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].at(0) > 27) ||"                                                         //single isolated muon, p_T > 27 GeV
            "(Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].size() > 0 &&\
            Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].at(0) > 27) ||"                                             //single isolated (tight?) electron, p_T > 27 GeV
            "(Truth_Muon_n > 0 && Truth_Muon_Pt.at(0) > 52) ||"                                                                  //single muon, p_T > 52 GeV
            "(Truth_Electron_n > 0 && Truth_Electron_Pt.at(0) > 61)"                                                             //single electron, p_T > 61 GeV
            , "trigger simulation|ATLAS single lepton trigger requirements");
      d_curr = &dtrig;


      rmin = 0; rmax = 600; bins = 80; 
      h = (TH1F*)d_curr->Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Truth_Lead_Lepton_Pt")->Clone("Truth_Lead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Lead_Lepton_Pt"));
      auto d_leadlep = d_curr->Filter(string("(Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].size() > 0 && Truth_Muon_Pt[Truth_Muon_Isol<Truth_Muon_Isol_Max].at(0) > ")+to_string(muon_leadpt)+string(") ||"
                  "(Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].size() > 0 && Truth_Electron_Pt[Truth_Electron_Isol<Truth_Electron_Isol_Max].at(0) > ")+to_string(electron_leadpt)
            +string(")"),string("lead lepton pt|requires highest pt lepton satisfies p_T > ")+to_string(muon_leadpt));       //(isolated) lepton lead p_T cut
      d_curr = &d_leadlep;

      rmin = 0; rmax = 600; bins = 80; 
      h = (TH1F*)d_curr->Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Truth_Sublead_Lepton_Pt")->Clone("Truth_Sublead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Sublead_Lepton_Pt"));
      auto d_sublep = d_curr->Filter(subptfuncstring,
         string("sublead lepton pt|requires second highest pt lepton satisfies p_T > ")+to_string(muon_subpt)+string(" if muon or p_T > ")+to_string(electron_subpt)+string(" if electron"));                                            //(isolated) lepton sub p_T cut
      d_curr = &d_sublep;

      //-------new subdivision---------------------------
      rmin = 0; rmax = 8; bins = 8; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{l}",""),bins,rmin,rmax},"Truth_Lepton_n")->Clone("Truth_Lepton_n");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Lepton_n"));
      auto d_4lep = d_curr->Filter("Truth_Lepton_n > 3", "four leptons present|at least four leptons created in event");    
      d_curr = &d_4lep;

      rmin = 0; rmax = 6; bins = 6; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{l}",""),bins,rmin,rmax},"Truth_Leptons_Contained_ATLAS")->Clone("Truth_Leptons_Contained_ATLAS");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Leptons_Contained_ATLAS"));
      auto d_eta = d_curr->Filter("Truth_Leptons_Contained_ATLAS > 3",
                                  string("four leptons contained|requires |eta| < ")+to_string(ATLAS_e_etamax)+string(" for electrons and |eta| < ")
                                  +to_string(ATLAS_mu_etamax)+string(" for muons, with at least 4 leptons satisfying this."));
      d_curr = &d_eta;

      //rmin = 0; rmax = 8; bins = 8; 
      //h = (TH1F*)d_curr->Histo1D({"",histname("n_{l}",""),bins,rmin,rmax},"Truth_Isolated_Lepton_n")->Clone("Truth_Isolated_Lepton_n");
      //h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Isolated_Lepton_n"));
      //auto d_lepisol = d_curr->Filter("Truth_Isolated_Lepton_n > 3", "four leptons isolated|at least four leptons in event satisfying isolation requirements i.e. sum of energy within deltaR=0.3 < 5%"); 
      //d_curr = &d_lepisol;
      //-------------------------------------------------

      rmin = 0; rmax = 4; bins = 4; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{ll}",""),bins,rmin,rmax},"Truth_Lepton_Pairs")->Clone("Truth_Lepton_Pairs");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Lepton_Pairs"));
      auto d_leppairs = d_curr->Filter("Truth_Lepton_Pairs > 1", "two lepton pairs|two Z candidates, i.e. opposite sign same flavour pairs");  //two pairs of opposite sign isolated leptons
      d_curr = &d_leppairs;

      rmin = 0; rmax = 140; bins = 140; 
      h = (TH1F*)d_curr->Histo1D({"",histname("m_{ll}","GeV"),bins,rmin,rmax},"Truth_Dilepton_M")->Clone("Truth_Dilepton_M");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Truth_Dilepton_M"));
      //--------new subdivision--------------------------
      auto d_oneshell = d_curr->Filter(string("fabs(Truth_Dilepton_M.at(0)-Z_MASS) < ")+to_string(z_shell_proximity),
            string("one on-shell Z|at least one of the lepton pairs satisfying |m_ll-m_Z| < ")+to_string(z_shell_proximity));
      d_curr = &d_oneshell;
      //-------------------------------------------------
      auto d_twoshell = d_curr->Filter(string("fabs(Truth_Dilepton_M.at(1)-Z_MASS) < ")+to_string(z_shell_proximity),
            string("two on-shell Zs|both lepton pairs satisfying |m_ll-m_Z| < ")+to_string(z_shell_proximity));
      d_curr = &d_twoshell;

      //-------new subdivision---------------------------
      rmin = 0; rmax = 6; bins = 6; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{j}",""),bins,rmin,rmax},"Jet_n")->Clone("Jet_n");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet_n"));
      auto d_twojet = d_curr->Filter("Jet_n > 1","two reconstructed jets|at least two jets detected in event");                                            //two jets
      d_curr = &d_twojet;
      //-------------------------------------------------

//      rmin = 0; rmax = 6; bins = 6; 
//      h = (TH1F*)d5a.Histo1D({"",histname("n_{j}",""),bins,rmin,rmax},"Jet_Good_n")->Clone("Jet_Good_n");
//      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet_Good_n"));
      auto d_isoljet = d_curr->Filter("Jet_Good_n > 1","jets isolated from electrons|two jets in event with no electrons within deltaR=0.3");               //two jets isolated from electrons
      d_curr = &d_isoljet;

      rmin = 0; rmax = 600; bins = 60; 
      h = (TH1F*)d_curr->Histo1D({"",histname("p_{T,j}","GeV"),bins,rmin,rmax},"Sublead_Jet_Pt")->Clone("Sublead_Jet_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Sublead_Jet_Pt"));
      auto d_subjet = d_curr->Filter(string("Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1) > ")+to_string(jet_subpt),
            string("sublead jet pt|jet with second highest pt satisfies p_T > ")+to_string(jet_subpt));                                                                       //sublead jet pt
      d_curr = &d_subjet;

      //rmin = -6; rmax = 6; bins = 40; 
      //h = (TH1F*)d6.Histo1D({"",histname("#eta_j","GeV"),bins,rmin,rmax},"Jet_Eta")->Clone("Jet_Eta");
      //h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet_Eta"));
      //auto d7 = d6.Filter(string("double max=")+to_string(jet_etamax)+string("; return \
      //      (fabs(Jet_Eta.at(0)) < max && fabs(Jet_Eta.at(1)) < max);"),string("jets within HCAL|jet eta < ")+to_string(jet_etamax));                         //both jets contained within HCAL

      rmin = 0; rmax = 2000; bins = 25; 
      h = (TH1F*)d_curr->Histo1D({"",histname("m_{jj}","GeV"),bins,rmin,rmax},"Dijet_M")->Clone("Dijet_M");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dijet_M"));
      auto d_EWKmjj = d_curr->Filter(string("Dijet_M > ")+to_string(dijet_mass_loose),
            string("EWK dijet mass|requires m_jj > ")+to_string(dijet_mass_loose));                                 //loose dijet mass constraint for electroweak selection
      d_curr = &d_EWKmjj;

      rmin = 0; rmax = 12; bins = 40; 
      h = (TH1F*)d_curr->Histo1D({"",histname("|#Delta_{jj}|",""),bins,rmin,rmax},"Jet12_Eta_Diff")->Clone("Jet12_Eta_Diff");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet12_Eta_Diff"));
      auto d_psrapdiff = d_curr->Filter(string("Jet12_Eta_Diff > ")+to_string(jet_eta_diff),
            string("jet pseudorapidity difference|requires |delta eta_jj| < ")+to_string(jet_eta_diff));    //jet psudorapidity difference
      d_curr = &d_psrapdiff;

      rmin = 0; rmax =2000; bins = 25; 
      h = (TH1F*)d_curr->Histo1D({"",histname("m_{jj}","GeV"),bins,rmin,rmax},"Dijet_M")->Clone("Dijet_M_2");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dijet_M_2"));
      auto d_VBSmjj = d_curr->Filter(string("Dijet_M > ")+to_string(dijet_mass_tight),
            string("VBS dijet mass|requires m_jj > ")+to_string(dijet_mass_tight));                                 //tight dijet mass constraint for VBS selection
      d_curr = &d_VBSmjj;

      //         .Filter("Jet_Eta.at(0)*Jet_Eta.at(1) < 0","jets opposite sign eta");                                        //jets have opposite sign eta

      for (auto &&cut : d.Report())
      {
         int split = cut.GetName().find("|"); 
         std::string cutName(cut.GetName().substr(0,split)), cutDesc(cut.GetName().substr(split+1));
         fprintf(of,"%-30spass=%-10llu all=%-10llu -- %-6.3f %%\tscaled=%-10f\n%s\n", (cutName+std::string(":")).c_str(), cut.GetPass(), cut.GetAll(), cut.GetEff(),cut.GetPass()*scale,
               (cutDesc==" ")?"":(std::string("(")+cutDesc+std::string(")\n")).c_str() );
      } 
      fprintf(of,"\n");
   }
}
