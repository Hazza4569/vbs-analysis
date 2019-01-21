#include "../utils/constants.h"
using floats = ROOT::VecOps::RVec<float>;
void selection(string date, bool CMS = false)
{
   string exper = CMS ? "CMS" : "ATLAS";
   FILE *of = fopen((string("/home/user108/y4p/cutflows/selection_cutflow_")+date+exper+string(".dat")).c_str(),"w");
   new TCanvas();
   gStyle->SetOptStat(1111111);
   vector<TH1*> output_hists;
   for (std::string file : {string("ZZjj_")+exper+string("_500K"),string("inclusive_ATLAS_500K")})
   {
      //setup:
      double target_luminosity = 35.9; //As CMS analysis
      std::string pre("/home/user108/y4p/root_output/"),post("_o.root");

      ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
      int n_events = d.Histo1D("Cross_Section")->GetEntries();//*(d.Take<int>("Event_Count").begin());
      double cross_section = *(d.Take<double>("Cross_Section").begin());
      Double_t scale = target_luminosity*cross_section/n_events;

      int bins(1); double rmin(0), rmax(1); char histstr[100];
      auto histname = [&bins,&rmin,&rmax,&file](string xlabel, string units){
         char rtn[500]; bool nu=(units=="");
         snprintf(rtn,500,"%s;%s %s%s%s; Events / %g %s",file.c_str(),xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
               (rmax-rmin)/bins,units.c_str());
         return rtn;
      };
      int cutnum = 0;
      auto savefile = [&](string var){
         char rtn[500];
         snprintf(rtn,500,"/home/user108/y4p/graph_logs/%s/selection1_%s%02d_%s_%s.pdf",date.c_str(),CMS?"CMS_":"",cutnum++,var.c_str(),file.c_str());
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

      auto d_new = d.Define("Lead_Lepton_Pt","Isolated_Lepton_Pt.size() < 1 ? -1 : Isolated_Lepton_Pt.at(0)")
         .Define("Sublead_Lepton_Pt","Isolated_Lepton_Pt.size() < 2 ? -1 : Isolated_Lepton_Pt.at(1)")
         .Define("Sublead_Jet_Pt","Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].size() < 2 ? -1 : Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1)");

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
            "int mu_size( Muon_Pt[Muon_Isol<Muon_Isol_Max].size() ), e_size( Electron_Pt[Electron_Isol<Electron_Isol_Max].size() );"
            "double mu_lead( (mu_size > 0) ? Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) : 0 ), e_lead( (e_size > 0) ? Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) : 0 );"
            "double lead, H0_sub, H0_subreq, H1_subreq, H1_isol_max; ROOT::VecOps::RVec<float> *H1, *H1_isol;"
            //H0 is that the subleading lepton will be of opposite type to the leading.
            "if ( e_lead > mu_lead ) {"
            "   lead = e_lead; H0_sub = mu_lead; H0_subreq = muon_subpt; H1_subreq = electron_subpt;"
            "   H1 = &Electron_Pt; H1_isol = &Electron_Isol; H1_isol_max = Electron_Isol_Max;"
            "} else {"
            "   lead = mu_lead; H0_sub = e_lead; H0_subreq = electron_subpt; H1_subreq = muon_subpt;"
            "   H1 = &Muon_Pt; H1_isol = &Muon_Isol; H1_isol_max = Muon_Isol_Max;"
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
            "(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 &&\
            Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > 27) ||"                                                         //single isolated muon, p_T > 27 GeV
            "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 &&\
            Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > 27) ||"                                             //single isolated (tight?) electron, p_T > 27 GeV
            "(Muon_n > 0 && Muon_Pt.at(0) > 52) ||"                                                                  //single muon, p_T > 52 GeV
            "(Electron_n > 0 && Electron_Pt.at(0) > 61)"                                                             //single electron, p_T > 61 GeV
            , "trigger simulation|ATLAS single lepton trigger requirements");
      d_curr = &dtrig;


      rmin = 0; rmax = 600; bins = 80; 
      h = (TH1F*)d_curr->Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Lead_Lepton_Pt")->Clone("Lead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Lead_Lepton_Pt"));
      auto d_leadlep = d_curr->Filter(string("(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 && Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > ")+to_string(muon_leadpt)+string(") ||"
                  "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 && Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > ")+to_string(electron_leadpt)
            +string(")"),string("lead lepton pt|requires highest pt lepton satisfies p_T > ")+to_string(muon_leadpt));       //(isolated) lepton lead p_T cut
      d_curr = &d_leadlep;

      rmin = 0; rmax = 600; bins = 80; 
      h = (TH1F*)d_curr->Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Sublead_Lepton_Pt")->Clone("Sublead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Sublead_Lepton_Pt"));
      auto d_sublep = d_curr->Filter(subptfuncstring,
         string("sublead lepton pt|requires second highest pt lepton satisfies p_T > ")+to_string(muon_subpt)+string(" if muon or p_T > ")+to_string(electron_subpt)+string(" if electron"));                                            //(isolated) lepton sub p_T cut
      d_curr = &d_sublep;

      //-------new subdivision---------------------------
      rmin = 0; rmax = 8; bins = 8; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{l}",""),bins,rmin,rmax},"Lepton_n")->Clone("Lepton_n");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Lepton_n"));
      auto d_4lep = d_curr->Filter("Lepton_n > 3", "four reconstructed leptons|at least four leptons detected in event");    
      d_curr = &d_4lep;

      rmin = 0; rmax = 8; bins = 8; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{l}",""),bins,rmin,rmax},"Isolated_Lepton_n")->Clone("Isolated_Lepton_n");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Isolated_Lepton_n"));
      auto d_lepisol = d_curr->Filter("Isolated_Lepton_n > 3", "four leptons isolated|at least four leptons in event satisfying isolation requirements i.e. sum of energy within deltaR=0.3 < 5%"); 
      d_curr = &d_lepisol;
      //-------------------------------------------------

      rmin = 0; rmax = 4; bins = 4; 
      h = (TH1F*)d_curr->Histo1D({"",histname("n_{ll}",""),bins,rmin,rmax},"Lepton_Pairs")->Clone("Lepton_Pairs");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Lepton_Pairs"));
      auto d_leppairs = d_curr->Filter("Lepton_Pairs > 1", "two lepton pairs|two Z candidates, i.e. opposite sign same flavour pairs");  //two pairs of opposite sign isolated leptons
      d_curr = &d_leppairs;

      rmin = 0; rmax = 140; bins = 140; 
      h = (TH1F*)d_curr->Histo1D({"",histname("m_{ll}","GeV"),bins,rmin,rmax},"Dilepton_M")->Clone("Dilepton_M");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dilepton_M"));
      //--------new subdivision--------------------------
      auto d_oneshell = d_curr->Filter(string("fabs(Dilepton_M.at(0)-Z_MASS) < ")+to_string(z_shell_proximity),
            string("one on-shell Z|at least one of the lepton pairs satisfying |m_ll-m_Z| < ")+to_string(z_shell_proximity));
      d_curr = &d_oneshell;
      //-------------------------------------------------
      auto d_twoshell = d_curr->Filter(string("fabs(Dilepton_M.at(1)-Z_MASS) < ")+to_string(z_shell_proximity),
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

      h = (TH1F*)d_EWKmjj.Histo1D({"","",15,100,1600},"Dijet_M")->Clone("");
      h->Scale(scale);
      output_hists.push_back(h);
      h = (TH1F*)d_VBSmjj.Histo1D({"","",16,400,2000},"Dijet_M")->Clone("");
      h->Scale(scale);
      output_hists.push_back(h);
   }
   auto stack_hists = [](string name, string label, TH1* hsig, TH1* hbg){
      auto hs = new THStack();

      hbg->Add(hsig,-1);
      for (int iBin = 1; iBin <= hbg->GetXaxis()->GetNbins(); iBin++) if (hbg->GetBinContent(iBin) < 0) hbg->SetBinContent(iBin,0);

      hsig->SetFillColor(6);
      hbg->SetFillColor(kAzure-6);
      hbg->SetLineColorAlpha(kBlack,1);
      hsig->SetLineColorAlpha(kBlack,1);

      hs->Add(hbg);
      hs->Add(hsig);
      return hs;
   };
   auto h_ewk = stack_hists("ZZjj Dijet Mass",";m_{jj} [GeV];Events / 100 GeV",output_hists[0],output_hists[2]);
   h_ewk->Draw("hist");

   new TCanvas();
   auto h_vbs = stack_hists("VBS Dijet Mass",";m_{jj} [GeV];Events / 100 GeV",output_hists[1],output_hists[3]);
   h_vbs->Draw("hist");
}
