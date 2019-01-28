#include "../utils/constants.h"
using floats = ROOT::VecOps::RVec<float>;
void selection_fancyout(string date, bool CMS = false)
{
   string exper = CMS ? "CMS" : "ATLAS";
   FILE *of = fopen((string("/home/user108/y4p/cutflows/selection_cutflow_")+date+exper+string(".dat")).c_str(),"w");
   gStyle->SetOptStat(1111111);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(1,0);
   gStyle->SetLabelSize(.033,"xyz");
   gStyle->SetTitleSize(.038,"xyz");
   gStyle->SetLabelFont(42,"xyz");
   gStyle->SetTitleOffset(1.42,"y"); 
   gStyle->SetCanvasDefW(800); 
   gStyle->SetCanvasDefH(800); 
   TLatex ltx;
   ltx.SetTextSize(.04);
   ltx.SetTextFont(42);
   new TCanvas();
   gPad->SetLeftMargin(0.13);
   vector<TH1*> output_hists;
   std::string file = "";
   int cutnum = 0;
   auto savefile = [&](string var){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/selection1_%s%02d_%s_%s.pdf",date.c_str(),CMS?"CMS_":"",cutnum++,var.c_str(),file.c_str());
      return rtn;
   };

   for (std::string iFile : {string("ZZjj_")+exper+string("_500K"),string("inclusive_ATLAS_500K")})
   {
      file = iFile;
      //setup:
      double target_luminosity = 35.9; //As CMS analysis
      std::string pre("/home/user108/y4p/root_output/"),post("_o.root");

      ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
      int n_events = d.Histo1D("Cross_Section")->GetEntries();//*(d.Take<int>("Event_Count").begin());
      double cross_section = *(d.Take<double>("Cross_Section").begin());
      double unweighted_scale = target_luminosity*cross_section/n_events;
      double weightsum = d.Sum("Event_Weight").GetValue();
      double scale = n_events*unweighted_scale/weightsum;

      printf("%s sum of weights: %e\n",file.c_str(),weightsum);

      int bins(1); double rmin(0), rmax(1); char histstr[100];
      auto histname = [&bins,&rmin,&rmax,&file](string xlabel, string units){
         char rtn[500]; bool nu=(units=="");
         snprintf(rtn,500,";%s %s%s%s; Events / %g %s"/*,file.c_str()*/,xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
               (rmax-rmin)/bins,units.c_str());
         return rtn;
      };
      TH1F *h;
      std::vector<double> cut_event_weights;

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

      auto d_plus = d_new.Define("Relative_Event_Weight",string("Event_Weight/")+to_string(weightsum));

      auto d_start = d_plus.Filter("true");
      auto *d_curr = &d_start;
      cutnum = 0;

      auto make_cut = [&d,&d_curr,&scale,&unweighted_scale,&histname,&savefile,&bins,&rmin,&rmax,&cut_event_weights,&ltx,&file]
      (string filterFunc, string filterName, bool makeHist=false,string plotColumn="", string graphName="", int l_bins=0, double l_rmin=0, double l_rmax=0, string xlabel="", string unit="",
        string latexdesc = "", double latexX = 0.64, double latexY = 0.85)
      {
         if (makeHist)
         {
            bins = l_bins; rmin=l_rmin; rmax=l_rmax;
            TH1F* h = (TH1F*)d_curr->Histo1D({"",histname(xlabel,unit),bins,rmin,rmax},plotColumn/*,"Event_Weight"*/)->Clone(graphName.c_str());
            h->Scale(unweighted_scale); h->Draw("hist"); 
            h->SetLineWidth(2); h->SetLineColor(1);
            gPad->GetCanvas()->Update();
            double ltxX(latexX*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin());
            double ltxY(latexY*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin());
            string sampleType = (file.substr(0,3) == "inc") ? "Incl." : "EW";
            ltx.DrawLatex(ltxX,ltxY,(string("#splitline{#splitline{#bf{#it{Selection}} }{")+latexdesc+string("}}{#splitline{")+
            sampleType+string(" ZZjj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}}")).c_str());
            gPad->SaveAs(savefile(graphName));
         }
         auto d_rtn = d_curr->Filter(filterFunc,filterName);
         cut_event_weights.push_back(d_rtn.Sum("Event_Weight").GetValue());
         d_curr = &d_rtn;
         return d_rtn;
      };
      auto d_none = make_cut("true","initial data");

      //Basic preselection:
      auto dtrig = make_cut(
            //SINGLE LEPTONS:
            "(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 &&\
            Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > 27) ||"                                                         //single isolated muon, p_T > 27 GeV
            "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 &&\
            Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > 27) ||"                                             //single isolated (tight?) electron, p_T > 27 GeV
            "(Muon_n > 0 && Muon_Pt.at(0) > 52) ||"                                                                  //single muon, p_T > 52 GeV
            "(Electron_n > 0 && Electron_Pt.at(0) > 61)"                                                             //single electron, p_T > 61 GeV
            , "trigger simulation|ATLAS single lepton trigger requirements");

      //int i = d_curr; //ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> *

      auto d_leadlep = make_cut(
         string("(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 && Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > ")+to_string(muon_leadpt)+string(") ||"
            "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 && Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > ")+to_string(electron_leadpt)
            +string(")"),
         string("lead lepton pt|requires highest pt lepton satisfies p_T > ")+to_string(muon_leadpt),
         true,"Lead_Lepton_Pt","Lead_Lepton_Pt",80, 0, 600,"p_{T,l}","GeV","Lead lepton p_{T}"
      );    //(isolated) lepton lead p_T cut

      auto d_sublep = make_cut(subptfuncstring,
         string("sublead lepton pt|requires second highest pt lepton satisfies p_T > ")+to_string(muon_subpt)+string(" if muon or p_T > ")+to_string(electron_subpt)+string(" if electron"),
         true,"Sublead_Lepton_Pt","Sublead_Lepton_Pt",80,0,600,"p_{T,l}","GeV","Sublead lepon p_{T}"); //(isolated) lepton sub p_T cut

      auto d_4lep = make_cut("Lepton_n > 3", "four reconstructed leptons|at least four leptons detected in event",
         true,"Lepton_n","Lepton_n",8,0,8,"n_{l}","","Num. leptons");    

      auto d_lepisol = make_cut("Isolated_Lepton_n > 3", "four leptons isolated|at least four leptons in event satisfying isolation requirements i.e. sum of energy within deltaR=0.3 < 5%",
         true,"Isolated_Lepton_n","Isolated_Lepton_n",8,0,8,"n_{l}","","Lepton isolation",0.05); 

      auto d_leppairs = make_cut("Lepton_Pairs > 1", "two lepton pairs|two Z candidates, i.e. opposite sign same flavour pairs",
         true,"Lepton_Pairs","Lepton_Pairs",4,0,4,"n_{ll}","","Lepton pairs",0.05);  //two pairs of opposite sign isolated leptons

      auto d_oneshell = make_cut(string("fabs(Dilepton_M.at(0)-Z_MASS) < ")+to_string(z_shell_proximity),
            string("one on-shell Z|at least one of the lepton pairs satisfying |m_ll-m_Z| < ")+to_string(z_shell_proximity),
            true,"Dilepton_M","Dilepton_M",140,0,140,"m_{ll}","GeV","Dilepton mass",0.05);

      auto d_twoshell = make_cut(string("fabs(Dilepton_M.at(1)-Z_MASS) < ")+to_string(z_shell_proximity),
            string("two on-shell Zs|both lepton pairs satisfying |m_ll-m_Z| < ")+to_string(z_shell_proximity));

      auto d_twojet = make_cut("Jet_n > 1","two reconstructed jets|at least two jets detected in event",
         true,"Jet_n","Jet_n",6,0,6,"n_{j}","","Num. jets");                                            //two jets

      auto d_isoljet = make_cut("Jet_Good_n > 1","jets isolated from electrons|two jets in event with no electrons within deltaR=0.3");               //two jets isolated from electrons

      auto d_subjet = make_cut(string("Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1) > ")+to_string(jet_subpt),
            string("sublead jet pt|jet with second highest pt satisfies p_T > ")+to_string(jet_subpt),
            true,"Sublead_Jet_Pt","Sublead_Jet_Pt",60,0,600,"p_{T,j}","GeV","Sublead jet p_{T}");                                                                       //sublead jet pt

      auto d_EWKmjj = make_cut(string("Dijet_M > ")+to_string(dijet_mass_loose),
            string("EWK dijet mass|requires m_jj > ")+to_string(dijet_mass_loose),
            true,"Dijet_M","Dijet_M",25,0,2000,"m_{jj}","GeV","Dijet mass");                                 //loose dijet mass constraint for electroweak selection

      auto d_psrapdiff = make_cut(string("Jet12_Eta_Diff > ")+to_string(jet_eta_diff),
            string("jet pseudorapidity difference|requires |delta eta_jj| > ")+to_string(jet_eta_diff),
            true,"Jet12_Eta_Diff","Jet12_Eta_Diff",40,0,12,"|#Delta#eta_{jj}|","","Jet #Delta#eta");    //jet psudorapidity difference

      auto d_VBSmjj = make_cut(string("Dijet_M > ")+to_string(dijet_mass_tight),
            string("VBS dijet mass|requires m_jj > ")+to_string(dijet_mass_tight),
            true,"Dijet_M","Dijet_M_2",25,0,2000,"m_{jj}","GeV","Dijet mass");                                 //tight dijet mass constraint for VBS selection

      //         .Filter("Jet_Eta.at(0)*Jet_Eta.at(1) < 0","jets opposite sign eta");                                        //jets have opposite sign eta
      int iCut = 0;
      for (auto &&cut : d.Report())
      {
         int split = cut.GetName().find("|"); 
         std::string cutName(cut.GetName().substr(0,split)), cutDesc(cut.GetName().substr(split+1));
         double weight = cut_event_weights.at(iCut);
         fprintf(of,"%-30spass=%-10llu all=%-10llu -- %-6.3f %%\tnaive=%-10f weights=%.3e\tscaled=%-10f\n%s\n", (cutName+std::string(":")).c_str(), cut.GetPass(), cut.GetAll(), cut.GetEff(),
               cut.GetPass()*unweighted_scale,weight,/*cut.GetPass()**/weight*scale,
               (cutDesc==" ")?"":(std::string("(")+cutDesc+std::string(")\n")).c_str() );
         iCut++;
      } 
      fprintf(of,"\n");

      h = (TH1F*)d_EWKmjj.Histo1D({"","",15,100,1600},"Dijet_M")->Clone("");
      h->Scale(unweighted_scale);
      output_hists.push_back(h);
      h = (TH1F*)d_EWKmjj.Histo1D({"","",14,0,7},"Jet12_Eta_Diff")->Clone("");
      h->Scale(unweighted_scale);
      output_hists.push_back(h);

      cout << d_VBSmjj.Sum("Event_Weight").GetValue() << endl;
      //d_curr->Snapshot("EventTree",pre+file+std::string("_filtered.root"));
   }
   file = "ATLAS_Combined_500K";
   auto stack_hists = [](string name, string label, TH1* hsig, TH1* hbg){
      auto hs = new THStack(name.c_str(),label.c_str());

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
   h_ewk->Draw("hist"); gPad->SetLogy(1); h_ewk->SetMaximum(120);
            gPad->GetCanvas()->Update();
            double ltxX(.64*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin());
            double ltxY(.9*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin());
            ltx.DrawLatex(ltxX,95,"#splitline{#bf{#it{EW Selection}}}{#splitline{X#rightarrow jj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}}");
   gPad->SaveAs(savefile("EWK_Dijet_M"));

   new TCanvas();
   gPad->SetLeftMargin(0.13);
   auto h_vbs = stack_hists("",";|#Delta#eta_{jj}|;Events / 0.5",output_hists[1],output_hists[3]);
   h_vbs->Draw("hist"); h_vbs->SetMaximum(22);
            gPad->GetCanvas()->Update();
            ltxX = 0.64*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin();
            ltxY = 0.9*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin();
            ltx.DrawLatex(ltxX,ltxY,"#splitline{#bf{#it{EW Selection}}}{#splitline{X#rightarrow jj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}}");
   gPad->SaveAs(savefile("EWK_DelEtJJ"));
}
