#include "../utils/constants.h"
using floats = ROOT::VecOps::RVec<float>;
void selection(string date)
{
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
         snprintf(rtn,500,"%s;%s %s%s%s; Events / %g %s",file.c_str(),xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
               (rmax-rmin+1)/bins,units.c_str());
         return rtn;
      };
      auto savefile = [&](string var){
         char rtn[500];
         snprintf(rtn,500,"/home/user108/y4p/graph_logs/%s/selection1_%s_%s.pdf",date.c_str(),var.c_str(),file.c_str());
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

      //Basic preselection:
      auto d0 = d_new.Filter(
                                                                                                                     //SINGLE LEPTONS:
            "(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 &&\
            Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > 27) ||"                                                         //single isolated muon, p_T > 27 GeV
            "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 &&\
            Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > 27) ||"                                             //single isolated (tight?) electron, p_T > 27 GeV
            "(Muon_n > 0 && Muon_Pt.at(0) > 52) ||"                                                                  //single muon, p_T > 52 GeV
            "(Electron_n > 0 && Electron_Pt.at(0) > 61) ||"                                                          //single electron, p_T > 61 GeV
                                                                                                                     //TWO LEPTONS:
            "(Muon_n > 1 && Muon_Pt.at(1) > 15) ||"                                                                  //two muons, each p_T > 15 GeV
            "(Muon_n > 1 && Muon_Pt.at(0) > 23 && Muon_Pt.at(1) > 9) ||"                                             //two muons, p_T > 23, 9 GeV
            "(Electron_n > 1 && Electron_Pt.at(1) > 18) ||"                                                          //two (v loose?) electrons, each p_T > 18 GeV
            "(Electron_n > 1 && Muon_n > 1 && ((Muon_Pt.at(0) > Electron_Pt.at(0) && Muon_Pt.at(0) > 25 && Electron_Pt.at(0) > 8) ||\
            (Electron_Pt.at(0) >= Muon_Pt.at(0) && Electron_Pt.at(0) > 25 && Muon_Pt.at(0) > 8)) )||"                //one electron one muon, p_T > 8, 25 GeV
            "(Electron_n > 1 && Muon_n > 1 && ((Muon_Pt.at(0) > Electron_Pt.at(0) && Muon_Pt.at(0) > 18 && Electron_Pt.at(0) > 15) ||\
            (Electron_Pt.at(0) >= Muon_Pt.at(0) && Electron_Pt.at(0) > 18 && Muon_Pt.at(0) > 15)) )||"               //one electron one muon, p_T > 18, 15 GeV
            "(Electron_n > 1 && Muon_n > 1 && ((Muon_Pt.at(0) > Electron_Pt.at(0) && Muon_Pt.at(0) > 27 && Electron_Pt.at(0) > 9) ||\
            (Electron_Pt.at(0) >= Muon_Pt.at(0) && Electron_Pt.at(0) > 27 && Muon_Pt.at(0) > 9)) )||"                //one electron one muon, p_T > 27, 9 GeV
                                                                                                                     //THREE LEPTONS:
            "(Electron_n > 2 && Electron_Pt.at(0) > 25 && Electron_Pt.at(2) > 13) ||"                                //three (loose?) electrons, p_T > 25, 13, 13 GeV
            "(Muon_n > 2 && Muon_Pt.at(2) > 7) ||"                                                                   //three muons, each p_T > 7 GeV
            "(Muon_n > 2 && Muon_Pt.at(0) > 21 && Muon_Pt.at(2) > 5) ||"                                             //three muons, p_T > 21, 5, 5 GeV
            "(Muon_n > 1 && Electron_n > 0 && ((Muon_Pt.at(0) > Electron_Pt.at(0) && Muon_Pt.at(0) > 13 && Muon_Pt.at(1) > 11 && Electron_Pt.at(0) > 11) ||\
            (Electron_Pt.at(0)>=Muon_Pt.at(0) && Electron_Pt.at(0) > 13 && Muon_Pt.at(1) > 11)) )||"                 //two muons one (loose?) electron, p_T > 11, 11, 13 GeV
            "(Muon_n > 0 && Electron_n > 1 && ((Muon_Pt.at(0) > Electron_Pt.at(1) && Muon_Pt.at(0) > 13 && Electron_Pt.at(1) > 11 && Electron_Pt.at(0) > 13) ||\
            (Electron_Pt.at(1)>=Muon_Pt.at(0) && Electron_Pt.at(1) > 13 && Muon_Pt.at(0) > 11)) )||"                 //two (loose?) electrons one muon, p_T > 13, 13, 11 GeV
                                                                                                                     //SINGLE JET
            "(Jet_n > 0 && Jet_Pt.at(0) > 435)"                                                                      //Jet (R=0.4), p_T > 435 GeV
            , "trigger simulation");
      
      rmin = 0; rmax = 599; bins = 80; 
      h = (TH1F*)d0.Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Lead_Lepton_Pt")->Clone("Lead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Lead_Lepton_Pt"));
      auto d1 = d0.Filter(string("(Muon_Pt[Muon_Isol<Muon_Isol_Max].size() > 0 && Muon_Pt[Muon_Isol<Muon_Isol_Max].at(0) > ")+to_string(muon_leadpt)+string(") ||"
            "(Electron_Pt[Electron_Isol<Electron_Isol_Max].size() > 0 && Electron_Pt[Electron_Isol<Electron_Isol_Max].at(0) > ")+to_string(electron_leadpt)
            +string(")"),"lead lepton pt");                                                                           //(isolated) lepton lead p_T cut

      rmin = 0; rmax = 599; bins = 80; 
      h = (TH1F*)d1.Histo1D({"",histname("p_{T,l}","GeV"),bins,rmin,rmax},"Sublead_Lepton_Pt")->Clone("Sublead_Lepton_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Sublead_Lepton_Pt"));
      auto d2 = d1.Filter([&muon_subpt,&electron_subpt](floats &mu_pt, floats &e_pt, floats &mu_isol, floats &e_isol, float mu_max, float e_max){
            int mu_size( mu_pt[mu_isol<mu_max].size() ), e_size( e_pt[e_isol<e_max].size() );
            double mu_lead( (mu_size > 0) ? mu_pt[mu_isol<mu_max].at(0) : 0 ), e_lead( (e_size > 0) ? e_pt[e_isol<e_max].at(0) : 0 );
            double lead, H0_sub, H0_subreq, H1_subreq, H1_isol_max; floats *H1, *H1_isol;
            //H0 is that the subleading lepton will be of opposite type to the leading.
            if ( e_lead > mu_lead ) {
               lead = e_lead; H0_sub = mu_lead; H0_subreq = muon_subpt; H1_subreq = electron_subpt;
               H1 = &e_pt; H1_isol = &e_isol; H1_isol_max = e_max;
            } else {
               lead = mu_lead; H0_sub = e_lead; H0_subreq = electron_subpt; H1_subreq = muon_subpt;
               H1 = &mu_pt; H1_isol = &mu_isol; H1_isol_max = mu_max; }
            for (auto &h1pt : (*H1)[*H1_isol<H1_isol_max]) if (h1pt > H0_sub && h1pt != lead) return (h1pt > H1_subreq);
            return (H0_sub > H0_subreq);
         },{"Muon_Pt","Electron_Pt","Muon_Isol","Electron_Isol","Muon_Isol_Max","Electron_Isol_Max"},
         "sublead lepton pt");                                                                                       //(isolated) lepton sub p_T cut

      rmin = 0; rmax = 4; bins = 5; 
      h = (TH1F*)d2.Histo1D({"",histname("n_{ll}",""),bins,rmin,rmax},"Lepton_Pairs")->Clone("Lepton_Pairs");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Lepton_Pairs"));
      auto d3 = d2.Filter("Lepton_Pairs > 1", "two lepton pairs");                                                    //two pairs of opposite sign isolated leptons

      rmin = 0; rmax = 140; bins = 141; 
      h = (TH1F*)d3.Histo1D({"",histname("m_{ll}","GeV"),bins,rmin,rmax},"Dilepton_M")->Clone("Dilepton_M");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dilepton_M"));
      auto d4 = d3.Filter(string("fabs(Dilepton_M.at(1)-Z_MASS) < ")+to_string(z_shell_proximity),
            "two lepton pairs on Z shell");                                                                          //both pairs on shell

      rmin = 0; rmax = 6; bins = 7; 
      h = (TH1F*)d4.Histo1D({"",histname("n_{j}",""),bins,rmin,rmax},"Jet_Good_n")->Clone("Jet_Good_n");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet_Good_n"));
      auto d5 = d4.Filter("Jet_Good_n > 1","two good jets");                                                         //two jets isolated from electrons with no b-tag

      rmin = 0; rmax = 599; bins = 80; 
      h = (TH1F*)d5.Histo1D({"",histname("p_{T,j}","GeV"),bins,rmin,rmax},"Sublead_Jet_Pt")->Clone("Sublead_Jet_Pt");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Sublead_Jet_Pt"));
      auto d6 = d5.Filter(string("Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1) > ")+to_string(jet_subpt),
            "sublead jet pt");                                                                                       //sublead jet pt

      rmin = -6; rmax = 6; bins = 130; 
      h = (TH1F*)d6.Histo1D({"",histname("#eta_j","GeV"),bins,rmin,rmax},"Jet_Eta")->Clone("Jet_Eta");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet_Eta"));
      auto d7 = d6.Filter(string("double max=")+to_string(jet_etamax)+string("; return \
            (fabs(Jet_Eta.at(0)) < max && fabs(Jet_Eta.at(1)) < max);"),"jets within HCAL");                         //both jets contained within HCAL

      rmin = 0; rmax = 1999; bins = 100; 
      h = (TH1F*)d7.Histo1D({"",histname("m_{jj}","GeV"),bins,rmin,rmax},"Dijet_M")->Clone("Dijet_M");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dijet_M"));
      auto d_EWK = d7.Filter(string("Dijet_M > ")+to_string(dijet_mass_loose),"EWK dijet mass");                                 //loose dijet mass constraint for electroweak selection

      rmin = 0; rmax = 12; bins = 13; 
      h = (TH1F*)d_EWK.Histo1D({"",histname("|#Delta_{jj}|",""),bins,rmin,rmax},"Jet12_Eta_Diff")->Clone("Jet12_Eta_Diff");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Jet12_Eta_Diff"));
      auto d8 = d_EWK.Filter(string("Jet12_Eta_Diff > ")+to_string(jet_eta_diff),"jet pseudorapidity difference");    //jet psudorapidity difference
   
      rmin = 0; rmax = 1999; bins = 100; 
      h = (TH1F*)d8.Histo1D({"",histname("m_{jj}","GeV"),bins,rmin,rmax},"Dijet_M")->Clone("Dijet_M_2");
      h->Scale(scale); h->Draw("hist"); gPad->SaveAs(savefile("Dijet_M_2"));
      auto d9 = d8.Filter(string("Dijet_M > ")+to_string(dijet_mass_tight),"VBS dijet mass");                                 //tight dijet mass constraint for VBS selection

//         .Filter("Jet_Eta.at(0)*Jet_Eta.at(1) < 0","jets opposite sign eta");                                        //jets have opposite sign eta

      

      //d.Report()->Print();
      for (auto &&cut : d.Report())
      {
         printf("%-30spass=%-10llu\t\tall=%-10llu\t--\t%-6.3f %%\t\tscaled=%-10f\n", (cut.GetName()+std::string(":")).c_str(), cut.GetPass(), cut.GetAll(), cut.GetEff(),cut.GetPass()*scale );
      } 
      printf("\n");
   }
}
