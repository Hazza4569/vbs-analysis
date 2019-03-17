using floats = ROOT::VecOps::RVec<float>;
#include "../utils/constants.h"
void angular_2(string date, string calcmethod="default", string extralabel="", bool verbose=false)
{
   //setup files:
   double target_luminosity = 3000; //Hi-Lumi
   double mjj = 1.037639e+03;
   double njj = 2.556276e+00; //as per optimisation

   std::string pre("/home/user108/y4p/root_output/"),post_uncut("_o.root"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_1M";
   string strBkg = "inclusive_ATLAS_5M";

   ROOT::RDataFrame d_uncut("EventTree",(pre+strSig+post_uncut).c_str());
   ROOT::RDataFrame d_bguncut("EventTree",(pre+strBkg+post_uncut).c_str());
   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());

   string abnwcut = "fabs(Event_Weight) < 40";

   double scale_sig = target_luminosity *
      (*d_uncut.Take<double>("Cross_Section").begin()) / 
      d_uncut.Filter(abnwcut).Sum("Event_Weight").GetValue() ;

   double scale_bkg = target_luminosity *
      (*d_bguncut.Take<double>("Cross_Section").begin()) / 
      d_bguncut.Filter(abnwcut).Sum("Event_Weight").GetValue() ;

   auto h = (TH1F*) d_uncut.Histo1D("Cross_Section")->Clone();
   printf("%.3f",h->GetEntries()/h->GetMean());

   std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
   auto d_sigcut = d_sig.Filter(cut);
   auto d_bkgcut = d_bkg.Filter(cut);

   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/angular_%s%s_%s.%s",date.c_str(),calcmethod.c_str(),extralabel.c_str(),var.c_str(),ft.c_str());
      return string(rtn);
   };

   auto combine_pairs = [](auto a, auto b) { return make_pair( a.first+a.second, b.first+b.second ); };
   auto pair_fourvec = [](auto a) { return (a.first+a.second); };
   //auto pair_Zorder = [](auto l, auto r) { double lm((l.first+l.second).M()), rm((r.first+r.second).M()); return fabs(lm-Z_MASS) < fabs(rm-Z_MASS);  };

   auto cos_theta_star = [&verbose,&combine_pairs,&pair_fourvec,&calcmethod]
                         (floats &mu_pt, floats &mu_eta, floats &mu_phi, float mu_m, floats &mu_q,
                          floats &e_pt , floats &e_eta , floats &e_phi , float e_m , floats &e_q )
                        ->ROOT::VecOps::RVec<float>
   {
      floats rtn;
      // *********************************** LEPTON PAIR RECO *********************************** //
      // Find two Z bosons from the leptons. The pairings chosen will be those which have the closest
      // reconstructed masses to the nominal Z mass.
      
      //make loopable over lepton types: vectors of RVecs
      vector<ROOT::VecOps::RVec<float>> l_pts{mu_pt,e_pt}, l_etas{mu_eta,e_eta},l_phis{mu_phi,e_phi}, l_qs{mu_q,e_q};
      vector<float> l_ms{mu_m,e_m};
      vector<int> l_ns{(int)mu_pt.size(),(int)e_pt.size()}; //not checking for inconsistent lengths
      
      //events with < 4 leptons aren't usable
      if ( l_ns[0]+l_ns[1] < 4 ) return rtn;

      //set up ordered set
      struct pair_Zorder { bool operator() (const pair<TLorentzVector,TLorentzVector>& l, const pair<TLorentzVector,TLorentzVector>& r) const
         { double lm((l.first+l.second).M()), rm((r.first+r.second).M()); return fabs(lm-Z_MASS) < fabs(rm-Z_MASS);  }}; 

      set< pair<TLorentzVector,TLorentzVector>, pair_Zorder >  l_pairs;

      //loop over lepton type:
      for (int ilep = 0; ilep < 2; ilep++)
      {
         //form lists of each charge:       (not checking for abnormal charges here, this should be checked [poss. manually] in advance.)
         vector<TLorentzVector> l_negatives, l_positives;
         for ( int i = 0; i < l_ns[ilep]; i++ )
         {
            TLorentzVector curr_lepton;
            curr_lepton.SetPtEtaPhiM(l_pts[ilep].at(i), l_etas[ilep].at(i), l_phis[ilep].at(i), l_ms[ilep]);
            ( (l_qs[ilep].at(i)<0) ? l_negatives : l_positives ).push_back(curr_lepton);
         }
         //add pairs into ordered set:
         for (auto l_pos : l_positives) for (auto l_neg : l_negatives) l_pairs.insert(pair<TLorentzVector,TLorentzVector>(l_neg,l_pos));
      }

      if ( l_pairs.size() < 2 ) return rtn;
      //NOTE: this does not completely check that physically possible pairs have been made (i.e. no double use of leptons).
      //this should be done with filters beforehand. Checks here are just for safe memory access not good physics.

      //extract best pairs: 
      //Z1 -> l1(-), l2(+)   Z2 -> l3(-), l4(+)
      vector<pair<TLorentzVector,TLorentzVector>> pairs; 
      for (auto it=l_pairs.begin(); pairs.size() < 2; it++) pairs.push_back(*it);
      // **************************************************************************************** //

      TVector3 diboson_boost = -pair_fourvec( combine_pairs(pairs.at(0),pairs.at(1)) ).BoostVector();

      //loop over the two lepton pairs:
      for (int pair_iter = 0; pair_iter < 2; pair_iter++)
      {
         auto cur_pair = pairs.at(pair_iter);
         auto alt_pair = pairs.at(1-pair_iter);

         //find theta*:
         auto pair_boost = -pair_fourvec(cur_pair).BoostVector();
         TLorentzVector lneg(cur_pair.first), lpos(cur_pair.second),
                        Z_cur(pair_fourvec(cur_pair)), Z_alt(pair_fourvec(alt_pair));

         double costhetastar;

         //first calculation
         //Z_cur.Boost(diboson_boost);
         //costhetastar = cos( Z_cur.Angle( lneg.Vect() ) );

         if (calcmethod == "default")
         {
            lneg.Boost(pair_boost);
            lpos.Boost(pair_boost);

            Z_alt.Boost(pair_boost);
            costhetastar = cos( (-Z_alt.Vect()).Angle( lneg.Vect() ) );
         }
         else if (calcmethod == "labZ")
         {
            lneg.Boost(pair_boost);
            lpos.Boost(pair_boost);

            costhetastar = cos( (-Z_cur.Vect()).Angle( lneg.Vect() ) );
         }

         rtn.push_back(costhetastar);
      }
      return rtn;
   };

   auto definitions = [&cos_theta_star](auto df)
   {
      return df
         .Define("anglediff",cos_theta_star,{"Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Mass","Truth_Muon_Charge",
                        "Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Mass","Truth_Electron_Charge"})
         .Define("anglediff1","anglediff.at(0)")
         .Define("anglediff2","anglediff.at(1)").Define("Z_Weight","ROOT::VecOps::RVec<float>(2,Event_Weight)");
   };

   auto d_uncutd = definitions(d_uncut).Filter("true");
   auto d_bguncutd = definitions(d_bguncut).Filter("true");
   auto d_sigcutd = definitions(d_sigcut).Filter("true");
   auto d_bkgcutd = definitions(d_bkgcut).Filter("true");

   //int nbins(11);
   //double xbins[nbins];
   //for (int i=0; i<=nbins; i++) xbins[i] = pow(10,-nbins+i);
   
   auto plot_angular = [&savefile](auto df, double scale, string title="", int order=0, double rangelim=1, string which="separate",string drawcom = "")
   {
      auto d_new = df.Filter("fabs(Event_Weight) < 40","Remove anomalous weightings");


      int nbins(16);

      TF1* f = new TF1("f","[0]*3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + (1.-[1]-[2])*(1.+2.*[3]*x+x**2) )/8.",-rangelim,rangelim);
      f->SetParameters(1,0.33,0.33,1);
      f->SetParNames("Scale","f_{L}","f_{+}","A");
      f->SetParLimits(1,0,1);
      f->SetParLimits(2,0,1);
      f->FixParameter(3,0.16/*8./(1-4* 0.23)*/);

      gStyle->SetOptStat(/*111111*/0); gStyle->SetOptFit(0);
      gStyle->SetStatX(0.7); gStyle->SetStatY(0.9);
      auto histtitle = [](int bins, string htitle, string afterthetastar="")
      {
         char rtn[500];
         snprintf(rtn,500,"%s;cos#theta*%s;Events / %g",htitle.c_str(),afterthetastar.c_str(),2./bins);
         return string(rtn);
      };

      if (which == "separate")
      {
         TH1F* h1 = (TH1F*) d_new.Histo1D({"",histtitle(nbins,title,"_{1}").c_str(),nbins,-1,1},"anglediff1","Event_Weight")->Clone(""); 
         h1->Scale(scale);
         TH1F* h2 = (TH1F*) d_new.Histo1D({"",histtitle(nbins,title,"_{2}").c_str(),nbins,-1,1},"anglediff2","Event_Weight")->Clone(""); 
         h2->Scale(scale);

         //      TF1* f = new TF1("f","[0]*3*( [1]*2*(1-x**2) + [2]*(1+x)**2 + (1-[1]-[2])*(1-x)**2 )/8.",-1,1); //      f->SetParameters(1,0.33,0.33);
         //      f->SetParNames("Scale","f_{L}","f_{+}");
         //      f->SetParLimits(1,0,1);
         //      f->SetParLimits(2,0,1);

         h1->SetMinimum(0);
         h1->Draw(drawcom.c_str());
         f->FixParameter(0,2.*h1->Integral()/nbins); //integral over bin width.
         h1->Fit("f");
         f->Draw("same");
         gPad->SaveAs(savefile(to_string(order)+title+string("_Z1"),"pdf").c_str());

         h2->SetMinimum(0);
         h2->Draw(drawcom.c_str());
         f->FixParameter(0,2.*h2->Integral()/nbins); //integral over bin width.
         h2->Fit("f");
         f->Draw("same");
         gPad->SaveAs(savefile(to_string(order)+title+string("_Z2"),"pdf").c_str());
      }
      else if (which == "combined")
      {
         TH1F* h1 = (TH1F*) d_new.Histo1D({"",histtitle(nbins,title).c_str(),nbins,-1,1},"anglediff","Z_Weight")->Clone(""); 
         h1->Scale(scale);

         h1->SetMinimum(0);
         h1->Draw(drawcom.c_str());
         f->FixParameter(0,2.*h1->Integral()/nbins); //integral over bin width.
         h1->Fit("f");
         f->Draw("same");

         //Extra lines:
         TF1 *f2 = new TF1("f2","[0]*3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + [4]*(1.+2.*[3]*x+x**2) )/8.",-rangelim,rangelim);
         for (int i=0; i < 4; i++) f2->SetParameter(i,f->GetParameter(i));
         f2->SetParameter(4,1-f->GetParameter(1)-f->GetParameter(2));
         TF1 *f3 = new TF1("f3","f2",-rangelim,rangelim);
         TF1 *f4 = new TF1("f4","f2",-rangelim,rangelim);
         f2->SetParameter(1,0); f2->SetParameter(4,0); f2->Draw("same"); f2->SetLineColor(kSpring);
         f3->SetParameter(2,0); f3->SetParameter(1,0); f3->Draw("same"); f3->SetLineColor(kCyan);
         f4->SetParameter(2,0); f4->SetParameter(4,0); f4->Draw("same"); f4->SetLineColor(kMagenta);
         for (int i=0; i < 4; i++) printf("%f ",f4->GetParameter(i));
         printf("\n");

         gPad->SaveAs(savefile(to_string(order)+title+string("_both"),"pdf").c_str());
         FILE* of = fopen(savefile(to_string(order)+title+string("_both"),"dat").c_str(),"w");
         fprintf(of,"Scale: %f\nf_L  :%f +/- %f\nf_+  :%f +/- %f\nf_-  :%f\nA    :%f\n\nStats:\nEntries: %f\nMean: %e\n Std. Dev: %e\nOverflow: %e\nUnderflow: %e\nIntegral: %e",
                  f->GetParameter(0),
                  f->GetParameter(1), f->GetParError(1),
                  f->GetParameter(2), f->GetParError(2),
                  1 - f->GetParameter(1) - f->GetParameter(2),
                  f->GetParameter(3),
                  h1->GetEntries(), h1->GetMean(), h1->GetStdDev(), h1->GetBinContent(h1->GetNbinsX()+1), h1->GetBinContent(0), h1->Integral());
         fclose(of);
      }
   }; 

   new TCanvas();
   //plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1"), scale_sig, "truth_preselection" ,0);
   //plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;") .Filter("fabs(Truth_Dilepton_M.at(0) - 91.19) < 3") , scale_sig, "truth_Z1tight" ,1);
   //plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;") .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 3","Z mass (3 GeV)") , scale_sig, "truth_Z2tight" ,1);
   //plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;") .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 8") , scale_sig, "truth_Z2med" ,1);
//   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
//         .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 5","mass within 5") .Filter("Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > 600 && Jet12_Eta_Diff > 3.6 ","jet cuts")
//         , scale_sig, "truth_mainasuggested" ,0,1,"combined","hist");
//   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
//         .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 5","mass within 5")
//         .Filter("Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > 600 && Jet12_Eta_Diff > 3.6 ","jet cuts")
//         .Filter("Tetralepton_M > 300 && Truth_Lepton_Pt.at(3) > 20","lepton cuts")
//         , scale_sig, "truth_mainasuggested+" ,1,.75,"combined","hist");
   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
         .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 5","mass within 5") .Filter("Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > 500 && Jet12_Eta_Diff > 2","jet cuts")
         , scale_sig, "truth_jetcuts" ,0,1,"combined","hist");
   plot_angular(d_bguncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
         .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 5","mass within 5") .Filter("Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > 500 && Jet12_Eta_Diff > 2","jet cuts")
         , scale_bkg, "bg_truth_jetcuts" ,0,1,"combined","hist");

   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
         .Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 5","mass within 5") .Filter("Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > 500 && Jet12_Eta_Diff > 2","jet cuts")
         .Filter("Truth_Lepton_Pt.at(3) > 20 && Truth_Tetralepton_M > 300").Filter("for (auto eta : Truth_Lepton_Eta) if (eta > 2.5) return false; return true;")
         , scale_sig, "truth_jetandleptcuts" ,0,1,"combined","hist");

//   //Look at effect of tight mass cuts on other distributions:
//   gStyle->SetStatX(1); gStyle->SetStatY(1);
//   auto d_this = d_uncutd.Filter("Truth_Lepton_Pairs > 1").Define("Truth_Dilepton_Eta1","Truth_Dilepton_Eta.at(0)").Define("Truth_Dilepton_Pt1","Truth_Dilepton_Pt.at(0)")
//                  .Define("Truth_Dilepton_Phi1","Truth_Dilepton_Phi.at(0)").Define("Truth_Dilepton_M1","Truth_Dilepton_M.at(0)").Filter("fabs(Event_Weight) < 40","Remove anomalous weightings");
//   auto plot_pt  = [&](auto d) { TH1F* htmp = (TH1F*)d.Histo1D({"",";#p_{T,Z_{1}};",300,0,600},"Truth_Dilepton_Pt1", "Event_Weight")->Clone(""); htmp->Scale(scale_sig); htmp->Draw(""); };
//   auto plot_eta = [&](auto d) { TH1F* htmp = (TH1F*)d.Histo1D({"",";#eta_{Z_{1}};",40,-10,10},"Truth_Dilepton_Eta1","Event_Weight")->Clone(""); htmp->Scale(scale_sig); htmp->Draw(""); };
//   auto plot_phi = [&](auto d) { TH1F* htmp = (TH1F*)d.Histo1D({"",";#phi_{Z_{1}};",99,0,3.14},"Truth_Dilepton_Phi1","Event_Weight")->Clone(""); htmp->Scale(scale_sig); htmp->Draw(""); };
//   auto plot_m   = [&](auto d) { TH1F* htmp = (TH1F*)d.Histo1D({"",";#m_{Z_{1}};",100,0,140},  "Truth_Dilepton_M1",  "Event_Weight")->Clone(""); htmp->Scale(scale_sig); htmp->Draw(""); };
//
//   plot_pt(d_this); gPad->SaveAs(savefile("dists_pt_pre","pdf").c_str());
//   plot_eta(d_this); gPad->SaveAs(savefile("dists_eta_pre","pdf").c_str());
//   plot_phi(d_this); gPad->SaveAs(savefile("dists_phi_pre","pdf").c_str());
//   plot_m(d_this); gPad->SaveAs(savefile("dists_m_pre","pdf").c_str());
//   auto d_that = d_this.Filter("fabs(Truth_Dilepton_M.at(1) - 91.19) < 3");
//   plot_pt(d_that); gPad->SaveAs(savefile("dists_pt_post","pdf").c_str());
//   plot_eta(d_that); gPad->SaveAs(savefile("dists_eta_post","pdf").c_str());
//   plot_phi(d_that); gPad->SaveAs(savefile("dists_phi_post","pdf").c_str());
//   plot_m(d_that); gPad->SaveAs(savefile("dists_m_post","pdf").c_str());
//
   d_uncut.Report()->Print();
   d_bguncut.Report()->Print();
}

