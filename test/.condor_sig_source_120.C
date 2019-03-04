void _condor_sig_source_120(string date, string pMethod = "poisson_seq", double target_luminosity=35.9)
{
   auto calcsig = [](double n_s, double n_b, string method = "s_over_root_b",bool verbose=false)
   {
      if (method == "s_over_root_b") return std::pair<double,double>(n_s/sqrt(n_b),-1);
      if (method == "poisson_seq")
      {
         if (verbose) printf("Sequential poisson significance calculation:\n");
         double n_tot = n_s+n_b;
         long double mean(0), previ(-2), pois_sum(0);
         for (int i = 1; pois_sum < 0.9 || fabs(mean-previ) > 0.001 || mean < 0; i++)
         {
            previ = mean;
            double pois = ROOT::Math::poisson_pdf(i,n_tot);
            mean +=
               ROOT::Math::gaussian_quantile(0.5*(1+ROOT::Math::poisson_cdf(i-1,n_b)),1)
               * pois;
            pois_sum += pois;
            if (verbose) printf("%.4Lf\n",mean);
            if (mean > INT_MAX) return std::pair<double,double>(previ,-1.);
         }
         return std::pair<double,double>(mean,-1);
      }
      if (method == "poisson_MC")
      {
         int n_toys = 10000000;
         double rmin = -10;
         double rmax = 10;
         int bins = 10000;
         TH1F *toy_significances = new TH1F("toy_significances","",bins,rmin,rmax);
         TRandom2 gen(1);
         for(int i = 0; i < n_toys; i++)
         {
            int test_statistic_i = gen.Poisson(n_s+n_b); //minus one needed so that p-value is probability of greater than *or equal to*
                                                           //otherwise the cdf function will exclude the current value
            double p_value_i = (test_statistic_i==0)?1:ROOT::Math::poisson_cdf_c(test_statistic_i-1,n_b);
            double significance_i = ROOT::Math::gaussian_quantile(1-.5*p_value_i,1);
            toy_significances->Fill(significance_i);
         }
         toy_significances->Draw();
         return std::pair<double,double>(toy_significances->GetMean(),toy_significances->GetStdDev());
      }
      return std::pair<double,double>(-1.,-1.);
   };

   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/significance_%.1f_%s.%s",date.c_str(),target_luminosity,var.c_str(),ft.c_str());
      return string(rtn);
   };

   //setup files:
   //double target_luminosity = 35.9; //As CMS analysis
   std::string pre("/home/user108/y4p/root_output/"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_1M";
   string strBkg = "inclusive_ATLAS_5M";

   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());

   double scale_sig = target_luminosity *
      (*d_sig.Take<double>("Cross_Section").begin()) / 
      (*d_sig.Take<int>("Event_Count").begin()) ;

   double scale_bkg = target_luminosity *
      (*d_bkg.Take<double>("Cross_Section").begin()) / 
      (*d_bkg.Take<int>("Event_Count").begin()) ;

   //define grid:
   double mjj_min = 100;
   double mjj_max = 2100;//1600;
   double mjj_stp = 50;
   int mjj_n = double(mjj_max-mjj_min)/mjj_stp + 1;

   double njj_min = 0;
   double njj_max = 6;
   double njj_stp = 0.25;
   int njj_n = double(njj_max-njj_min)/njj_stp + 1;

   TH2F* significance_map = new TH2F("significance_map",";Required m_{jj} [GeV];Required |#Delta#eta_{jj}|;Significance",
         mjj_n, mjj_min-.5*mjj_stp, mjj_max+.5*mjj_stp,
         njj_n, njj_min-.5*njj_stp, njj_max+.5*njj_stp);
   TH2F* n_map_sig        = new TH2F("n_map_sig",";Required m_{jj} [GeV];Required |#Delta#eta_{jj}|;Events",
         mjj_n, mjj_min-.5*mjj_stp, mjj_max+.5*mjj_stp,
         njj_n, njj_min-.5*njj_stp, njj_max+.5*njj_stp);
   TH2F* n_map_bkg        = new TH2F("n_map_bkg",";Required m_{jj} [GeV];Required |#Delta#eta_{jj}|;Events",
         mjj_n, mjj_min-.5*mjj_stp, mjj_max+.5*mjj_stp,
         njj_n, njj_min-.5*njj_stp, njj_max+.5*njj_stp);

   std::vector<TH1F*> stats_limit;
   
   std::vector<int> n_statslim = {5,10,20}; //draw line marking where each of these regions end

   for (int i = 0; i < n_statslim.size(); i++)
      stats_limit.push_back(new TH1F((std::string("stats_limit")+to_string(i)).c_str(),"",mjj_n, mjj_min-.5*mjj_stp, mjj_max+.5*mjj_stp));

   //iterate:
   for (int i=1; i<=mjj_n; i++)
   {
      std::vector<bool> limit_reached(stats_limit.size(),false);
      for (int j=1; j<=njj_n; j++)
      {
         double mjj = mjj_min + (i-1)*mjj_stp;
         double njj = njj_min + (j-1)*njj_stp;

         std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
         int n_sig = *d_sig.Filter(cut).Count();
         int n_bkg = *d_bkg.Filter(cut).Count();

         n_map_sig->SetBinContent(i,j,n_sig);
         n_map_bkg->SetBinContent(i,j,n_bkg);

         for (int k = 0; k < n_statslim.size(); k++)
         {
            if ( !limit_reached[k] && (n_bkg < n_statslim[k] || n_sig < n_statslim[k]) )
            {
               stats_limit[k]->SetBinContent(i,njj-.5*njj_stp);
               limit_reached[k] = true;
            }
         }

         double signif = calcsig(n_sig*scale_sig, n_bkg*scale_bkg, pMethod).first;
         //if (signif < 0) printf("neg: %.0f %.0f\n",mjj,njj);
         //if (njj > 7) printf("%4.0f %2.0f %.4f %.4f  |  %.4f\n",mjj,njj,n_sig*scale_sig,n_bkg*scale_bkg,signif);
         significance_map->SetBinContent(i, j, signif);
      }
      for (int k = 0; k < n_statslim.size(); k++) if (!limit_reached[k]) stats_limit[k]->SetBinContent(i,njj_max+.5*njj_stp);
   }

   significance_map->SaveAs(savefile("2Dmap_recipe","C").c_str());
   n_map_sig->SaveAs(savefile("sig_map_recipe","C").c_str());
   n_map_bkg->SaveAs(savefile("bkg_map_recipe","C").c_str());

   new TCanvas();
   gStyle->SetOptStat(0);
   significance_map->Draw("colz");
   for (int k = 0; k < n_statslim.size(); k++)
   {
      stats_limit[k]->SetLineColor(kBlack);
      stats_limit[k]->SetLineWidth(1);
      stats_limit[k]->Draw("histsame");
   }
   gPad->SaveAs(savefile("2Dmap","pdf").c_str());

   n_map_sig->Draw("colz"); gPad->SetLogz();
   gPad->SaveAs(savefile("sig_map","pdf").c_str());
   n_map_bkg->Draw("colz"); gPad->SetLogz();
   gPad->SaveAs(savefile("bkg_map","pdf").c_str());
}
