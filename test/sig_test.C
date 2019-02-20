void sig_test()
{
   auto calcsig = [](double n_s, double n_b, string method = "s_over_root_b",bool verbose=false)
   {
      if (method == "s_over_root_b") return std::pair<double,double>(n_s/sqrt(n_b),-1);
      if (method == "poisson_seq")
      {
         if (verbose) printf("Sequential poisson significance calculation:\n");
         double n_tot = n_s+n_b;
         long double mean(0), previ(-2);
         for (int i = 1; i <= 20*ceil(n_tot) && (fabs(mean-previ) > 0.01 || mean < 0); i++)
         {
            previ = mean;
            mean +=
               ROOT::Math::gaussian_quantile(0.5*(1+ROOT::Math::poisson_cdf(i-1,n_b)),1)
               * ROOT::Math::poisson_pdf(i,n_tot);
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

   //setup files:
   double target_luminosity = 35.9; //As CMS analysis
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

   double mjj = 1550;
   double njj = 5.5;

   std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
   int n_sig = *d_sig.Filter(cut).Count();
   int n_bkg = *d_bkg.Filter(cut).Count();

   auto mcval = calcsig(n_sig*scale_sig,n_bkg*scale_bkg,"poisson_MC");

   printf("%21s %6s %6s\n%20s: %6d %6d\n%20s: %6e %6e\n\n%20s: %.4f\n%20s: %.4f\n%20s: %.4f +/- %.4f\n",
         //headers
         "","signal","background",
         "Simulated events",n_sig,n_bkg,
         "Scaled events",n_sig*scale_sig,n_bkg*scale_bkg,
         "Approx significance",calcsig(n_sig*scale_sig,n_bkg*scale_bkg,"s_over_root_b").first,
         "Poisson significance",calcsig(n_sig*scale_sig,n_bkg*scale_bkg,"poisson_seq").first,
         "MonteCarlo significance",mcval.first,mcval.second);
}
