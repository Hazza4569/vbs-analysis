void likelihood_significance(string date)
{
   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/significance_%s.%s",date.c_str(),var.c_str(),ft.c_str());
      return string(rtn);
   };

   //setup files:
   double target_luminosity = 35.9; //As CMS analysis
   std::string pre("/home/user108/y4p/root_output/"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_500K";
   string strBkg = "inclusive_ATLAS_500K";

   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());

   double scale_sig = target_luminosity *
                      (*d_sig.Take<double>("Cross_Section").begin()) / 
                      (*d_sig.Take<int>("Event_Count").begin()) ;
      
   double scale_bkg = target_luminosity *
                      (*d_bkg.Take<double>("Cross_Section").begin()) / 
                      (*d_bkg.Take<int>("Event_Count").begin()) ;

   double mjj = 1000; //GeV
   double njj = 3.5;

   std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
   int n_sig = *d_sig.Filter(cut).Count();
   int n_bkg = *d_bkg.Filter(cut).Count();

   double n_scale_sig = n_sig*scale_sig;
   double n_scale_bkg = n_bkg*scale_bkg;
   double n_tot = n_scale_sig + n_scale_bkg;
   printf("===Events===\nSignal:     %.3e\nBackground: %.3e\nTotal: %.3e\n\n",n_scale_sig,n_scale_bkg,n_tot);

   //not sure how best to calculate significance. A few attempts:

   //First: negaitve log likelihood ratio... using Poisson likelihoods
   auto nn_log_likelihood = [](double n_obs, double n_model) { return n_obs*log(n_model) - n_model; };
      //log of Poisson(n_obs; n_model) without factorial divisor (should cancel)

   //taking ratio to be L(mu=1)/L(mu=0) comes out as > 1, resulting in -ve test statistic...
   //makes no sense with chisquare dist. so removing factor -1 for now.
   double test_statistic_nllr = 2*(nn_log_likelihood(n_tot,n_tot) - nn_log_likelihood(n_tot,n_scale_bkg));
   printf("twice (not negative) log likelihood ratio: %.4e\n",test_statistic_nllr);

   //distributed as chisquare w/ 1 dof
   double p_value_nllr = ROOT::Math::chisquared_cdf_c(test_statistic_nllr,1);
   printf("resulting p-value: %.4e\n",p_value_nllr);

   double significance_nllr = ROOT::Math::gaussian_quantile(1-p_value_nllr,1);
   printf("significance: %.4f\n\n",significance_nllr);
   
   //now basic Poisson counting... null hypothesis that events are distributed according to
   //Poisson(n; b). Simply use n_obs as test statistic distributed thusly.
   //Only issue here is that while the mean of a poisson may be non-integer, any test statistic can't be...
   //can't integrate from non-integer limit as the function is not defined there. Ceiling or floor?

   //ceiling would seem to make sense, as the p-value is probability of getting a value 'equally or less likely'
   //to have come from null dist.. Can't do equally so only less.
   int test_statistic_pceil = ceil(n_tot);
   printf("ceiling of total num. events: %d\n",test_statistic_pceil);
   
   double p_value_pceil = ROOT::Math::poisson_cdf_c(test_statistic_pceil,n_scale_bkg);
   printf("resulting p-value: %.4e\n",p_value_pceil);

   double significance_pceil = ROOT::Math::gaussian_quantile(1-p_value_pceil,1);
   printf("significance: %.4f\n\n",significance_pceil);


   //however, analytically speaking a ceiling would always give a better significance than might be
   //realistic. Maybe best to go for a pessimistic estimate? 
   int test_statistic_pfloor = floor(n_tot);
   printf("floor of total num. events: %d\n",test_statistic_pfloor);
   
   double p_value_pfloor = ROOT::Math::poisson_cdf_c(test_statistic_pfloor,n_scale_bkg);
   printf("resulting p-value: %.4e\n",p_value_pfloor);

   double significance_pfloor = ROOT::Math::gaussian_quantile(1-p_value_pfloor,1);
   printf("significance: %.4f\n\n",significance_pfloor);


   //alternatively... hack the poisson by integrating over all values? almost certainly wrong
   //auto nn_poisson_hack = [](double n_obs, double n_model) { return pow(n_model,n_obs) * exp(-n_model); };
   
   
   //and then the basic estimate, which I don't think is really applicable for numbers this small?
   printf("basic [s/sqrt(b)] estimate of significance: %.4f\n",n_scale_sig/sqrt(n_scale_bkg));


   //foun it! need to use the Poisson method but run multiple times with toy experiments which
   //draw a random total num. events from a Poisson centred on our expected total.
   int n_toys = 10000000;
   double rmin = -5;
   double rmax = 10;
   int bins = 10000;
   auto histname = [&rmin,&rmax,&bins](string xlabel, string units){
      char rtn[500]; bool nu=(units=="");
      snprintf(rtn,500,";%s %s%s%s; Events / %g %s",xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
            (rmax-rmin)/bins,units.c_str());
      return string(rtn);
   };
   TH1F *toy_significances = new TH1F("toy_significances",histname("Significance","#sigma").c_str(),bins,rmin,rmax);
   TRandom2 gen(1);
   std::set<double> discrete_sigfs;

   for(int i = 0; i < n_toys; i++)
   {
      int test_statistic_i = gen.Poisson(n_tot);
      double p_value_i = ROOT::Math::poisson_cdf_c(test_statistic_i,n_scale_bkg);
      double significance_i = ROOT::Math::gaussian_quantile(1-p_value_i,1);

      toy_significances->Fill(significance_i);
      discrete_sigfs.insert(significance_i);
   }

   double prev_sig = -100;
   FILE* of = fopen(savefile("dsc_sigs","dat").c_str(),"w");
   for (double sigf : discrete_sigfs)
   {
      fprintf(of,"%.4f,%.4f\n",sigf,sigf-prev_sig);
      prev_sig=sigf;
   }

   toy_significances->Draw();
   gPad->SaveAs(savefile("optimum_significance_dist","pdf").c_str());
   printf("\nPoisson method with toy experiments:\nExpected significance = %.4f +/- %.4f\n",toy_significances->GetMean(),toy_significances->GetStdDev());


   //approximation by converting each integer of the poisson:
   //const int max_int = 18;
   int maximal = (n_tot > 10) ? ceil(2*n_tot) : ceil(10*n_tot);
   for (int max_int = max(maximal-30,0); max_int < maximal; max_int++)
   {
      long double mean = 0;
      double sum = 0;
      long double arr_sigfs[max_int+1];
      for (int i = 0; i <= max_int; i++)
      {
         arr_sigfs[i] = ROOT::Math::gaussian_quantile(ROOT::Math::poisson_cdf(i,n_scale_bkg),1);
         mean +=
            //value of significance for this integer:
            arr_sigfs[i]
            //probability of obtaining this integer:
            * ROOT::Math::poisson_pdf(i,n_tot);
      }
      for (int i = 0; i <= max_int; i++) sum += pow(arr_sigfs[i]-mean,2);
      double std_dev = sqrt(sum)/max_int;

      printf("%d term approximation = %.4Lf +/- %.4f (or %.4Lf)\n",max_int+1,mean,std_dev,sqrt(mean));
   }

   //further testing:
   auto calcsig = [](double n_s, double n_b, string method = "s_over_root_b",bool verbose=false)
   {
      if (method == "sig_over_root_b") return n_s/sqrt(n_b);
      if (method == "poisson_seq")
      {
         if (verbose) printf("Sequential poisson significance calculation:\n");
         double n_tot = n_s+n_b;
         long double mean(0), previ;
         for (int i = 0; i <= 20*ceil(n_tot) && (fabs(mean-previ) > 0.01 || mean < 0); i++)
         {
            previ = mean;
            mean +=
               ROOT::Math::gaussian_quantile(ROOT::Math::poisson_cdf(i,n_b),1)
               * ROOT::Math::poisson_pdf(i,n_tot);
            if (verbose) printf("%.4Lf\n",mean);
            if (mean > INT_MAX) return double(previ);
         }
         return double(mean);
      }
      return -1.;
   };

   printf("result: %f\n",calcsig(n_scale_sig,n_scale_bkg,"poisson_seq",true)); 
}
