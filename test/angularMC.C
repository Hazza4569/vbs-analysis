void angularMC(string date, string date_input, double mjj=600, double njj=2.5 )
{
   auto readfile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/angularsetup_%g_%g_%s.%s",date.c_str(),mjj,njj,var.c_str(),ft.c_str());
      return string(rtn);
   };
   cout << readfile("parameters","dat") << endl;
   //read in parameters obtained from simulation fits
   ifstream infile(readfile("parameters","dat").c_str());
   double f_L_bkg,f_L_bkg_err,f_p_bkg,f_p_bkg_err,f_L_cmb,f_L_cmb_err,f_p_cmb,f_p_cmb_err,n_bkg,n_tot;
   for (auto var : {&f_L_bkg,&f_L_bkg_err,&f_p_bkg,&f_p_bkg_err,&f_L_cmb,&f_L_cmb_err,&f_p_cmb,&f_p_cmb_err,&n_bkg,&n_tot})
   { double tmp; infile >> tmp; *var = tmp; }

   int nbins(10); //amount used for fitted hists

   //set up inverse transform for random numbers drawn from distribution
   //reduce eqn to for 1 + a1*cos + a2*cos^2
   //related by: a1 = 2(1+f_L)/(A(1-f_L-2f_+)), a2 = (1+f_L)/(1-3f_L).
   TF1 f_gen("f_gen","[0] - (3./(6+2*[2])) * x * ( 1 + [1]*x/2. + [2]*x**2/3. ) - (3./(6+2*[2]))*(1-[1]/2. + [2]/3.)",-1,1);
   //[0] = random number from uniform dist, [1] = a1, [2] = a2.
   double a1 = 2*0.16*(1-f_L_cmb-2*f_p_cmb)/(1+f_L_cmb);
   double a2 = (1-3*f_L_cmb)/(1+f_L_cmb);
   double params[] = {0.1,a1,a2};
   ROOT::Math::WrappedTF1 wf_gen(f_gen);
   wf_gen.SetParameters(params);

   ROOT::Math::BrentRootFinder brf; 
   TRandom2 gen(1);

   auto transform_rndm = [&brf,&wf_gen,&f_gen,&params](double r)
   {
      params[0] = r;
      wf_gen.SetParameters(params);
      brf.SetFunction(wf_gen,-1,1);
      brf.Solve();
      return brf.Root();
   };

   //set up fit for signal amongst background
   auto f_combination = new TF1("f_combination",
         "[0] * (\
         (1-[4])*(3.*( [5]*2.*(1.-x**2) + [6]*(1.-2.*[3]*x+x**2) + (1.-[5]-[6])*(1.+2.*[3]*x+x**2) )/8.)\
         + [4]*(3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + (1.-[1]-[2])*(1.+2.*[3]*x+x**2) )/8.))",
         -1,1);
   f_combination->SetParameters(1,.3,.3,1,1);
   f_combination->FixParameter(3,.16);
   f_combination->FixParameter(5,f_L_bkg);
   f_combination->FixParameter(6,f_p_bkg);

   //parameter limits: allow some deviation, all should really be between 0 and 1.
   f_combination->SetParLimits(1,0,1);
   f_combination->SetParLimits(2,0,1);
   f_combination->SetParLimits(4,0,1);

   //book histograms for parameters of interest
   auto h_fL = new TH1F("h_fL",";f_{L};Events",1000,-5,5);
   auto h_fp = new TH1F("h_fp",";f_{+};Events",1000,-5,5);
   auto h_fs = new TH1F("h_fs",";f_{sig};Events",1000,-5,5);

   //main MC loop for toy experiments
   int n_toys = 10000;
   for (int i_toy = 0; i_toy < n_toys; i_toy++)
   {
      auto h_toy = new TH1F("h_toy","",nbins,-1,1);

      //draw number of events recorded in this toy experiment
      int n_events = gen.Poisson(n_tot);
      //printf("e: %d\n",n_events);

      //draw events from pol. distribution
      for(int i_event = 0; i_event < n_events; i_event++) h_toy->Fill( transform_rndm(gen.Rndm()) );

      //fitting toy histogram
      //new TCanvas();
      f_combination->FixParameter(0,2.*h_toy->Integral()/nbins);
      h_toy->Fit("f_combination","Q");
      //h_toy->Draw("histe");
      //f_combination->Draw("same");
      //h_toy->SetMinimum(0);

      //fill parameter histos
      h_fL->Fill(f_combination->GetParameter(1));
      h_fp->Fill(f_combination->GetParameter(2));
      h_fs->Fill(f_combination->GetParameter(4));

      delete h_toy;
   }

   new TCanvas();
   h_fL->Draw();

   new TCanvas();
   h_fp->Draw();

   new TCanvas();
   h_fs->Draw();
}
