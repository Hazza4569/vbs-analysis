void significance(string date)
{
   FILE *of = fopen((string("/home/user108/y4p/graph_logs/test.txt")).c_str(),"w");
   fprintf(of,"tst\n");
   fclose(of);

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
      
   //define grid:
   double mjj_min = 100;
   double mjj_max = 1600;
   double mjj_stp = 10;
   int mjj_n = double(mjj_max-mjj_min)/mjj_stp + 1;
    
   double njj_min = 0;
   double njj_max = 6;
   double njj_stp = 0.1;
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

   //iterate:
   for (int i=1; i<=mjj_n; i++)
      for (int j=1; j<=njj_n; j++)
      {
         double mjj = mjj_min + (i-1)*mjj_stp;
         double njj = njj_min + (j-1)*njj_stp;

         std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
         int n_sig = *d_sig.Filter(cut).Count();
         int n_bkg = *d_bkg.Filter(cut).Count();

         n_map_sig->SetBinContent(i,j,n_sig);
         n_map_bkg->SetBinContent(i,j,n_bkg);
         significance_map->SetBinContent(i, j,  (n_sig*scale_sig)/sqrt(n_bkg*scale_bkg)   );
      }
   
   significance_map->SaveAs(savefile("2Dmap_recipe","C").c_str());
   n_map_sig->SaveAs(savefile("sig_map_recipe","C").c_str());
   n_map_bkg->SaveAs(savefile("bkg_map_recipe","C").c_str());

   new TCanvas();
   significance_map->Draw("colz");
   gPad->SaveAs(savefile("2Dmap","pdf").c_str());

   n_map_sig->Draw("colz"); gPad->SetLogz();
   gPad->SaveAs(savefile("sig_map","pdf").c_str());
   n_map_bkg->Draw("colz"); gPad->SetLogz();
   gPad->SaveAs(savefile("bkg_map","pdf").c_str());
}
