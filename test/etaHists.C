void etaHists()
{
   double target_luminosity = 150;
   std::string pre("/home/user108/y4p/root_output/"),post("_o.root");
   std::string output("/home/user108/y4p/tex/eta/graphs/"),out_post(".pdf");
   for (std::string file : {"ZZjj_10K","ZZ_10K","ZZjj_ATLAS_500K","inclusive_ATLAS_500K"})
   {
      ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
      int n_events = *(d.Take<int>("Event_Count").begin());
      double cross_section = *(d.Take<double>("Cross_Section").begin());
      Double_t scale = target_luminosity*cross_section/n_events;
      std::string eta("_Eta");
      for (std::string obj : {"Electron","Muon","Dielectron","Dimuon","Dilepton","Tetralepton",
                              "Truth_Lepton"})
      {
         std::string filter = "true";
         double eta_max = 6;
         int nBins = 61;

         if ( obj == "Electron" || obj == "Muon" ) eta_max = 3;
         else if ( obj == "Tetralepton" )
         {
            eta_max = 0.2;
            filter = "Lepton_Pairs > 1";
         }
         double eta_min = -eta_max;
         char histstr[50];
         snprintf(histstr,50,"%s;#eta;Events / %.2f",file.c_str(),(eta_max-eta_min)/(nBins-1));

         TH1D* h = (TH1D*) d.Filter(filter).Histo1D({"",histstr,nBins,eta_min,eta_max},
                                    obj+eta) ->Clone((obj+eta).c_str());
         h->Scale(scale); 
         h->Draw("hist");
         gPad->SaveAs((output+file+obj+out_post).c_str());
      }
   }
}
