void weights_test()
{
   std::string pre("/home/user108/y4p/root_output/"),post("_o.root");
   std::string file("ZZjj_ATLAS_500K");

   ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());

   double target_luminosity = 35.9;
   int n_events = d.Histo1D("Cross_Section")->GetEntries();//*(d.Take<int>("Event_Count").begin());
   double cross_section = *(d.Take<double>("Cross_Section").begin());
   double unweighted_scale = target_luminosity*cross_section/n_events;

   double weightsum = d.Sum("Event_Weight").GetValue();
   double scale = n_events*unweighted_scale/weightsum;
      auto d_new = d.Define("Lead_Lepton_Pt","Isolated_Lepton_Pt.size() < 1 ? -1 : Isolated_Lepton_Pt.at(0)")
         .Define("Sublead_Lepton_Pt","Isolated_Lepton_Pt.size() < 2 ? -1 : Isolated_Lepton_Pt.at(1)")
         .Define("Sublead_Jet_Pt","Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].size() < 2 ? -1 : Jet_Pt[Jet_DeltaR>=Jet_DeltaR_Min&&Jet_Flav<Jet_Flav_Max].at(1)");
   auto d_filtered = d_new.Filter("Lepton_Pairs>1");

  //Histogram with weights 
   new TCanvas();
   TH1F *hWeight = (TH1F*)d_filtered.Histo1D(/*{"",";m_{ll} [GeV]; Events / 2 GeV",5000,-500,1000},*/"Muon_Pt","Event_Weight")
      ->Clone("hWeight");
   hWeight->Scale(scale);
   hWeight->Draw("hist");
}
