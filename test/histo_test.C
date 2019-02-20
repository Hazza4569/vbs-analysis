void histo_test()
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


   auto d_filtered = d.Filter("Lepton_Pairs>1");

  //Histogram with no scaling
   //TH1F *h = new TH1F("h","",,,);
   new TCanvas();
   TH1F *h = (TH1F*)d_filtered.Histo1D({"",";m_{4l} [GeV]; Events / 2 GeV",100,0,200},"Tetralepton_M")->Clone("h");
   h->Draw();

  //Histogram combining
   new TCanvas();
   TH1F *h1 = (TH1F*)d_filtered.Histo1D({"",";m_{ee} [GeV]; Events / 2 GeV",100,0,200},"Dielectron_M")->Clone("h1");
   h1->Scale(unweighted_scale);

   TH1F *h2 = (TH1F*)d_filtered.Histo1D({"",";m_{#mu#mu} [GeV]; Events / 2 GeV",100,0,200},"Dimuon_M")->Clone("h2");
   h2->Scale(unweighted_scale);

   TH1F *hCombined = new TH1F("hCombined",";m_{ll} [GeV]; Events / 2 GeV",100,0,200);

   hCombined->Add(h1);
   hCombined->Add(h2);
   hCombined->Draw("hist");

  //Histogram stacking
   new TCanvas();
   TH1F *h1b = (TH1F*)d_filtered.Histo1D({"",";m_{ee} [GeV]; Events / 2 GeV",100,0,200},"Dielectron_M")->Clone("h1");
   h1b->Scale(unweighted_scale); h1b->SetFillColor(6); 

   TH1F *h2b = (TH1F*)d_filtered.Histo1D({"",";m_{#mu#mu} [GeV]; Events / 2 GeV",100,0,200},"Dimuon_M")->Clone("h2");
   h2b->Scale(unweighted_scale); h2b->SetFillColor(7);

   THStack *hStacked = new THStack("hStacked",";m_{ll} [GeV]; Events / 2 GeV");

   hStacked->Add(h1b); 
   hStacked->Add(h2b); 
   hStacked->Draw("hist");

  //Histogram with weights 
   new TCanvas();
   TH1F *hWeight = (TH1F*)d_filtered.Histo1D({"",";m_{ll} [GeV]; Events / 2 GeV",100,0,200},"Tetralepton_M","Event_Weight")
      ->Clone("hWeight");
   hWeight->Scale(scale);
   hWeight->Draw("hist");
}
