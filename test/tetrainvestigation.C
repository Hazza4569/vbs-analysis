void tetraleptoninvestigation(string file)
{
   double target_luminosity = 150;
   std::string pre("/home/user108/y4p/root_output/"),post("_o.root");

   ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
   int n_events = *(d.Take<int>("Event_Count").begin());
   double cross_section = *(d.Take<double>("Cross_Section").begin());
   Double_t scale = target_luminosity*cross_section/n_events;

   int bins(1); double rmin(0), rmax(1); char histstr[100];
   auto histname = [&bins,&rmin,&rmax](string xlabel, string units){
      char rtn[200]; bool nu=(units=="");
      snprintf(rtn,200,";%s %s%s%s; Events / %g %s",xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
         (rmax-rmin)/bins,units.c_str());
      return rtn;
   };

   auto d_new = d.Filter("Lepton_Pairs > 1"); 
   
   rmin = 0; rmax = 400; bins = 80; 
   TH1D* hPt = (TH1D*)d_new.Histo1D({"",histname("p_{T,4l}","GeV"),bins,rmin,rmax},"Tetralepton_Pt")->Clone("Tetralepton_Pt");
   rmin = -8; rmax = 8; bins = 100; 
   TH1D* hEta = (TH1D*)d_new.Histo1D({"",histname("#eta_{4l}",""),bins,rmin,rmax},"Tetralepton_Eta")->Clone("Tetralepton_Eta");
   rmin = -3.5; rmax = 3.5; bins = 70; 
   TH1D* hPhi = (TH1D*)d_new.Histo1D({"",histname("#phi_{4l}","rad"),bins,rmin,rmax},"Tetralepton_Phi")->Clone("Tetralepton_Phi");
   rmin = 0; rmax = 400; bins = 100; 
   TH1D* hM = (TH1D*)d_new.Histo1D({"",histname("m_{4l}","GeV"),bins,rmin,rmax},"Tetralepton_M")->Clone("Tetralepton_M");
   rmin = -3; rmax = 3; bins = 40; 
   TH1D* hY = (TH1D*)d_new.Histo1D({"",histname("y_{4l}",""),bins,rmin,rmax},"Tetralepton_Rapidity")->Clone("Tetralepton_Rapidity");

   for (auto &h : {hPt,hEta,hPhi,hM,hY})
   {
      h->Scale(scale);
      new TCanvas();
      h->Draw("hist");
   }
}
