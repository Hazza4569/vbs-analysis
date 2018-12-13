void jet_check(string date)
{
   for (std::string file : {"ZZjj_ATLAS_500K","inclusive_ATLAS_500K"})
   {
      double target_luminosity = 150;
      std::string pre("/home/user108/y4p/root_output/"),post("_o.root");

      ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());
      int n_events = *(d.Take<int>("Event_Count").begin());
      double cross_section = *(d.Take<double>("Cross_Section").begin());
      Double_t scale = target_luminosity*cross_section/n_events;

      int bins(1); double rmin(0), rmax(1); char histstr[100];
      auto histname = [&bins,&rmin,&rmax,&file](string xlabel, string units){
         char rtn[500]; bool nu=(units=="");
         snprintf(rtn,500,"%s;%s %s%s%s; Events / %g %s",file.c_str(),xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
               (rmax-rmin)/bins,units.c_str());
         return rtn;
      };

      auto d_new = d.Filter("Jet_Good_n > 1"); 

      rmin = 0; rmax = 600; bins = 80; 
      TH1D* hPtjj = (TH1D*)d_new.Histo1D({"",histname("p_{T,jj}","GeV"),bins,rmin,rmax},"Dijet_Pt")->Clone("Dijet_Pt");
      rmin = -8; rmax = 8; bins = 100; 
      TH1D* hEtajj = (TH1D*)d_new.Histo1D({"",histname("#eta_{jj}",""),bins,rmin,rmax},"Dijet_Eta")->Clone("Dijet_Eta");
      rmin = -3.5; rmax = 3.5; bins = 70; 
      TH1D* hPhijj = (TH1D*)d_new.Histo1D({"",histname("#phi_{jj}","rad"),bins,rmin,rmax},"Dijet_Phi")->Clone("Dijet_Phi");
      rmin = 0; rmax = 2000; bins = 100; 
      TH1D* hMjj = (TH1D*)d_new.Histo1D({"",histname("m_{jj}","GeV"),bins,rmin,rmax},"Dijet_M")->Clone("Dijet_M");
      rmin = -5; rmax = 5; bins = 50; 
      TH1D* hYjj = (TH1D*)d_new.Histo1D({"",histname("y_{jj}",""),bins,rmin,rmax},"Dijet_Rapidity")->Clone("Dijet_Rapidity");
      rmin = 0; rmax = 500; bins = 80; 
      TH1D* hPtj = (TH1D*)d_new.Histo1D({"",histname("p_{T,j}","GeV"),bins,rmin,rmax},"Jet_Pt")->Clone("Jet_Pt");
      rmin = -5; rmax = 5; bins = 100; 
      TH1D* hEtaj = (TH1D*)d_new.Histo1D({"",histname("#eta_{j}",""),bins,rmin,rmax},"Jet_Eta")->Clone("Jet_Eta");
      rmin = -3.5; rmax = 3.5; bins = 70; 
      TH1D* hPhij = (TH1D*)d_new.Histo1D({"",histname("#phi_{j}","rad"),bins,rmin,rmax},"Jet_Phi")->Clone("Jet_Phi");
      rmin = 0; rmax = 100; bins = 100; 
      TH1D* hMj = (TH1D*)d_new.Histo1D({"",histname("m_{j}","GeV"),bins,rmin,rmax},"Jet_M")->Clone("Jet_M");
      rmin = -5; rmax = 5; bins = 100; 
      TH1D* hYj = (TH1D*)d_new.Histo1D({"",histname("y_{j}",""),bins,rmin,rmax},"Jet_Rapidity")->Clone("Jet_Rapidity");
      rmin = 0; rmax = 10; bins = 100; 
      TH1D* hDR = (TH1D*)d.Histo1D({"",histname("#Delta R_{j}",""),bins,rmin,rmax},"Jet_DeltaR")->Clone("Jet_DeltaR");

      for (auto &h : {hPtjj,hEtajj,hPhijj,hMjj,hYjj,hPtj,hEtaj,hPhij,hMj,hYj,hDR})
      {
         h->Scale(scale);
         h->Draw("hist");
         gPad->SaveAs((string("/home/user108/y4p/graph_logs/")+date+string("/")+string(h->GetName())+file+string(".pdf")).c_str());
      }
   }
}
