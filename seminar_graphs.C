void seminar_graphs()
{
   ROOT::RDataFrame d0("EventTree","/home/user108/y4p/root_output/ZZjj_ATLAS_500K_o.root");
   int n_events = d0.Histo1D("Cross_Section")->GetEntries();//*(d.Take<int>("Event_Count").begin());
   double cross_section = *(d0.Take<double>("Cross_Section").begin());
   double target_luminosity = 35.9;
   double unweighted_scale = target_luminosity*cross_section/n_events;

   ROOT::RDataFrame d1("EventTree","/home/user108/y4p/root_output/ZZjj_ATLAS_1M_filtered.root");

   auto d_tmp0 = d1.Filter("true");
   auto d_tmp = d0.Filter("true");
   auto d_curr = &d_tmp;

   int bins(1); double rmin(0), rmax(1); char histstr[100];
   auto make_hist = [&](string var, int bins, double rmin, double rmax, string xlabel, string units=""){
      char ttl[500]; bool nu=(units=="");
      snprintf(ttl,500,";%s %s%s%s; Events / %g %s",xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
            (rmax-rmin)/bins,units.c_str());

      TH1F* h = (TH1F*)d_curr->Histo1D({"",ttl,bins,rmin,rmax},var)->Clone("");
      h->Scale(unweighted_scale); 
      h->SetLineWidth(2); h->SetLineColor(1);
      return h;
   };

   gStyle->SetOptStat(0);
   gStyle->SetPalette(1,0);
   gStyle->SetLabelSize(.04,"xyz");
   gStyle->SetTitleSize(.04,"xyz");
   gStyle->SetLabelFont(42,"xyz");
   //gStyle->SetFrameLineWidth(2);
   //gStyle->SetLineWidth(2);
   //
   TLatex ltx;
   ltx.SetTextSize(.04);
   ltx.SetTextFont(42);

   new TCanvas();
   auto d_pair = d0.Filter("Lepton_Pairs > 0"); d_curr = &d_pair;
   auto h = make_hist("Dilepton_M",140,0,140,"m_{ll}","GeV");
   h->Draw("hist");
   ltx.DrawLatex(10,8.2,"#splitline{Z#rightarrow ll}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}");
   gPad->SaveAs("/home/user108/y4p/seminargraphs/mll1.pdf");

   new TCanvas();
   auto d_twopair = d0.Filter("Lepton_Pairs > 1"); d_curr = &d_twopair;
   h = make_hist("Tetralepton_M",80,50,250,"m_{llll}","GeV");
   h->Draw("hist");
   h->GetYaxis()->SetTitleOffset(1.2);
   ltx.DrawLatex(60,.38,"#splitline{X#rightarrow llll}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}");
   gPad->SaveAs("/home/user108/y4p/seminargraphs/mllll1.pdf");

   new TCanvas();
   auto d_twojet = d0.Filter("Jet_n > 1"); d_curr = &d_twojet;
   h = make_hist("Dijet_M",16,0,1600,"m_{jj}","GeV");
   h->Draw("hist");
   //h->GetYaxis()->SetTitleOffset(1.2);
   ltx.DrawLatex(1200,70,"#splitline{X#rightarrow jj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}");
   gPad->SaveAs("/home/user108/y4p/seminargraphs/mjj1.pdf");

   new TCanvas();
   d_curr = &d_tmp0;
   h = make_hist("Tetralepton_M",20,100,800,"m_{llll}","GeV");
   h->Draw("hist");
   gPad->GetCanvas()->Update();
   double ltxX(.73*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin());
   double ltxY(.84*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin());
   ltx.DrawLatex(ltxX,ltxY,"#splitline{#splitline{#bf{#it{VBS Enhanced}}}{X#rightarrow llll}}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}");
   gPad->SaveAs("/home/user108/y4p/seminargraphs/mllll2.pdf");
}
