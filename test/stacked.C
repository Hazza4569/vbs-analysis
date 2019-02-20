void stacked(std::string date)
{
   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/significance_%s.%s",date.c_str(),var.c_str(),ft.c_str());
      return string(rtn);
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

   auto d_sig_tmp = d_sig.Filter("true");
   auto d_bkg_tmp = d_bkg.Filter("true");

   auto *d_curr = &d_sig_tmp;

   int bins(1); double rmin(0), rmax(1); char histstr[100];
   auto make_hist = [&](string var, int bins, double rmin, double rmax, string xlabel, string units=""){
      char ttl[500]; bool nu=(units=="");
      snprintf(ttl,500,";%s %s%s%s; Events / %g %s",xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
            (rmax-rmin)/bins,units.c_str());

      TH1F* h = (TH1F*)d_curr->Histo1D({"",ttl,bins,rmin,rmax},var)->Clone("");
      h->SetLineWidth(1); h->SetLineColor(1);
      return h;
   };

   TH1F* h_sig_M = make_hist("Dijet_M",15,100,1600,"m_{jj}","GeV"); h_sig_M->Scale(scale_sig);
   TH1F* h_sig_e = make_hist("Jet12_Eta_Diff",14,0,7,"|#Delta#eta_{jj}|",""); h_sig_e->Scale(scale_sig);

   d_curr = &d_bkg_tmp;
   TH1F* h_bkg_M = make_hist("Dijet_M",15,100,1600,"m_{jj}","GeV"); h_bkg_M->Scale(scale_bkg);
   TH1F* h_bkg_e = make_hist("Jet12_Eta_Diff",14,0,7,"|#Delta#eta_{jj}|",""); h_bkg_e->Scale(scale_bkg);

   gStyle->SetOptStat(0);
   gStyle->SetPalette(1,0);
   gStyle->SetLabelSize(.033,"xyz");
   gStyle->SetTitleSize(.038,"xyz");
   gStyle->SetLabelFont(42,"xyz");
   gStyle->SetTitleOffset(1.42,"y"); 
   gStyle->SetCanvasDefW(800); 
   gStyle->SetCanvasDefH(800); 
   //gStyle->SetFrameLineWidth(2);
   //gStyle->SetLineWidth(2);
   //
   TLatex ltx;
   ltx.SetTextSize(.04);
   ltx.SetTextFont(42);

   auto stack_hists = [](string name, string label, TH1* hsig, TH1* hbg){
      auto hs = new THStack(name.c_str(),label.c_str());

      hbg->Add(hsig,-1);
      for (int iBin = 1; iBin <= hbg->GetXaxis()->GetNbins(); iBin++) if (hbg->GetBinContent(iBin) < 0) hbg->SetBinContent(iBin,0);

      hsig->SetFillColor(6);
      hbg->SetFillColor(kAzure-6);
      hbg->SetLineColorAlpha(kBlack,1);
      hsig->SetLineColorAlpha(kBlack,1);

      hs->Add(hbg);
      hs->Add(hsig);
      return hs;
   };

   auto h_ewk = stack_hists("ZZjj Dijet Mass",";m_{jj} [GeV];Events / 100 GeV",h_sig_M,h_bkg_M);
   h_ewk->Draw("hist"); gPad->SetLogy(1); h_ewk->SetMaximum(120);
            gPad->GetCanvas()->Update();
            double ltxX(.64*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin());
            double ltxY(.9*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin());
            ltx.DrawLatex(ltxX,95,"#splitline{#bf{#it{EW Selection}}}{#splitline{X#rightarrow jj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}}");
   gPad->SaveAs(savefile("EWK_Dijet_M","pdf").c_str());

   new TCanvas();
   gPad->SetLeftMargin(0.13);
   auto h_vbs = stack_hists("",";|#Delta#eta_{jj}|;Events / 0.5",h_sig_e,h_bkg_e);
   h_vbs->Draw("hist"); h_vbs->SetMaximum(22);
            gPad->GetCanvas()->Update();
            ltxX = 0.64*(gPad->GetUxmax()-gPad->GetUxmin())+gPad->GetUxmin();
            ltxY = 0.9*(gPad->GetUymax()-gPad->GetUymin())+gPad->GetUymin();
            ltx.DrawLatex(ltxX,ltxY,"#splitline{#bf{#it{EW Selection}}}{#splitline{X#rightarrow jj}{#splitline{L_{int} = 35.9 fb^{-1}}{#sqrt{s} = 13 TeV}}}");
   gPad->SaveAs(savefile("EWK_DelEtJJ","pdf").c_str());
}
