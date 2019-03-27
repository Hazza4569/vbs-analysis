void stacked(std::string date)
{
   gStyle->SetOptStat(0);
   //TColor::InvertPalette();
   gStyle->SetLabelSize(.034,"xyz");
   gStyle->SetTitleSize(.036,"xyz");
   gStyle->SetLabelFont(42,"xyz");
   gStyle->SetCanvasDefW(900); 
   gStyle->SetCanvasDefH(800); 
   gStyle->SetTitleOffset(1.1,"y"); 
   gStyle->SetTitleOffset(1,"x"); 
   gStyle->SetLineWidth(2);
   TLatex ltx;
   ltx.SetTextSize(.03);
   ltx.SetTextFont(42);

   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/significance_%s.%s",date.c_str(),var.c_str(),ft.c_str());
      return string(rtn);
   };

   //setup files:
   double target_luminosity = 35.9; //As CMS analysis
   //setup files:
   std::string pre("/home/user108/y4p/root_output/"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_1M";
   string strBkg = "inclusive_ATLAS_5M";
   string strBkg2 = "ZZgg_ATLAS_600K";

   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());
   ROOT::RDataFrame d_bkg2("EventTree",(pre+strBkg2+post).c_str());

   double scale_sig = target_luminosity *
      (*d_sig.Take<double>("Cross_Section").begin()) / 
      (*d_sig.Take<float>("Weights_Sum").begin()) ;
//     (*d_sig.Take<int>("Event_Count").begin()) ;

   double scale_bkg = target_luminosity *
      (*d_bkg.Take<double>("Cross_Section").begin()) / 
      (*d_bkg.Take<float>("Weights_Sum").begin()) ;
//      (*d_bkg.Take<int>("Event_Count").begin()) ;

   double scale_bkg2 = target_luminosity *
      (*d_bkg2.Take<double>("Cross_Section").begin()) / 
      (*d_bkg2.Take<float>("Weights_Sum").begin()) ;
//      (*d_bkg2.Take<int>("Event_Count").begin()) ;

   auto d_sig_tmp = d_sig.Filter("true");
   auto d_bkg_tmp = d_bkg.Filter("true");
   auto d_bkg2_tmp = d_bkg2.Filter("true");

   auto getevents = [](TH1* h){ return h->Integral()+h->GetBinContent(0)+h->GetBinContent(h->GetNbinsX()); };

   auto *d_curr = &d_sig_tmp;

   int bins(1); double rmin(0), rmax(1); char histstr[100];
   auto make_hist = [&](string var, int bins, double rmin, double rmax, string xlabel, string units=""){
      char ttl[500]; bool nu=(units=="");
      snprintf(ttl,500,";%s %s%s%s; Events / %g %s",xlabel.c_str(),nu?"":"[",units.c_str(),nu?"":"]",
            (rmax-rmin)/bins,units.c_str());

      TH1F* h = (TH1F*)d_curr->Histo1D({"",ttl,bins,rmin,rmax},var,"Event_Weight")->Clone("");
      h->SetLineWidth(1); h->SetLineColor(1);
      return h;
   };

   TH1F* h_sig_M = make_hist("Dijet_M",15,100,1600,"m_{jj}","GeV"); h_sig_M->Scale(scale_sig);
   TH1F* h_sig_e = make_hist("Jet12_Eta_Diff",14,0,7,"|#Delta#eta_{jj}|",""); h_sig_e->Scale(scale_sig);

   d_curr = &d_bkg_tmp;
   TH1F* h_bkg_M = make_hist("Dijet_M",15,100,1600,"m_{jj}","GeV"); h_bkg_M->Scale(scale_bkg);
   TH1F* h_bkg_e = make_hist("Jet12_Eta_Diff",14,0,7,"|#Delta#eta_{jj}|",""); h_bkg_e->Scale(scale_bkg);

   d_curr = &d_bkg2_tmp;
   TH1F* h_bkg2_M = make_hist("Dijet_M",15,100,1600,"m_{jj}","GeV"); h_bkg2_M->Scale(scale_bkg2);
   TH1F* h_bkg2_e = make_hist("Jet12_Eta_Diff",14,0,7,"|#Delta#eta_{jj}|",""); h_bkg2_e->Scale(scale_bkg2);

   //printf("sig: %f bg1: %f bg2: %f\n",getevents(h_sig_e),getevents(h_bkg_e),getevents(h_bkg2_e));

   auto stack_hists = [](string name, string label, TH1* hsig, TH1* hbg, TH1* hbg2){
      auto hs = new THStack(name.c_str(),label.c_str());

      //hbg->Add(hsig,-1);
      //for (int iBin = 1; iBin <= hbg->GetXaxis()->GetNbins(); iBin++) if (hbg->GetBinContent(iBin) < 0) hbg->SetBinContent(iBin,0);

      auto cOrange = new TColor(1411,1,127./150,0);
      hsig->SetFillColor(kRed-4);
      hbg->SetFillColor(cOrange->GetNumber());
      hbg2->SetFillColor(kGreen-4);
      hsig->SetLineColorAlpha(kBlack,1);
      hbg->SetLineColorAlpha(kBlack,1);
      hbg2->SetLineColorAlpha(kBlack,1);
      hsig->SetLineWidth(2);
      hbg->SetLineWidth(2);
      hbg2->SetLineWidth(2);

      hs->Add(hbg2);
      hs->Add(hbg);
      hs->Add(hsig);
      return hs;
   };

   new TCanvas();
   gPad->SetTickx(); gPad->SetTicky();
   gPad->SetLeftMargin(0.1);

   auto h_ewk = stack_hists("ZZjj Dijet Mass",";m_{jj} [GeV];Events / 100 GeV",h_sig_M,h_bkg_M,h_bkg2_M);
   h_ewk->Draw("hist"); gPad->SetLogy(1); h_ewk->SetMaximum(30);
   gPad->GetCanvas()->Update();

   ltx.SetTextSize(.044);
   ltx.SetTextAlign(12);
   ltx.DrawLatex(230,28,"#bf{ #it{Electroweak Selection}}");
   ltx.SetTextSize(.036);
   ltx.DrawLatex(250,19.5,"Dijet mass");
   ltx.DrawLatex(250,13,"#sqrt{s} = 13 TeV""   ""#intLdt = 35.9 fb^{-1}");

   double legx(.62/*.25*/),legy(.7);
   auto leg = new TLegend(legx,legy,legx+.45,legy+.15);
   leg->AddEntry(h_sig_M,"EW ZZjj","f");
   leg->AddEntry(h_bkg_M,"qq #rightarrow ZZ","f");
   leg->AddEntry(h_bkg2_M,"gg #rightarrow ZZ","f");
   leg->SetFillColorAlpha(kWhite, 0);
   leg->SetLineColorAlpha(kBlack, 0);
   leg->SetTextSize(.036);
   leg->Draw();

   gPad->SaveAs(savefile("EWK_Dijet_M","pdf").c_str());

   new TCanvas();
   gPad->SetTickx(); gPad->SetTicky();
   gPad->SetLeftMargin(0.1);

   auto h_vbs = stack_hists("",";|#Delta#eta_{jj}|;Events / 0.5",h_sig_e,h_bkg_e,h_bkg2_e);
   h_vbs->Draw("hist"); h_vbs->SetMaximum(8.5);
   gPad->GetCanvas()->Update();

   ltx.SetTextSize(.044);
   ltx.SetTextAlign(12);
   ltx.DrawLatex(0.32,8.3,"#bf{ #it{Electroweak Selection}}");
   ltx.SetTextSize(.036);
   ltx.DrawLatex(.4,7.83,"Jet Pseudorapidity Difference");
   ltx.DrawLatex(.4,7,"#sqrt{s} = 13 TeV""   ""#intLdt = 35.9 fb^{-1}");

   auto leg2 = new TLegend(legx,legy,legx+.45,legy+.15);
   leg2->AddEntry(h_sig_e,"EW ZZjj","f");
   leg2->AddEntry(h_bkg_e,"qq #rightarrow ZZ","f");
   leg2->AddEntry(h_bkg2_e,"gg #rightarrow ZZ","f");
   leg2->SetFillColorAlpha(kWhite, 0);
   leg2->SetLineColorAlpha(kBlack, 0);
   leg2->SetTextSize(.036);
   leg2->Draw();

   gPad->SaveAs(savefile("EWK_DelEtJJ","pdf").c_str());
}
