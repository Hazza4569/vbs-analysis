void histGen(std::string file)
{
   ROOT::RDataFrame d("EventTree",file);
   double target_luminosity = 150;//inv fb

   std::string nfile = file.substr(0,file.find("_o"));
   std::string path("graphs/"); file = path+nfile;
   int n_events = *(d.Take<int>("Event_Count").begin());
   double cross_section = *(d.Take<double>("Cross_Section").begin());
   Double_t scale = target_luminosity*cross_section/n_events;

   TH1D* h_mu = (TH1D*)d.Histo1D({"",(nfile+std::string(";m_{#mu#mu} [GeV];Events / Gev")).c_str(),141,0,140},"m_mumu")->Clone("dimuon");
   h_mu->Scale(scale);
   h_mu->Draw("hist");
   gPad->SaveAs((file+std::string("_mu.pdf")).c_str());

   TH1D *h_e = (TH1D*)d.Histo1D({"",(nfile+std::string(";m_{ee} [GeV];Events / Gev")).c_str(),141,0,140},"m_ee")->Clone("dielectron");
   h_e->Scale(scale);
   h_e->Draw("hist");
   gPad->SaveAs((file+std::string("_e.pdf")).c_str());

   TH1D *h_l = (TH1D*)d.Histo1D({"",(nfile+std::string(";m_{ll} [GeV];Events / Gev")).c_str(),141,0,140},"m_ll")->Clone("dilepton");
   h_l->Scale(scale);
   h_l->Draw("hist");
   gPad->SaveAs((file+std::string("_l.pdf")).c_str());
   
   auto d_new = d.Define("m1","(Lepton_Pairs>0)?m_ll.at(0):-3.14").Define("m2","(Lepton_Pairs>1)?m_ll.at(1):-3.14")
                 .Define("m1_sq","pow(m1,2)").Define("m2_sq","pow(m2,2)");
   auto d_two = d_new.Filter("Lepton_Pairs>1");

   TH2D *h_2d = (TH2D*)d_two.Histo2D({"",(nfile+std::string(";m_{ll} [GeV];m_{l'l'} [GeV]")).c_str(),141,0,140,141,0,140},"m1","m2")->Clone("2D dilepton");
   h_2d->Scale(scale);
   h_2d->Draw("colz");
   gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
   gPad->SaveAs((file+std::string("_2d.pdf")).c_str());

   TH2D *h_2d_sq = (TH2D*)d_two.Histo2D({"",(nfile+std::string(";m_{ll}^{2} [GeV^{2}];m_{l'l'}^{2} [GeV^{2}]")).c_str(),141,0,19600,141,0,18600},"m1_sq","m2_sq")->Clone("2D dilepton squared");
   h_2d_sq->Scale(scale);
   h_2d_sq->Draw("colz");
   gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
   gPad->SaveAs((file+std::string("_2d_sq.pdf")).c_str());
}
