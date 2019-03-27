void angularMC_setup(string date, double mjj=600, double njj=2.5)
{
   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/angularsetup_%g_%g_%s.%s",date.c_str(),mjj,njj,var.c_str(),ft.c_str());
      return string(rtn);
   };

   //define cuts for this investigation
   char cuts[1000];
   snprintf(cuts,1000,"Truth_Lepton_Pairs > 1 && fabs(Truth_Dilepton_M.at(1) - 91.19) < 5 && Jet_n > 1 && Jet_Pt.at(1) > 20 && Dijet_M > %f && Jet12_Eta_Diff > %f && fabs(Event_Weight) < 40",mjj,njj);
   
   //Load signal and background files 
   double target_luminosity = 3000; //fb^-1
   string sigfile("ZZjj_ATLAS_1M"),bkgfile("inclusive_ATLAS_5M");
   string pre(/*"/home/user108/y4p/root_output/"*/"/scratch/rootfiles/"),post("_angle.root");
   ROOT::RDataFrame d_sig_init("EventTree",pre+sigfile+post); auto d_sig = d_sig_init.Filter("fabs(Event_Weight) < 40");
   ROOT::RDataFrame d_bkg_init("EventTree",pre+bkgfile+post); auto d_bkg = d_bkg_init.Filter("fabs(Event_Weight) < 40");

   double scale_sig = target_luminosity *
      (*d_sig.Take<double>("Cross_Section").begin()) / 
      d_sig.Sum("Event_Weight").GetValue() ;
   double scale_bkg = target_luminosity *
      (*d_bkg.Take<double>("Cross_Section").begin()) / 
      d_bkg.Sum("Event_Weight").GetValue() ;
//perform cuts
   auto d_sig_cut = d_sig.Define("Z_weight","ROOT::VecOps::RVec<float>(2,Event_Weight)").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
                         .Filter(cuts,"sigcuts");
   auto d_bkg_cut = d_bkg.Define("Z_weight","ROOT::VecOps::RVec<float>(2,Event_Weight)").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
                         .Filter(cuts,"bkgcuts");
   
   //extract histograms, put into background and combined histograms.
   int nbins(10); //amount used for fitted hists
   auto histtitle = [](int bins, string htitle)
   {
      char rtn[500];
      snprintf(rtn,500,"%s;cos#theta*;Events / %g",htitle.c_str(),2./bins);
      return string(rtn);
   };
   auto h_bkg = (TH1F*)d_bkg_cut.Histo1D({"",histtitle(nbins,"background").c_str(),nbins,-1,1},"Truth_Costhetastar","Z_weight")->Clone("h_bkg");
   h_bkg->Scale(scale_bkg);
   auto h_sig = (TH1F*)d_sig_cut.Histo1D({"",histtitle(nbins,"signal").c_str(),nbins,-1,1},"Truth_Costhetastar","Z_weight")->Clone("h_cmb");
   h_sig->Scale(scale_sig);
   printf("S: %e, B: %e, S/B: %e\n",h_sig->Integral(),h_bkg->Integral(),h_sig->Integral()/h_bkg->Integral());
   TH1F* h_cmb = new TH1F("",histtitle(nbins,"combined").c_str(),nbins,-1,1);
   h_cmb->Add(h_bkg);
   h_cmb->Add(h_sig);
   THStack *hstak = new THStack("",histtitle(nbins,"stacked").c_str());
   hstak->Add(h_bkg);
   hstak->Add(h_sig);

   //drawing styles:
   gStyle->SetStatX(0.7); gStyle->SetStatY(0.55);
   gStyle->SetOptStat(1111110); gStyle->SetOptFit(1);

   //fit background and extract parameters (and errors - feeds systematics) [INCLUDES PLOT TO AID DEBUGGING/UNDERSTANDING]
   auto f_simple = new TF1("f_simple","[0]*3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + (1.-[1]-[2])*(1.+2.*[3]*x+x**2) )/8.",-1,1);
   //auto f_simple = new TF1("f_simple","[0]*3.*( [1]*2.*(1.-x**2) + (1-[1])*[2]*(1.-2.*[3]*x+x**2) + (1-[1])*(1-[2])*(1.+2.*[3]*x+x**2) )/8.",-1,1);
   f_simple->SetParameters(1,.3,.3,1); 
   f_simple->FixParameter(3,.16);
   
   f_simple->FixParameter(0,2.*h_bkg->Integral()/nbins);
   new TCanvas();
   h_bkg->Fit("f_simple","Q");
   h_bkg->Draw("histe");
   f_simple->Draw("same");
   h_bkg->SetMinimum(0);
   gPad->SaveAs(savefile("background","pdf").c_str());
   double f_L_bkg = f_simple->GetParameter(1);
   double f_L_bkg_err = f_simple->GetParError(1);
   double f_p_bkg = f_simple->GetParameter(2);
   double f_p_bkg_err = f_simple->GetParError(2);

   //fit signal and extract parameters (and errors? feed systematics)
   auto f_simple1 = new TF1("f_simple1","f_simple",-1,1);
   f_simple1->SetParameters(1,.3,.3,1); 
   f_simple1->FixParameter(3,.16);
   f_simple1->FixParameter(0,2.*h_sig->Integral()/nbins);
   new TCanvas();
   h_sig->Fit("f_simple1","Q");
   h_sig->Draw("histe");
   f_simple1->Draw("same");
   h_sig->SetMinimum(0);
   gPad->SaveAs(savefile("signal","pdf").c_str());
   double f_L_sig = f_simple1->GetParameter(1);
   double f_L_sig_err = f_simple1->GetParError(1);
   double f_p_sig = f_simple1->GetParameter(2);
   double f_p_sig_err = f_simple1->GetParError(2);

   //fit combined and extract parameters (and errors? feed systematics)
   auto f_simple2 = new TF1("f_simple2","f_simple",-1,1);
   f_simple2->SetParameters(1,.3,.3,1); 
   f_simple2->FixParameter(3,.16);
   f_simple2->FixParameter(0,2.*h_cmb->Integral()/nbins);
   new TCanvas();
   h_cmb->Fit("f_simple2","Q");
   h_cmb->Draw("histe");
   f_simple2->Draw("same");
   h_cmb->SetMinimum(0);
   gPad->SaveAs(savefile("combination","pdf").c_str());
   double f_L_cmb = f_simple2->GetParameter(1);
   double f_L_cmb_err = f_simple2->GetParError(1);
   double f_p_cmb = f_simple2->GetParameter(2);
   double f_p_cmb_err = f_simple2->GetParError(2);

   FILE* of = fopen(savefile("parameters","dat").c_str(),"w");
   for (auto var : {f_L_bkg,f_L_bkg_err,f_p_bkg,f_p_bkg_err,f_L_cmb,f_L_cmb_err,f_p_cmb,f_p_cmb_err, h_bkg->Integral(), h_cmb->Integral(),f_L_sig,f_L_sig_err,f_p_sig,f_p_sig_err})
      fprintf(of,"%f ",var);
   fclose(of);

   //fit to get values we expect from MC
   auto f_combination = new TF1("f_combination",
         "[0] * (\
         (1-[4])*(3.*( [5]*2.*(1.-x**2) + [6]*(1.-2.*[3]*x+x**2) + (1.-[5]-[6])*(1.+2.*[3]*x+x**2) )/8.)\
         + [4]*(3.*( [1]*2.*(1.-x**2) + (1-[1])*[2]*(1.-2.*[3]*x+x**2) + (1-[1])*(1-[2])*(1.+2.*[3]*x+x**2) )/8.))",
         -1,1);
   f_combination->SetParameters(1,.3,.3,1,.5);
   f_combination->FixParameter(3,.16);
   f_combination->FixParameter(5,f_L_bkg);
   f_combination->FixParameter(6,f_p_bkg);
   f_combination->FixParameter(2,f_p_sig/(1-f_L_sig));
   f_combination->FixParameter(0,2.*h_cmb->Integral()/nbins);
   f_combination->FixParameter(4,h_sig->Integral()/h_cmb->Integral());
   new TCanvas();
   h_cmb->Fit("f_combination","Q");
   hstak->Draw("histe");
   f_combination->Draw("same");
   hstak->SetMinimum(0);
   gPad->SaveAs(savefile("population_fit","pdf").c_str());

   FILE* of2 = fopen(savefile("stacked_fit","dat").c_str(),"w");
   //fprintf(of2,"f_L  :%f +/- %f\nf_s  :%f +/- %f",
   //f_combination->GetParameter(1), f_combination->GetParError(1),
   //f_combination->GetParameter(4), f_combination->GetParError(4));
   fprintf(of2,"f_L  :%f +/- %f\nf_p  :%f +/- %f\nf_s  :%f +/- %f",
   f_combination->GetParameter(1), f_combination->GetParError(1),
   f_combination->GetParameter(2), f_combination->GetParError(2),
   f_combination->GetParameter(4), f_combination->GetParError(4));
   fclose(of2);
}
