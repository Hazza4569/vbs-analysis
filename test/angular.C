using floats = ROOT::VecOps::RVec<float>;
void angular(string date, bool verbose=false)
{
   //setup files:
   double target_luminosity = 3000; //Hi-Lumi
   double mjj = 1.037639e+03;
   double njj = 2.556276e+00; //as per optimisation

   std::string pre("/home/user108/y4p/root_output/"),post_uncut("_o.root"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_1M";
   string strBkg = "inclusive_ATLAS_5M";

   ROOT::RDataFrame d_uncut("EventTree",(pre+strSig+post_uncut).c_str());
   ROOT::RDataFrame d_bguncut("EventTree",(pre+strBkg+post_uncut).c_str());
   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());

   double scale_sig = target_luminosity *
      (*d_uncut.Take<double>("Cross_Section").begin()) / 
      d_uncut.Sum("Event_Weight").GetValue() ;

   double scale_bkg = target_luminosity *
      (*d_bguncut.Take<double>("Cross_Section").begin()) / 
      d_bguncut.Sum("Event_Weight").GetValue() ;

   std::string cut = string("Dijet_M > ")+to_string(mjj)+string(" && Jet12_Eta_Diff > ")+to_string(njj);
   auto d_sigcut = d_sig.Filter(cut);
   auto d_bkgcut = d_bkg.Filter(cut);

   auto savefile = [&](string var,string ft){
      char rtn[800];
      snprintf(rtn,800,"/home/user108/y4p/graph_logs/%s/angular_%s.%s",date.c_str(),var.c_str(),ft.c_str());
      return string(rtn);
   };

   auto prox_4vecs = [](TLorentzVector a, TLorentzVector b, bool relative=true)
   {
      double X = fabs( (a.X() - b.X())/( relative? b.X() : 1 ) );
      double Y = fabs( (a.Y() - b.Y())/( relative? b.Y() : 1 ) );
      double Z = fabs( (a.Z() - b.Z())/( relative? b.Z() : 1 ) );
      double T = fabs( (a.T() - b.T())/( relative? b.T() : 1 ) );
      return max(max(max(X,Y),Z),T);
   };

   auto comp_4vecs = [&prox_4vecs](TLorentzVector a, TLorentzVector b, double relerr=1e-7)
   {
      return prox_4vecs(a,b) < relerr;
   };
   int not_founds = 0;
   int close_seconds = 0;
   auto cos_theta_star = [&prox_4vecs,&not_founds,&close_seconds,&verbose](floats &lep_pt, floats &lep_eta, floats &lep_phi, floats &lep_m, floats &lep_q,
                                       floats &Z_pt,   floats &Z_eta,   floats &Z_phi,   floats &Z_m)   
   {
      int Z_n = Z_pt.size();
      int lep_n = lep_pt.size();

      //if ( min( lep_pt[lep_q < 0].size() , lep_pt[lep_q > 0].size() ) < Z_n ) printf("===========================================size mismatch=======================================\n");

      //sort leptons into two sets of indices
//      std::set<int> pos_indi,neg_indi;
//      for (int i = 0; i < lep_n; i++)
//      {
//         if (lep_q.at(i) > 0) pos_indi.insert(i);
//         else if (lep_q.at(i) < 0) neg_indi.insert(i);
//         else printf("Error in cos_theta_star: seemingly uncharged lepton.\n");
//      }

      floats rtn;
      if (Z_n < 2) return rtn;   

      std::vector<TLorentzVector> Zs;
      for (int i = 0; i < 2; i++) { TLorentzVector Z; Z.SetPtEtaPhiM(Z_pt.at(i),Z_eta.at(i),Z_phi.at(i),Z_m.at(i));
         Zs.push_back(Z); 
      }
      TVector3 diboson_boost = -(Zs.at(0)+Zs.at(1)).BoostVector();
      if ( fabs(Zs.at(1).M() - 91.19) < fabs(Zs.at(0).M() - 91.19) ) printf("***************************************Z ORDERING?????*********************************\n");

      //main loop
      for (int Z_i = 0; rtn.size() < Z_n; Z_i++)
      {
         TLorentzVector Z; Z.SetPtEtaPhiM(Z_pt.at(Z_i),Z_eta.at(Z_i),Z_phi.at(Z_i),Z_m.at(Z_i));
         TLorentzVector pos, neg, pos_best, neg_best;
         //find which leptons make up this Z
         bool found = false;
         double first_proximal(INT_MAX), second_proximal(INT_MAX);
         int loops(0);
         for (int pos_i=0; pos_i < lep_n; pos_i++) for (int neg_i=0; neg_i < lep_n; neg_i++)
         {
            //if ( pos_i == neg_i) continue;
            if ( abs((lep_m.at(pos_i) - lep_m.at(neg_i))/lep_m.at(pos_i)) > 1e-2 ) continue;
            if ( lep_q.at(pos_i) <= 0 || lep_q.at(neg_i) >= 0 ) continue;
            //if ( lep_q.at(pos_i)*lep_q.at(neg_i) >= 0 ) {printf("================triggered===================================================\n");continue;}
            pos.SetPtEtaPhiM(lep_pt.at(pos_i),lep_eta.at(pos_i),lep_phi.at(pos_i),lep_m.at(pos_i));
            neg.SetPtEtaPhiM(lep_pt.at(neg_i),lep_eta.at(neg_i),lep_phi.at(neg_i),lep_m.at(neg_i));
            
            loops++; double prox = prox_4vecs(pos+neg,Z,true);
            
            if ( prox < first_proximal ) 
            {
               second_proximal = first_proximal;
               first_proximal = prox;
               pos_best = pos; neg_best = neg;
               found = true;
               //break;
            }
            else
            {
               if (!found && verbose) printf("worse than max int? %f\n",prox);
               if (prox < second_proximal) second_proximal = prox;
            }
         } 
         if (loops == 0 && verbose) printf("loop nones\n");
         pos = pos_best; neg = neg_best;
         if (first_proximal > 1e-2)
         {
            not_founds++; if (verbose) printf("==================================NOT FOUND=========================\n");
            if (verbose) printf("reconstructed: %.3e %.3e %.3e %.3e\nstored       : %.3e %.3e %.3e %.3e\nfound %d\n\n",
                  (neg+pos).X(),(neg+pos).Y(),(neg+pos).Z(),(neg+pos).T(),
                  Z.X(),Z.Y(),Z.Z(),Z.T(), (int)found); 
         }
         if (second_proximal < 1e-2)
         {
            close_seconds++;
            if (verbose) printf("Close Second: %e\n", second_proximal);
         }
         //rtn.push_back(first_proximal);

         //find theta*:
         if (verbose) printf("neg: %.3e %.3e %.3e %.3e\npos: %.3e %.3e %.3e %.3e\nZ  : %.3e %.3e %.3e %.3e\nZ2 : %.3e %.3e %.3e %.3e\nBoost...\n",
                neg.X(),neg.Y(),neg.Z(),neg.T(),
                pos.X(),pos.Y(),pos.Z(),pos.T(),
                Z.X(),Z.Y(),Z.Z(),Z.T(),
                Zs.at(1-Z_i).X(),Zs.at(1-Z_i).Y(),Zs.at(1-Z_i).Z(),Zs.at(1-Z_i).T()); 
         auto lep_boost = -(Z).BoostVector();
         pos.Boost(lep_boost);
         neg.Boost(lep_boost);
         Z.Boost(diboson_boost);
         if (verbose) printf("neg: %.3e %.3e %.3e %.3e\npos: %.3e %.3e %.3e %.3e\nZ  : %.3e %.3e %.3e %.3e\nZ2 : %.3e %.3e %.3e %.3e\n\n",
                neg.X(),neg.Y(),neg.Z(),neg.T(),
                pos.X(),pos.Y(),pos.Z(),pos.T(),
                Z.X(),Z.Y(),Z.Z(),Z.T(),
                Zs.at(1-Z_i).X(),Zs.at(1-Z_i).Y(),Zs.at(1-Z_i).Z(),Zs.at(1-Z_i).T()); 

         double posangle(Z.Angle(pos.Vect())), negangle(Z.Angle(neg.Vect()));       
         if (verbose) printf("%f %f %f\n",posangle,negangle,posangle+negangle);
         rtn.push_back(cos(negangle));
         //rtn.push_back(cos(min(posangle,negangle)));
      }
      return rtn;
   };

   auto definitions = [&cos_theta_star](auto df)
   {
      return df
         .Define("anglediff",cos_theta_star,{"Truth_Isolated_Lepton_Pt","Truth_Isolated_Lepton_Eta","Truth_Isolated_Lepton_Phi","Truth_Isolated_Lepton_M","Truth_Isolated_Lepton_Charge",
                                             "Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
         .Define("anglediff1","anglediff.at(0)")
         .Define("anglediff2","anglediff.at(1)");
   };

   auto d_uncutd = definitions(d_uncut).Filter("true");
   auto d_sigcutd = definitions(d_sigcut).Filter("true");
   auto d_bkgcutd = definitions(d_bkgcut).Filter("true");

   //int nbins(11);
   //double xbins[nbins];
   //for (int i=0; i<=nbins; i++) xbins[i] = pow(10,-nbins+i);
   
   auto plot_angular = [&cos_theta_star,&savefile](auto df, double scale, string title="", int order=0)
   {
      auto d_new = df.Filter("Event_Weight < 2e6 && Event_Weight > -1e4","Remove anomalous weightings");

      //not sure why this is necessary, but the dataframe seemed to be remembering the definitions between function calls
      //and throwing complaints.
      //vector<string> cols = d_old.GetColumnNames();
      //printf("%s\n",title.c_str());
      //for (auto s : cols) printf("%s\n",s.c_str());
      //bool alreadydefd = ( std::find(cols.begin(), cols.end(), "anglediff") != cols.end() );
      //printf("%d\n", (int) alreadydefd);
   
      //auto d_new = (alreadydefd) ? d_old.Filter("true") : d_old
      //   .Filter("true");

      int nbins(30);
      TH1F* h1 = (TH1F*) d_new.Histo1D({"",(title+string(";cos#theta*_{1};Events / 0.08")).c_str(),nbins,-1,1},"anglediff1","Event_Weight")->Clone(""); 
      h1->Scale(scale);
      TH1F* h2 = (TH1F*) d_new.Histo1D({"",(title+string(";cos#theta*_{2};Events / 0.08")).c_str(),nbins,-1,1},"anglediff2","Event_Weight")->Clone(""); 
      h2->Scale(scale);

//      TF1* f = new TF1("f","[0]*3*( [1]*2*(1-x**2) + [2]*(1+x)**2 + (1-[1]-[2])*(1-x)**2 )/8.",-1,1);
//      f->SetParameters(1,0.33,0.33);
//      f->SetParNames("Scale","f_{L}","f_{+}");
//      f->SetParLimits(1,0,1);
//      f->SetParLimits(2,0,1);
      TF1* f = new TF1("f","[0]*3.*( [1]*2.*(1.-x**2) + [2]*(1.-2.*[3]*x+x**2) + (1.-[1]-[2])*(1.+2.*[3]*x+x**2) )/8.",-1,1);
      f->SetParameters(1,0.33,0.33,1);
      f->SetParNames("Scale","f_{L}","f_{+}","A");
      f->SetParLimits(1,0,1);
      f->SetParLimits(2,0,1);
      f->FixParameter(3,0.16/*8./(1-4* 0.23)*/);

      gStyle->SetOptStat(1111110); gStyle->SetOptFit(1);
      gStyle->SetStatX(0.7); gStyle->SetStatY(0.9);

      h1->Draw();
      f->FixParameter(0,2.*h1->Integral()/nbins); //integral over bin width.
      h1->Fit("f");
      gPad->SaveAs(savefile(to_string(order)+title+string("_Z1"),"pdf").c_str());

      h2->Draw();
      f->FixParameter(0,2.*h2->Integral()/nbins); //integral over bin width.
      h2->Fit("f");
      gPad->SaveAs(savefile(to_string(order)+title+string("_Z2"),"pdf").c_str());
   }; 
   
   new TCanvas();
   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1").Filter("bool rtn=true; for (auto &id : Truth_Lepton_Parent_ID) if (abs(id) == 15) rtn=false; return rtn;")
                         .Filter("fabs(Truth_Dilepton_M.at(0) - 91.19) < 3") , scale_sig, "truth_Z1tight" ,0);
//   plot_angular(d_uncutd .Filter("Truth_Lepton_Pairs > 1","two truth Z") , scale_sig, "truth_preselection" ,0);
//   plot_angular(d_sigcutd.Filter("Truth_Lepton_Pairs > 1","two truth Z") , scale_sig, "truth_postselection",2);
//   plot_angular(d_uncutd .Filter("Lepton_Pairs > 1","two Z") , scale_sig, "signal_preselection" ,1);
//   plot_angular(d_sigcutd.Filter("Lepton_Pairs > 1","two Z") , scale_sig, "signal_postselection",3);
//
//   plot_angular(d_uncutd.Filter("Truth_Lepton_Pairs > 1","two Z").Filter(cut) , scale_sig, "truth_vbsenhancementonly",4);
//   plot_angular(d_uncutd.Filter("Truth_Lepton_Pairs > 1","two Z").Filter("fabs(Truth_Dilepton_M.at(0) - 91.2) < 10") , scale_sig, "truth_Z1masscut",5);
//   plot_angular(d_uncutd.Filter("Truth_Lepton_Pairs > 1","two Z").Filter("fabs(Truth_Dilepton_M.at(1) - 91.2) < 20") , scale_sig, "truth_Z2masscut",6);
//   plot_angular(d_uncutd.Filter("Truth_Lepton_Pairs > 1","two Z").Filter("fabs(Truth_Dilepton_M.at(0) - 91.2) < 10").Filter(cut) , scale_sig, "truth_Z1massandVBS",7);
   //printf("Not found: %d\n",not_founds);
   //printf("Close Seconds: %d\n",close_seconds);
   //d_uncut.Report()->Print();
}

   //TRUTH:
//   TH1F* h = (TH1F*) d_uncut
//      .Filter("Truth_Lepton_Pairs > 1","two Z")
//      .Filter("Event_Weight < 2e6 && Event_Weight > -1e4","Remove anomalous weightings")
//      //.Filter("fabs(Truth_Dilepton_Eta.at(0)) < .1 && fabs(Truth_Dilepton_Eta.at(1)) < .1","Low eta Zs")
//      //.Filter("bool rtn=true; for (auto &n : Truth_Lepton_Eta) if(n < 1.5) rtn=false; return rtn;","Low eta leptons")
//      //.Define("Z_Weight","ROOT::VecOps::RVec<float>(2,Event_Weight)")
//      .Define("anglediff",cos_theta_star,{"Truth_Isolated_Lepton_Pt","Truth_Isolated_Lepton_Eta","Truth_Isolated_Lepton_Phi","Truth_Isolated_Lepton_M","Truth_Isolated_Lepton_Charge",
//                                          "Truth_Dilepton_Pt","Truth_Dilepton_Eta","Truth_Dilepton_Phi","Truth_Dilepton_M"})
//      .Define("anglediff1","anglediff.at(0)")
//      .Define("anglediff2","anglediff.at(1)")
//      //N.B: Without making the above cut, one might assume that requiring 2 Zs (as is done in cos_theta_star) would do the job anyway.
//      //This is not the case because the Dilepton properties are not requiring two different leptons, so 3 lepton events could have two Dilepton entries.
//
//      .Histo1D({"",";cos#theta*;Events / 0.0025",nbins,-1,1},"anglediff1","Event_Weight")->Clone("Truth"); 
//      //.Histo1D({"",";closest Z;",nbins,xbins},"anglediff")->Clone("h");
   
//   h->Scale(scale_uncut);
//
//   TF1* f = new TF1("f","[0]*3*( [1]*2*(1-x**2) + [2]*(1+x)**2 + (1-[1]-[2])*(1-x)**2 )/8.",-1,1);
//   f->SetParameters(1,0.33,0.33);
//   f->SetParNames("Scale","f_{L}","f_{+}");
//   f->SetParameter(0,2.*h->Integral()/nbins); //integral over bin width.
//   f->SetParLimits(1,0,1);
//   f->SetParLimits(2,0,1);

//   TF1* f_translated = new TF1("f_translated","[0]*3*( [1]*2*(1-(x*cos([3])-sqrt(1-x**2)*sin([3]))**2) + [2]*(1+(x*cos([3])-sqrt(1-x**2)*sin([3])))**2 + (1-[1]-[2])*(1-(x*cos([3])-sqrt(1-x**2)*sin([3])))**2 )/8.",-1,1);
//   f_translated->SetParameters(1,0.33,0.33,0);
//   f_translated->SetParNames("Scale","f_{L}","f_{+}","#delta");
//   f_translated->SetParameter(0,2.*h->Integral()/nbins); //integral over bin width.
//   f_translated->SetParLimits(1,0,1);
//   f_translated->SetParLimits(2,0,1);
   
//   gStyle->SetOptStat(1111111); gStyle->SetOptFit(1);
//   new TCanvas();
//   //gPad->SetLogx();
//   h->Draw();
//   h->Fit("f");
//
//   TH1F *hres = new TH1F("Truth_residuals",";cos#theta*;Residuals / 0.0025",nbins,-1,1);
//
//   for (int i=0; i <= nbins; i++)
//   {
//      hres->SetBinContent(i,h->GetBinContent(i) - f->Eval(h->GetBinCenter(i)));
//      hres->SetBinError(i,h->GetBinError(i));
//   }
//
//   new TCanvas();
//   hres->Draw("E");
