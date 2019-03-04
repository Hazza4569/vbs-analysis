using floats = ROOT::VecOps::RVec<float>;
void angular()
{
   //setup files:
   double target_luminosity = 35.9; //As CMS analysis
   std::string pre("/home/user108/y4p/root_output/"),post_uncut("_o.root"),post("_justewk_filtered.root");

   string strSig = "ZZjj_ATLAS_1M";
   string strBkg = "inclusive_ATLAS_5M";

   ROOT::RDataFrame d_uncut("EventTree",(pre+strSig+post_uncut).c_str());
   ROOT::RDataFrame d_sig("EventTree",(pre+strSig+post).c_str()), d_bkg("EventTree",(pre+strBkg+post).c_str());

   double scale_sig = target_luminosity *
      (*d_sig.Take<double>("Cross_Section").begin()) / 
      (*d_sig.Take<int>("Event_Count").begin()) ;

   double scale_bkg = target_luminosity *
      (*d_bkg.Take<double>("Cross_Section").begin()) / 
      (*d_bkg.Take<int>("Event_Count").begin()) ;

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
   auto cos_theta_star = [&prox_4vecs,&not_founds,&close_seconds](floats &lep_pt, floats &lep_eta, floats &lep_phi, floats &lep_m, floats &lep_q,
                                       floats &Z_pt,   floats &Z_eta,   floats &Z_phi,   floats &Z_m)   
   {
      int Z_n = Z_pt.size();
      int lep_n = lep_pt.size();

      if ( min( lep_pt[lep_q < 0].size() , lep_pt[lep_q > 0].size() ) < Z_n ) printf("===========================================size mismatch=======================================\n");

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
      TVector3 diboson_boost = -(Zs.at(0)+Zs.at(1)).Vect()*(1./(Zs.at(0)+Zs.at(1)).E());

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
            
            loops++;
            double prox = prox_4vecs(pos+neg,Z,true);
            
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
               if (!found) printf("worse than max int? %f\n",prox);
               if (prox < second_proximal) second_proximal = prox;
            }
         } 
         if (loops == 0) printf("loop nones\n");
         pos = pos_best; neg = neg_best;
         if (first_proximal > 1e-2)
         {
            not_founds++; printf("==================================NOT FOUND=========================\n");
            printf("reconstructed: %.3e %.3e %.3e %.3e\nstored       : %.3e %.3e %.3e %.3e\nfound %d\n\n",
                  (neg+pos).X(),(neg+pos).Y(),(neg+pos).Z(),(neg+pos).T(),
                  Z.X(),Z.Y(),Z.Z(),Z.T(), (int)found); 
         }
         if (second_proximal < 1e-2)
         {
            close_seconds++;
            printf("Close Second: %e\n", second_proximal);
         }
         rtn.push_back(first_proximal);
         //find theta*:
         //printf("neg: %.3e %.3e %.3e %.3e\npos: %.3e %.3e %.3e %.3e\nBoost...\n",
                //neg.X(),neg.Y(),neg.Z(),neg.T(),
                //pos.X(),pos.Y(),pos.Z(),pos.T()); 
         //auto lep_boost = -(pos+neg).BoostVector();
         //pos.Boost(lep_boost);
         //neg.Boost(lep_boost);
         //printf("neg: %.3e %.3e %.3e %.3e\npos: %.3e %.3e %.3e %.3e\n\n",
                //neg.X(),neg.Y(),neg.Z(),neg.T(),
                //pos.X(),pos.Y(),pos.Z(),pos.T()); 
         //Z.Boost(diboson_boost);

         ////double posangle(Z.Angle(pos.Vect())), negangle(Z.Angle(pos.Vect()));       
         //rtn.push_back(cos(min(posangle,negangle)));
      }
      return rtn;
   };

   int nbins(11);
   double xbins[nbins];
   for (int i=0; i<=nbins; i++) xbins[i] = pow(10,-nbins+i);

   TH1F* h = (TH1F*) d_uncut
      .Filter("Lepton_Pairs > 1","two Z")
      .Define("anglediff",cos_theta_star,{"Isolated_Lepton_Pt","Isolated_Lepton_Eta","Isolated_Lepton_Phi","Isolated_Lepton_M","Isolated_Lepton_Charge","Dilepton_Pt","Dilepton_Eta","Dilepton_Phi","Dilepton_M"})
      //N.B: Without making the above cut, one might assume that requiring 2 Zs (as is done in cos_theta_star) would do the job anyway.
      //This is not the case because the Dilepton properties are not requiring two different leptons, so 3 lepton events could have two Dilepton entries.

      //.Histo1D({"",";cos#theta*;Freq",200,-1,1},"anglediff")->Clone("h"); 
      .Histo1D({"",";closest Z;",nbins,xbins},"anglediff")->Clone("h");

   //TF1* f = new TF1("f","[0]*3*( [1]*2*(1-x**2) + [2]*(1+x)**2 + [3]*(1-x)**2 )/8.",-1,1);
   //f->SetParameters(1000,0.33,0.33,0.33);
   
   new TCanvas();
   gPad->SetLogx();
   h->Draw();
   //h->Fit("f");
   printf("Not found: %d\n",not_founds);
   printf("Close Seconds: %d\n",close_seconds);
   d_uncut.Report()->Print();
}
