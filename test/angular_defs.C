using floats = ROOT::VecOps::RVec<float>;
#include "../utils/constants.h"
void angular_defs(string file, bool verbose=false)
{
   std::string pre("/home/user108/y4p/root_output/"),post("_o.root"),newpost("_angle.root");

   ROOT::RDataFrame d("EventTree",(pre+file+post).c_str());

   string calcmethod = "default";

   auto combine_pairs = [](auto a, auto b) { return make_pair( a.first+a.second, b.first+b.second ); };
   auto pair_fourvec = [](auto a) { return (a.first+a.second); };
   auto cos_theta_star = [&verbose,&combine_pairs,&pair_fourvec,&calcmethod]
                         (floats &mu_pt, floats &mu_eta, floats &mu_phi, float mu_m, floats &mu_q,
                          floats &e_pt , floats &e_eta , floats &e_phi , float e_m , floats &e_q )
                        ->ROOT::VecOps::RVec<float>
   {
      floats rtn;
      // *********************************** LEPTON PAIR RECO *********************************** //
      // Find two Z bosons from the leptons. The pairings chosen will be those which have the closest
      // reconstructed masses to the nominal Z mass.
      
      //make loopable over lepton types: vectors of RVecs
      vector<ROOT::VecOps::RVec<float>> l_pts{mu_pt,e_pt}, l_etas{mu_eta,e_eta},l_phis{mu_phi,e_phi}, l_qs{mu_q,e_q};
      vector<float> l_ms{mu_m,e_m};
      vector<int> l_ns{(int)mu_pt.size(),(int)e_pt.size()}; //not checking for inconsistent lengths
      
      //events with < 4 leptons aren't usable
      if ( l_ns[0]+l_ns[1] < 4 ) return rtn;

      //set up ordered set
      struct pair_Zorder { bool operator() (const pair<TLorentzVector,TLorentzVector>& l, const pair<TLorentzVector,TLorentzVector>& r) const
         { double lm((l.first+l.second).M()), rm((r.first+r.second).M()); return fabs(lm-Z_MASS) < fabs(rm-Z_MASS);  }}; 

      set< pair<TLorentzVector,TLorentzVector>, pair_Zorder >  l_pairs;

      //loop over lepton type:
      for (int ilep = 0; ilep < 2; ilep++)
      {
         //form lists of each charge:       (not checking for abnormal charges here, this should be checked [poss. manually] in advance.)
         vector<TLorentzVector> l_negatives, l_positives;
         for ( int i = 0; i < l_ns[ilep]; i++ )
         {
            TLorentzVector curr_lepton;
            curr_lepton.SetPtEtaPhiM(l_pts[ilep].at(i), l_etas[ilep].at(i), l_phis[ilep].at(i), l_ms[ilep]);
            ( (l_qs[ilep].at(i)<0) ? l_negatives : l_positives ).push_back(curr_lepton);
         }
         //add pairs into ordered set:
         for (auto l_pos : l_positives) for (auto l_neg : l_negatives) l_pairs.insert(pair<TLorentzVector,TLorentzVector>(l_neg,l_pos));
      }

      if ( l_pairs.size() < 2 ) return rtn;
      //NOTE: this does not completely check that physically possible pairs have been made (i.e. no double use of leptons).
      //this should be done with filters beforehand. Checks here are just for safe memory access not good physics.

      //extract best pairs: 
      //Z1 -> l1(-), l2(+)   Z2 -> l3(-), l4(+)
      vector<pair<TLorentzVector,TLorentzVector>> pairs; 
      for (auto it=l_pairs.begin(); pairs.size() < 2; it++) pairs.push_back(*it);
      // **************************************************************************************** //

      TVector3 diboson_boost = -pair_fourvec( combine_pairs(pairs.at(0),pairs.at(1)) ).BoostVector();

      //loop over the two lepton pairs:
      for (int pair_iter = 0; pair_iter < 2; pair_iter++)
      {
         auto cur_pair = pairs.at(pair_iter);
         auto alt_pair = pairs.at(1-pair_iter);

         //find theta*:
         auto pair_boost = -pair_fourvec(cur_pair).BoostVector();
         TLorentzVector lneg(cur_pair.first), lpos(cur_pair.second),
                        Z_cur(pair_fourvec(cur_pair)), Z_alt(pair_fourvec(alt_pair));

         double costhetastar;

         //first calculation
         //Z_cur.Boost(diboson_boost);
         //costhetastar = cos( Z_cur.Angle( lneg.Vect() ) );

         if (calcmethod == "default")
         {
            lneg.Boost(pair_boost);
            lpos.Boost(pair_boost);

            Z_alt.Boost(pair_boost);
            costhetastar = cos( (-Z_alt.Vect()).Angle( lneg.Vect() ) );
         }
         else if (calcmethod == "labZ")
         {
            lneg.Boost(pair_boost);
            lpos.Boost(pair_boost);

            costhetastar = cos( (-Z_cur.Vect()).Angle( lneg.Vect() ) );
         }

         rtn.push_back(costhetastar);
      }
      return rtn;
   };

   auto d_defd = d.Define("Truth_Costhetastar",cos_theta_star,{"Truth_Muon_Pt","Truth_Muon_Eta","Truth_Muon_Phi","Truth_Muon_Mass","Truth_Muon_Charge",
                        "Truth_Electron_Pt","Truth_Electron_Eta","Truth_Electron_Phi","Truth_Electron_Mass","Truth_Electron_Charge"});

   d_defd.Snapshot("EventTree",pre+file+newpost);
}
