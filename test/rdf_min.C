using floats = ROOT::VecOps::RVec<float>;
void rdf_min(std::string file)
{
   std::string base_path("/home/user108/y4p/"), ROOT_dir("root/"), ROOT_ft(".root");
   ROOT::RDataFrame d("EventTree",base_path+ROOT_dir+file+ROOT_ft);

   auto get_pair = [](floats &charge, floats &isol, float isolreq) { return (int)min(charge[charge==-1 && isol<isolreq].size(), charge[charge==1 && isol<isolreq].size()); };

   auto d_new =  d.Define("Electron_Isol_Max","(float)0.05").Define("Muon_Isol_Max","(float)0.05")
      .Define("Electron_Pairs",get_pair,{"Electron_Charge","Electron_Isol","Electron_Isol_Max"})
      .Define("Muon_Pairs",get_pair,{"Muon_Charge","Muon_Isol","Muon_Isol_Max"})
      .Define("Lepton_Pairs","Electron_Pairs+Muon_Pairs");

   auto d_twopair = d_new.Filter("Lepton_Pairs > 1","two_pairs");

   d_twopair.Filter("Electron_Pairs == 2 && Muon_Pairs == 0","4e");
   d_twopair.Filter("Muon_Pairs == 2 && Electron_Pairs == 0","4mu");
   d_twopair.Filter("Electron_Pairs == 1 && Muon_Pairs == 1","2e2mu");

   d.Report()->Print();
}
