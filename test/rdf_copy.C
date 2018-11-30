using floats = ROOT::VecOps::RVec<float>;
void rdf()
{
   ROOT::RDataFrame d("EventTree","root/ZZjj_10K.root",{"Muon_Pt","Electron_Pt"});
//   double cut_value = 0;
//   auto leadptcut = [&](floats  &v) { return (v.size() > 0) ? v.at(0) > cut_value : false; };
//   cut_value = 60; auto f2 = d.Filter(leadptcut,{"Electron_Pt"},"lead_electron_pt");
//   //auto f2 = f1.Filter("Electron_n > 3","lead_electron_pt");
//   
//   auto fnz = d.Filter("Muon_n > 0","muon_num");
//   auto f3 = fnz.Filter("Muon_Pt.at(0) > 150","muon_zeroth");
//   cut_value = 150; auto f1 = fnz.Filter(leadptcut,{"Muon_Pt"},"lead_muon_pt");

   //auto d_lepisol = d.Filter("Electron_Isol[true] > 0.0015","Electron Isolation").Filter("Muon_Isol[true] > 0.003","Muon Isolation");
   std::string object = "";
   auto get_pair = [](floats &charge, floats &isol, float isolreq) { return min(charge[charge==-1 && isol<isolreq].size(), charge[charge==-1 && isol<isolreq].size()); };

   auto d0 =  d.Define("Electron_Isol_Max","(float)0.0015"); auto d1= d0.Define("Muon_Isol_Max","(float)0.003");
   auto d2 = d1.Define("Muon_Pairs","int n_neg=Muon_Pt[Muon_Charge==-1].size(); return min(n_neg,Muon_n-n_neg);");
   auto d_new = d2.Define("Electron_Pairs",get_pair,{"Electron_Charge","Electron_Isol","Electron_Isol_Max"});

//   new TCanvas();
//   d_new.Histo1D({"muon","",6,0,5},"Muon_Pairs")->DrawCopy();
//   new TCanvas();
//   d_new.Histo1D({"electron","",6,0,5},"Electron_Pairs")->DrawCopy();

   d.Report()->Print();
}
