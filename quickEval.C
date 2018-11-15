#include "initialTest.C"

#include <vector>

void quickEval(std::string filename="ZZjj_10K", std::string metric="mds")
{
   initialTest *reconstructor = new initialTest();

   TFile * f = TFile::Open( ( std::string("/home/user108/y4p/root/") + filename
                              + std::string(".root")   ).c_str() );
   TTree* t = (TTree*)f->Get("EventTree");

   reconstructor->optimisingMetric = metric;
   t->Process(reconstructor);
}
