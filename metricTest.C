#include "initialTest.C"

#include <vector>

void metricTest()
{
   initialTest *reconstructor = new initialTest();
   reconstructor->signalRegionLowerBound = 180;
   reconstructor->signalRegionUpperBound = 230;
   reconstructor->drawOn = false;

   for ( std::string filename : {"ZZ_10K.root", "ZZjj_10K.root"} )
   {
      TFile * f = TFile::Open( (std::string("/home/user108/y4p/root/") + filename).c_str() );
      TTree* t = (TTree*)f->Get("EventTree");

      for ( std::string metric : { "mds", "smds" } )
      {
         reconstructor->optimisingMetric = metric;
         t->Process(reconstructor);
         printf( "File: %s, metric: %s  --  %d events found in signal region.\n",
                 filename.c_str(), metric.c_str(), reconstructor->nSignalEvents );
      }
   }

}
