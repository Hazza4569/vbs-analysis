#ifndef PairSet_hh
#define PairSet_hh

#include <set>

namespace utils
{
   class PairSet
   {
      public:
         PairSet(bool simple=false);
         ~PairSet();

         void AddPair(Int_t id_1, Int_t id_2, TLorentzVector combined_fourmomentum);
         void AddPair(std::pair< std::pair<Int_t,Int_t>, TLorentzVector> in_pair);

         void SetBenchmark(Double_t benchmark);

         //Get N pairs that best fit the chosen metric. If requireN is true then return
         //an empty set if a complete N pairs is not found, otherwise return all pairs found.
         PairSet GetBestNPairs(Int_t N, bool requireN = true);

         Double_t M(Int_t  n); //get the invariant mass of the nth pair

         TLorentzVector Fourmomentum(Int_t n); //get fourmomentum of nth pair

         Int_t GetNPairs();

         Double_t Metric(TLorentzVector P);

         void SetMetric(std::string metric);

         TH1F *hist_ranks;
         TH1F *hist_options;

      private:
         struct proximal {
            proximal(PairSet* ps ) { ps_ = ps; }
            bool operator() ( const std::pair< std::pair<Int_t,Int_t>, TLorentzVector >& lhs,
                              const std::pair< std::pair<Int_t,Int_t>, TLorentzVector >& rhs ) const
            {
               
               return ps_->Metric(lhs.second) < ps_->Metric(rhs.second);
            }

            PairSet* ps_;
         };

         bool CommonID( std::pair<Int_t, Int_t> ids_1, std::pair<Int_t, Int_t> ids_2 );

         std::list<std::vector<Int_t> > FindCombinations(Int_t total, Int_t n_nums); 

         void FindCombinationsUtil(  Int_t arr[], Int_t index, Int_t num, Int_t reducedNum,
               Int_t quant, std::vector<Int_t> * rank_group,
               std::list<std::vector<Int_t>>* rank_list );

         //set of pairs of ID's and their combined fourmomentum. Ordered based on
         //proximity to the benchmark.
         //Using vector as indices are helpful when finding best pairs. Will not
         //sort on entry, but before access instead.
         std::vector< std::pair< std::pair<Int_t,Int_t>, TLorentzVector> > set_;

         //benchmark: value to compare invariant masses to.
         Double_t benchmark_;
         std::string metric_;

         bool simple_; //if true then no counting histograms are created. To avoid memory issues.
   };

}

#endif
