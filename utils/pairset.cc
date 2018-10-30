#include "pairset.hh"

utils::PairSet::PairSet() : benchmark_(91.19), metric_("") {};

utils::PairSet::~PairSet() {};

void
utils::PairSet::AddPair(Int_t id_1, Int_t id_2, TLorentzVector combined_fourmomentum)
{ 
   std::pair<Int_t, Int_t> ids(id_1, id_2);

   set_.push_back( std::pair< std::pair<Int_t,Int_t>, TLorentzVector >( ids, combined_fourmomentum ) );
}

void
utils::PairSet::AddPair(std::pair< std::pair<Int_t,Int_t>, TLorentzVector> in_pair)
{
   set_.push_back(in_pair);
}

void
utils::PairSet::SetBenchmark(Double_t benchmark)
{ benchmark_ = benchmark; }

utils::PairSet
utils::PairSet::GetBestNPairs(Int_t N)
{
   if ( N > set_.size() || N < 1 )
   {
      return PairSet();
   }

   std::sort(set_.begin(), set_.end(), proximal(this));

   if ( N == 1 )
   {
      PairSet rtn;
      rtn.AddPair(set_[0]);
      return rtn;
   }

   Int_t sum = 0;
   for (Int_t i = 0; i < N; i++) sum += i;

   //will iteratively try the lowest ranked pairs (i.e. closest in mass to Z)
   //first. Need to test that they have unique particles, so not the same ID
   //in any of the ID pairs.
   //The sum represents the sum of ranks being tested, i.e. 0th and 1st rank, sum=1,
   //1st and 4th rank, sum=2. GetUniqueRankSets uses this to find N different rank
   //combinations. These combinations are then compared (if multiple) to find the set with the
   //lowest sum of invariant mass differences.
   for (; sum < set_.size(); sum++)
   {
      std::list< std::vector<Int_t> > rank_sets = FindCombinations( sum+N, N );

      std::vector<PairSet> uniques;
      //loop over sets of ranks
      for ( std::list< std::vector<Int_t> >::iterator it = rank_sets.begin(); it != rank_sets.end();
            it++)
      {
         bool failed = false;

         //check if ID pairs in ranks share a common value ID
         for (Int_t i = 0; i < N && !failed; i++)
         {
            for (Int_t j = i+1; j < N && !failed; j++)
            {
               failed = CommonID( set_.at(it->at(i)).first, set_.at(it->at(j)).first );
            }
         }
      
         if (!failed) 
         {
//            if ( uniques.size() == 0 ) std::cout << "Se found for event.\n";
            PairSet unique_pairs;
            for (Int_t i = 0; i < N; i++) unique_pairs.AddPair( set_.at(it->at(i)) );
            uniques.push_back( unique_pairs ); 
         }
      }

      if ( uniques.size() == 0 )
      {
         //std::cout << "Un-unique (sum " << sum << ")\n";
         continue;
      }  
      if ( uniques.size() == 1 )
      {
//         std::cout << "size1\n";
         return uniques[0];
      }
      else
      {
//         std::cout << "size2\n";
         Double_t smallest_metric = INT_MAX;
         Int_t smallest_index;
         for (Int_t i = 0; i < uniques.size(); i++)
         {
            PairSet iPairs = uniques[i];
            Double_t metric = 0;
            for (Int_t j = 0; j < N; j++)
            {
               metric += Metric(iPairs.Fourmomentum(j));
            }
            if ( metric < smallest_metric )
            {
               smallest_metric = metric;
               smallest_index = i;
            }
         }

         return uniques.at(smallest_index);
      }
   }
   return PairSet();
}

   bool
utils::PairSet::CommonID( std::pair<Int_t, Int_t> ids_1, std::pair<Int_t, Int_t> ids_2 )
{
   return ( ids_1.first == ids_2.first || ids_1.second == ids_2.second ||
         ids_1.first == ids_2.second || ids_1.second == ids_2.first );
}

//{
//   //FOLLOWING CODE DOES NOT DO AS IT IS SUPPOSED TO
//   //This is just so I can do some testing. Matty writing this method.
//   std::vector<Int_t> ranks;
//   int sum = 0;
//   for (int i = 0; i < n_nums-1; i++)
//   {
//      ranks.push_back(i);
//      sum += i;
//   }
//   ranks.push_back(abs(total-sum));
//
//   std::set< std::vector<Int_t> > rtn;
//   rtn.insert(ranks);
//
//   return rtn;   
//}

   Double_t
utils::PairSet::M(Int_t n)
{
   return set_.at(n).second.M();
}

   TLorentzVector
utils::PairSet::Fourmomentum(Int_t n)
{
   return set_.at(n).second;
}

   Int_t
utils::PairSet::GetNPairs()
{
   return set_.size();
}

   void
utils::PairSet::FindCombinationsUtil(  Int_t arr[], Int_t index, Int_t num, Int_t reducedNum,
      Int_t quant, std::vector<Int_t> * rank_group,
      std::list<std::vector<Int_t>>* rank_list ) 
{ 
   if ( reducedNum < 0 ) 
      return; 

   if ( reducedNum == 0 ) 
   { 
      if ( index == quant )
      {
         for (Int_t i = 0; i < index; i++) {
            rank_group->push_back( arr[i]-1 ); //-1 to conver to zero based index
         }
         //add finished vector to list then clears vector
         rank_list->push_back(*rank_group);
         rank_group->clear();
      }
      return; 
   } 

   Int_t prev = (index == 0)? 1 : arr[index-1]; 

   for (Int_t k = prev; k <= num ; k++) 
   { 
      for (Int_t j=0; j<index; j++)
      {
         if (arr[j] == k)
         {
            k++;
            break;
         }
      }

      arr[index] = k; 

      // call recursively with reduced number 
      FindCombinationsUtil(arr, index + 1, num, 
            reducedNum - k, quant, rank_group, rank_list); 
   } 
} 

   std::list< std::vector<Int_t> >
utils::PairSet::FindCombinations(Int_t total, Int_t n_nums)
{ 
   Int_t arr[total]; 

   //defines vector and list of vectors
   std::vector<Int_t> rank_group;
   std::list<std::vector<Int_t>> rank_list;

   //find all combinations 
   FindCombinationsUtil(arr, 0, total, total, n_nums, &rank_group, &rank_list); 

   return rank_list;
}

Double_t
utils::PairSet::Metric(TLorentzVector P)
{
   if ( metric_ == "mds" )
      return fabs(P.M() - benchmark_);
   else if ( metric_ == "smds" )
      return pow( P.M() - benchmark_, 2 );

   std::cout << "Metric \"" << metric_ << "\" is not valid.\n";
   return -1;
}

void
utils::PairSet::SetMetric(std::string metric)
{
   metric_ = metric;
}
