#ifndef Selector_hh
#define Selector_hh

namespace utils
{
   class Selector
   {
      public:
         Selector(std::string cut_type);
         ~Selector();

         bool Pass(  Float_t pt, Float_t eta,
                     Float_t phi, Float_t charge, Float_t isol      );

         bool GoodLeading();

         Int_t totalCount;
         Int_t isolCut;
         Int_t ptCut;

      private:
         std::string cut_type_;
         Float_t pt_leading_;
         Float_t pt_subleading_;
   };

}

#endif
