#ifndef Object_h
#define Object_h

namespace utils {
   class Object {
      public :
         Object(std::string type, FLoat_t pt, Float_t eta, Float_t phi,
                  Float_t mass, Float_t charge);

         

         std::string type_;
         Float_t pt_;
         Float_t eta_;
         Float_t phi_;
         Float_t mass_;
         Float_t charge_;

         TLorentzVector 4_mom_;
   };

}
