#include "selector.hh"

utils::Selector::Selector(std::string cut_type):pt_subleading_(0), pt_leading_(0), cut_type_(cut_type),
totalCount(0), isolCut(0), ptCut(0)
{}

utils::Selector::~Selector() {};

bool
utils::Selector::Pass(  Float_t pt, Float_t eta,
                        Float_t phi, Float_t charge, Float_t isol      )
{
   totalCount++;

   if ( cut_type_ == "Muon" )
   {
      //isolation criteria      
      if ( isol > 0.0027 )
      {
         isolCut++;
         return false;
      }
      //pt cut
      if ( pt < 5 )
      {
         ptCut++; 
         return false;       
      }
   }

   if ( cut_type_ == "Electron" )
   {
      //isolation criteria
      if ( isol > 0.0013 )
      {
         isolCut++;
         return false;
      }
      //pt cut
      if ( pt < 7 )
      {
         ptCut++
         return false;
      }
   }

   //assign leading/subleading
   if ( pt > pt_subleading_ )
   {
      if ( pt > pt_leading_ )
      {
         pt_subleading_ = pt_leading_;
         pt_leading_ = pt;
      }
      else pt_subleading_ = pt;
   }  

   return true;
}

bool
utils::Selector::GoodLeading()
{
   if ( cut_type_ == "Muon" )
   {
      //leading pt
      if ( pt_leading_ < 20 ) return false;
      //subleading pt
      if ( pt_subleading_ < 10 ) return false;
   }

   if ( cut_type_ == "Electron" ) 
   {
      //leading pt
      if ( pt_leading_ < 20 ) return false;
      //subleading pt
      if ( pt_subleading_ < 12 ) return false;
   }

   return true;
}
