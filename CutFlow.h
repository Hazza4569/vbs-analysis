#ifndef CutFlow_h
#define CutFlow_h

//Object should take signal and background sample, have one or more cuts applied, and show results
//of these cuts.

class CutFlow
{
   public:
      CutFlow(std::string signalFile, std::string bgFile);
      virtual ~CutFlow();

      void PreSelection( Int_t n_on_shell ); //makes 4lepton and other basic cuts to clean up data

      void JetEta(bool tag=false);
      void JetEtaDiff();
      void JetPt(Int_t n); //n=0 for all, 1 for leading, 2 for subleading etc.
      void DijetM();
      void ElectronPt(Int_t n); //same as JetPt with n
      void MuonPt(Int_t n); //same as JetPt with n

      std::vector<TLorentzVector> ReconstructVectors(std::string sample, std::string obj_type);
      bool JetIsElectron( std::string sample, Int_t jet_n );

      void ProcessAll( Int_t n_on_shell ); //par -1:no preselection 0,1,2:0,1,2 Zs required on shell
      void FileInit();

      Double_t target_luminosity = 150; //fb^{-1}

   private:
      struct pt_ordered
      {
         bool operator() ( const TLorentzVector &lhs, const TLorentzVector &rhs ) const
         {
            return lhs.Pt() > rhs.Pt();
         }
      };
      std::map< std::string, std::string > files_;
      std::map< std::string, TTree* > trees_;
      std::map< std::string, TTreeReader*> readers_;
      std::map< std::string, Int_t> num_events_;  
      std::map< std::string, Double_t> cross_section_; //fb
   
      std::map< std::string, TTreeReaderArray<Float_t>* > event_vars_;
      std::map< std::string, TTreeReaderValue<Int_t>* > object_numbers_;
      const std::string delim_ = "_";
      std::map<std::string, std::fstream*> fileEventNums_;
      std::map<std::string, std::string> eventFileName_;
      bool filedEvents_ = false;
      //std::map< std::string, std::vector<Float_t> > cut_vars;
};

#endif
