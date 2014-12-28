//____________________________________________________________________________
/*!

\program gtestRewght

\brief   A simple program to illustrate how to use the GENIE event reweighting.

\syntax  gtestRewght -f filename [-n nev]

         where 
         [] is an optional argument
         -f specifies a GENIE event file (GHEP format)
         -o specifies a GENIE output filename
         -n specifies the number of events to process (default: all)

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 19, 2010

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include <string>
#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"

// number of coefficient values to allow variation of
#define MAX_COEF 3

using namespace genie;
using namespace genie::rew;
using std::string;
using std::ostringstream;

void GetCommandLineArgs (int argc, char ** argv);
int  GetNumberOfWeights (float* coefmin, float* coefmax, float* coefinc, int kmaxinc);
bool IncrementCoefficients(float* coefmin, float* coefmax, float* coefinc,
                           int kmaxinc, float* twkvals,
                           GSystSet& syst, GReWeightNuXSecCCQE* rwccqe);

string gOptInpFilename;
string gOptOutFilename;
int    gOptNEvt;
int    gOptKmaxInc = 0;
float  gOptCoeffMin[MAX_COEF] = {0.};
float  gOptCoeffMax[MAX_COEF] = {0.};
float  gOptCoeffInc[MAX_COEF] = {0.};

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  TFile file(gOptInpFilename.c_str(),"READ");
  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  LOG("test", pNOTICE) << "Input tree header: " << *thdr;
  if(!tree){
    LOG("grwght1scan", pFATAL)
      << "Can't find a GHEP tree in input file: "<< file.GetName();
    gAbortingInErr = true;
    //PrintSyntax();
    exit(1);
  }
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  LOG("test", pNOTICE) << "Will process " << nev << " events";

  //
  // Create a GReWeight object and add to it a set of 
  // weight calculators
  //

  GReWeight rw;
  rw.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );

  //
  // Create a list of systematic params (more to be found at GSyst.h)
  // set non-default values and re-configure.
  // Weight calculators included above must be able to handle the tweaked params.
  // Each tweaking dial t modifies a physics parameter p as:
  // p_{tweaked} = p_{default} ( 1 + t * dp/p )
  // So setting a tweaking dial to +/-1 modifies a physics quantity
  // by +/- 1sigma.
  // Default fractional errors are defined in GSystUncertainty
  // and can be overriden.
  //

  GSystSet & syst = rw.Systematics();
  GReWeightNuXSecCCQE * rwccqe = 
    dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeZExp);

  //
  // Concrete weight calculators can be retrieved and fine-tuned.
  // For example:

  //GReWeightNuXSecCCQE * rwccqe = 
  //  dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
  //rwccqe -> RewNue    (false); 
  //rwccqe -> RewNuebar (false); 
  //rwccqe -> RewNumubar(false); 

  // Declare the weights and twkvals arrays 
  const int n_events = (const int) nev;
  const int n_params = (const int) gOptKmaxInc;
  const int n_points = GetNumberOfWeights(gOptCoeffMin,gOptCoeffMax,gOptCoeffInc,gOptKmaxInc);
  float weights  [n_events][n_points];
  float twkvals  [n_points][n_params];
  // Initialize
  for (int i = 0; i < n_points; i++)
  {
    for (int j = 0; j < n_events; j++) weights[j][i] = 1.;
    for (int j = 0; j < n_params; j++) twkvals[i][j] = gOptCoeffMin[j];
  }
  // set first values for weighting
  for (int j = 0; j < n_params; j++)
  {
    rwccqe->SetCurrZExpIdx(j);
    rwccqe->SetSystematic(kXSecTwkDial_ZExpCCQE,gOptCoeffMin[j]);
    //syst.Set(kXSecTwkDial_ZExpCCQE,gOptCoeffMin[j]);
  }

  int twk_prm_idx = -1;
  // point loop
  for (int j = 0; j < n_points; j++) {
    twk_prm_idx++;
    rw.Reconfigure();
    // Event loop
    for(int i = 0; i < nev; i++) {
      tree->GetEntry(i);
  
      EventRecord & event = *(mcrec->event);
      LOG("test", pNOTICE) << event;
  
      double wght = rw.CalcWeight(event);
      LOG("test", pNOTICE) << "Overall weight = " << wght;

      // add to arrays
      weights[i][twk_prm_idx] = wght;
  
      mcrec->Clear();
    } // events
    IncrementCoefficients(gOptCoeffMin,gOptCoeffMax,gOptCoeffInc,
                          gOptKmaxInc,twkvals[twk_prm_idx],syst,rwccqe);
  }   // points

  // Close event file
  file.Close();

  //
  // Save weights 
  //

  // Make an output tree for saving the weights. As only considering 
  // varying a single systematic use this for name of tree.
  TFile * wght_file = new TFile(gOptOutFilename.c_str(), "RECREATE");
  TTree * wght_tree = new TTree(GSyst::AsString(kXSecTwkDial_ZExpCCQE).c_str(),
                                "GENIE weights tree");
  int branch_eventnum = 0;
  TArrayF *  branch_weight_array   = new TArrayF(n_points);
  TArrayF ** branch_twkdials_array = new TArrayF* [n_params];

  wght_tree->Branch("eventnum", &branch_eventnum);
  wght_tree->Branch("weights",  &branch_weight_array);
  ostringstream twk_dial_brnch_name;
  for (int i = 0; i < n_params; i++) {
    twk_dial_brnch_name.str("");
    twk_dial_brnch_name << "twk_dial_param_" << i+1;
    branch_twkdials_array[i] = new TArrayF(n_points);
    wght_tree->Branch(twk_dial_brnch_name.str().c_str(), branch_twkdials_array[i]);
  }

  // Compatibility with Rwght1Scan
  int nfirst=0;
  int nlast=nev;
  ostringstream str_wght;
  for(int iev = nfirst; iev <= nlast; iev++) {
    //int idx = iev - nfirst;
    int idx = iev;
    branch_eventnum = iev;
    for(int ith_comb = 0; ith_comb < n_points; ith_comb++){
       str_wght.str(", tweaked parameter values");
       for (int i = 0; i < n_params; i++) {
          if (i > 0) str_wght << ", ";
          str_wght << i+1 << " -> " << twkvals[ith_comb][i];
       }
       LOG("grwght1scan", pDEBUG)
          << "Filling tree with wght = " << weights[idx][ith_comb]
          << str_wght.str();
       branch_weight_array   -> AddAt (weights [idx][ith_comb], ith_comb);
       for (int i = 0; i < n_params; i++) {
         branch_twkdials_array[i] -> AddAt (twkvals[ith_comb][i], ith_comb);
       }
    } // twk_dial loop
    wght_tree->Fill();
  }

  wght_file->cd();
  wght_tree->Write();
  delete wght_tree;
  wght_tree = 0;
  wght_file->Close();

  LOG("test", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("test", pINFO) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {  
    LOG("testRwghtAxFF", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("testRwghtAxFF", pFATAL) 
      << "Unspecified input filename - Exiting";
    exit(1);
  }

  // output weight file
  if(parser.OptionExists('o')) {
    LOG("testRwghtAxFF", pINFO) << "Reading requested output filename";
    gOptOutFilename = parser.ArgAsString('o');
  } else {
    LOG("testRwghtAxFF", pINFO) << "Setting default output filename";
    //ostringstream nm;
    //nm << "weights_" << GSyst::AsString(gOptSyst) << ".root";
    //gOptOutFilename = nm.str();
    gOptOutFilename = "test_rw_axff_zexp.root";
  }

  // number of events:
  if( parser.OptionExists('n') ) {  
    LOG("testRwghtAxFF", pINFO) << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("testRwghtAxFF", pINFO)
       << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }

  //if( parser.OptionExists('t') ) {
  //  LOG("testRwghtAxFF", pINFO) << "Reading Number of Tweaks";
  //  string coef = parser.ArgAsString('t');

  //}

  // coefficient ranges:
  if( parser.OptionExists('c') ) {
    LOG("testRwghtAxFF", pINFO) << "Reading Coefficient Ranges";
    string coef = parser.ArgAsString('c');
    
    // split into sections of min,max,inc(rement)
    vector<string> coefrange = utils::str::Split(coef, ",");
    assert(coefrange.size() % 3 == 0);
    gOptKmaxInc = coefrange.size() / 3;
    LOG("testRwghtAxFF", pINFO) << "Number of ranges to implement : " << gOptKmaxInc;
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      gOptCoeffMin[ik] = atof(coefrange[ik*3  ].c_str());
      gOptCoeffMax[ik] = atof(coefrange[ik*3+1].c_str());
      gOptCoeffInc[ik] = atof(coefrange[ik*3+2].c_str());
    }
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      LOG("testRwghtAxFF",pWARN)<<ik+1<<": "<< gOptCoeffMin[ik]
        <<","<< gOptCoeffMax[ik]<<","<< gOptCoeffInc[ik];
    }
  } else {
    LOG("testRwghtAxFF", pFATAL) 
      << "Unspecified weighting for parameters - Exiting";
    exit(1);
  }

}
//_________________________________________________________________________________
bool IncrementCoefficients(float* coefmin, float* coefmax, float* coefinc,
                           int kmaxinc, float* twkvals,
                           GSystSet& syst, GReWeightNuXSecCCQE* rwccqe) 
{
  if (kmaxinc < 1)
  {
    LOG("gtestRwghtAxFF",pERROR) << "No coefficients to increment";
    return false;
  } else {

  int ip = -1;
  bool stopflag;
  do
  {
    if (ip > -1)
    { // a previous iteration went over max
      twkvals[ip] = coefmin[ip];
      syst.Set(kXSecTwkDial_ZExpCCQE, twkvals[ip]);
    }
    stopflag = true;

    ip++;                             // increment index
    rwccqe->SetCurrZExpIdx(ip);
    if (ip == kmaxinc) return false;  // done with incrementing

    twkvals[ip] += coefinc[ip];
    syst.Set(kXSecTwkDial_ZExpCCQE, twkvals[ip]);
    if (twkvals[ip] > coefmax[ip]) stopflag=false; // went over

  } while (! stopflag); // loop

  return true;
  } // if kmaxinc < 1

  return false;
}
//_________________________________________________________________________________
int GetNumberOfWeights(float* coefmin, float* coefmax, float* coefinc, int kmaxinc)
{
  float range = 0.;
  int  num_pts = 1;
  for (int i=0;i<kmaxinc;i++)
  {
    range = coefmax[i]-coefmin[i];
    num_pts *= (int) (range > 0 ? range/coefinc[i] : 1);
  }
  return num_pts;
}
//_________________________________________________________________________________
