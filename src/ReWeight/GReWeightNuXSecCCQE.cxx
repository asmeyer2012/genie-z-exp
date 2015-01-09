//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.
 @ May 17, 2010 - CA
   Code extracted from GReWeightNuXSec and redeveloped in preparation for 
   the Summer 2010 T2K analyses.
 @ Oct 22, 2010 - CA
   Make static consts kModeMa and kModeNormAndMaShape public to aid
   external configuration.
 @ Nov 25, 2010 - CA
   Allow asymmetric errors
 @ Jan 29, 2013 - CA
   In SetSystematic() add protection against using systematic params which are
   inconsistent with the current reweighting mode [reported by Rik Gran].
 @ Dec 26, 2014 - AM
   Added Support for z-expansion axial form factor reweighting
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TFile.h>
#include <TNtupleD.h>
#include <cstdlib>
#include <sstream>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"
#include "ReWeight/GReWeightUtils.h"
#include "Registry/Registry.h"

//#define _G_REWEIGHT_CCQE_DEBUG_

using namespace genie;
using namespace genie::rew;
using std::ostringstream;

const int GReWeightNuXSecCCQE::kModeMa;
const int GReWeightNuXSecCCQE::kModeNormAndMaShape;
const int GReWeightNuXSecCCQE::kModeZExp;

//_______________________________________________________________________________________
GReWeightNuXSecCCQE::GReWeightNuXSecCCQE() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQE::~GReWeightNuXSecCCQE()
{
  delete [] fZExpTwkDial;
  delete [] fZExpDef;
  delete [] fZExpCurr;
  delete [] fZExpSigmaLo;
  delete [] fZExpSigmaHi;

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQE::IsHandled(GSyst_t syst)
{
   // read form factor model and compare to mode
   bool handle;

   switch(syst) {

     case ( kXSecTwkDial_NormCCQE    ) :
     case ( kXSecTwkDial_MaCCQEshape ) :
       if(fMode==kModeNormAndMaShape && strcmp(fFFModel.c_str(),kModelDipole) == 0)
       { 
          handle = true;  
       } else { 
          handle = false; 
       }
       break;

     case ( kXSecTwkDial_MaCCQE ) :
       if(fMode==kModeMa && strcmp(fFFModel.c_str(),kModelDipole) == 0)
       {
          handle = true;  
       } else { 
          handle = false; 
       }
       break;

     case ( kXSecTwkDial_ZExpCCQE ):
       if(fMode==kModeZExp && strcmp(fFFModel.c_str(),kModelZExp) == 0)
       {
          handle = true;
       } else {
          handle = false;
       }
       break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst))
  {
    LOG("ReW",pWARN) << "Systematic " << GSyst::AsString(syst) << " is not handled for algorithm "
      << fFFModel << " and mode " << fMode;
    return;
  }

  switch(syst) {
    case ( kXSecTwkDial_NormCCQE ) :
      fNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_MaCCQEshape ) :
    case ( kXSecTwkDial_MaCCQE ) :
      fMaTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_ZExpCCQE ) :
      //
      // GReWeight::Reconfigure only fixes one value
      // Fix others by explicitly looping over fZExpCurrIdx and
      //  setting them with (GReWeightNuXSecCCQE object pointer)->SetSystematic
      // 
      if (fZExpCurrIdx > -1 && fZExpCurrIdx < fZExpMaxCoef)
      {
        fZExpTwkDial[fZExpCurrIdx] = twk_dial;
      }
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;

  fZExpCurrIdx = 0;
  for (int i=0;i<fZExpMaxCoef;i++)
  {
    //sigmas for ZExp not reset
    fZExpTwkDial[i] = 0.;
    fZExpCurr   [i] = fZExpDef[i];
  } 

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::ResetZExpSigma(void)
{
  for (int i=0;i<fZExpMaxCoef;i++)
  {
    fZExpSigmaLo[i] = .1;
    fZExpSigmaHi[i] = .1;
  }
  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::SetCurrZExpIdx(int idx)
{
  fZExpCurrIdx = idx;
  GSystUncertainty * unc = GSystUncertainty::Instance();
  unc->SetUncertainty(kXSecTwkDial_ZExpCCQE,fZExpSigmaLo[idx],fZExpSigmaHi[idx]); 
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();
  // save current index to put it back after (if necessary)
  int curr_idx = fZExpCurrIdx; 

  if(fMode==kModeMa && strcmp(fFFModel.c_str(),kModelDipole) == 0) {   
     int    sign_matwk = utils::rew::Sign(fMaTwkDial);
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQE, sign_matwk);
     fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
     fMaCurr = TMath::Max(0., fMaCurr  );
  }
  else
  if(fMode==kModeNormAndMaShape && strcmp(fFFModel.c_str(),kModelDipole) == 0) {
     int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
     int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCQE,    sign_normtwk);
     double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQEshape, sign_mashtwk);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);
     fNormCurr = TMath::Max(0., fNormCurr);
     fMaCurr   = TMath::Max(0., fMaCurr  );
  }
  else
  if(fMode==kModeZExp && strcmp(fFFModel.c_str(),kModelZExp) == 0) {
     int    sign_twk = utils::rew::Sign(0.);
     double fracerr_zexp = 0.;
     // loop over all indices and update each
     for (int i=0;i<fZExpMaxCoef;i++)
     {
       this->SetCurrZExpIdx(i);
       fracerr_zexp = fracerr->OneSigmaErr(kXSecTwkDial_ZExpCCQE, sign_twk);
       fZExpCurr[i] = fZExpDef[i] * (1. + fZExpTwkDial[i] * fracerr_zexp);
     }
  }
  else {
    return;
  }

  Registry & r = *fXSecModelConfig;
  if (fMode==kModeMa || fMode==kModeNormAndMaShape)
  {
    r.Set(fMaPath, fMaCurr); 
  }
  else
  if (fMode==kModeZExp)
  {
    ostringstream alg_key;
    for (int i=0;i<fZExpMaxCoef;i++)
    {
      alg_key.str(""); // algorithm key for each coefficient
      alg_key << fZExpPath << "QEL-Z_A" << i+1;
      r.Set(alg_key.str(), fZExpCurr[i]);
    }
  }
  fXSecModel->Configure(r);
  this->SetCurrZExpIdx(curr_idx);

}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeight(const genie::EventRecord & event) 
{
  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  bool charm = event.Summary()->ExclTag().IsCharmEvent(); // skip CCQE charm
  if(charm) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  if(fMode==kModeMa && strcmp(fFFModel.c_str(),kModelDipole) == 0) {
     double wght = this->CalcWeightMa(event);
     return wght;
  } 
  else 
  if(fMode==kModeNormAndMaShape && strcmp(fFFModel.c_str(),kModelDipole) == 0) {
     double wght = 
         this->CalcWeightNorm    (event) *
         this->CalcWeightMaShape (event);
     return wght;
  }
  else
  if(fMode==kModeZExp && strcmp(fFFModel.c_str(),kModelZExp) == 0) {
     double wght = this->CalcWeightZExp(event);
     return wght;
  }

  return 1.;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcChisq()
{
  double chisq = 0.;
  if(fMode==kModeMa) {   
     chisq += TMath::Power(fMaTwkDial, 2.);
  }
  else
  if(fMode==kModeNormAndMaShape) { 
     chisq += TMath::Power(fNormTwkDial, 2.);
     chisq += TMath::Power(fMaTwkDial,   2.);
  }
  else
  if(fMode==kModeZExp) {
     for (int i=0;i<fZExpMaxCoef;i++)
     { 
       chisq += TMath::Power(fZExpTwkDial[i], 2.);
     }
  }
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQE::Init(void)
{
  AlgId id("genie::LwlynSmithQELCCPXSec","Default");

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
  fFFModel = fXSecModelConfig->GetAlg("FormFactorsAlg/AxialFormFactorModel").name;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  //this->SetMaPath("FormFactorsAlg/Ma");
  this->SetMaPath("FormFactorsAlg/AxialFormFactorModel/QEL-Ma");
  this->SetZExpPath("FormFactorsAlg/AxialFormFactorModel/");

  if (strcmp(fFFModel.c_str(),kModelDipole) == 0)
  {
    this->SetMode(kModeNormAndMaShape);
    fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
    fZExpMaxCoef = 0;
  } else
  if (strcmp(fFFModel.c_str(),kModelZExp) == 0)
  {
    this->SetMode(kModeZExp);
    fMaDef = 0.;
    fZExpMaxCoef = fXSecModelConfig->GetInt(fZExpPath + "QEL-Kmax");
  }

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;

  fZExpCurrIdx = 0;
  fZExpTwkDial = new double[fZExpMaxCoef];
  fZExpDef     = new double[fZExpMaxCoef];
  fZExpCurr    = new double[fZExpMaxCoef];
  fZExpSigmaLo = new double[fZExpMaxCoef];
  fZExpSigmaHi = new double[fZExpMaxCoef];

  ostringstream alg_key;
  for (int i=0;i<fZExpMaxCoef;i++)
  {
    alg_key.str("");
    alg_key << fZExpPath << "QEL-Z_A" << i+1;
    if (strcmp(fFFModel.c_str(),kModelZExp) == 0)
    { fZExpDef[i] = fXSecModelConfig->GetDouble(alg_key.str()); }
    else
    { fZExpDef[i] = 0.; }
    fZExpTwkDial[i] = 0.;
    fZExpCurr   [i] = fZExpDef[i];
    fZExpSigmaLo[i] = .1;
    fZExpSigmaHi[i] = .1;
  }

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  fTestFile = new TFile("./ccqe_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightNorm(const genie::EventRecord & /*event*/) 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightMa(const genie::EventRecord & event) 
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();
  interaction->ResetBit(kIAssumeFreeNucleon);

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightMaShape(const genie::EventRecord & event) 
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double new_integrated_xsec = fXSecModel    -> Integral(interaction);   
  assert(new_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/new_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (new) = " << new_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();
  interaction->ResetBit(kIAssumeFreeNucleon);

#ifdef _G_REWEIGHT_CCQE_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQE::CalcWeightZExp(const genie::EventRecord & event) 
{
  // very similar to CalcWeightMa
  bool tweaked = false;
  for (int i=0;i<fZExpMaxCoef;i++)
  {
    tweaked = tweaked || (TMath::Abs(fZExpTwkDial[i]) > controls::kASmallNum);
  }
  if(!tweaked) { return 1.0; }

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

  //LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
  //LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
  //LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
  //LOG("ReW", pDEBUG) << "new weight = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();
  interaction->ResetBit(kIAssumeFreeNucleon);

  return new_weight;
}
//_______________________________________________________________________________________
