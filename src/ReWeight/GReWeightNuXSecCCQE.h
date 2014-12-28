//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQE

\brief    Reweighting CCQE GENIE neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_H_

#include <map>
#include <string>

#include "ReWeight/GReWeightI.h"

using std::map;
using std::string;

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 static const char* kModelDipole = "genie::DipoleAxialFormFactorModel";
 static const char* kModelZExp   = "genie::ZExpAxialFormFactorModel";

 class GReWeightNuXSecCCQE : public GReWeightI 
 {
 public:
   static const int   kModeMa             = 0;
   static const int   kModeNormAndMaShape = 1;
   static const int   kModeZExp           = 2;

   GReWeightNuXSecCCQE();
  ~GReWeightNuXSecCCQE();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void SetMode     (int mode) { fMode       = mode; }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }
   void SetMaPath   (string p) { fMaPath     = p;    }

   void ResetZExpSigma (void);
   void SetCurrZExpIdx (int idx); // { fZExpCurrIdx = idx; }
   void SetCurrZExpSig (double siglo, double sighi)
                                 { fZExpSigmaLo[fZExpCurrIdx] = siglo;
                                   fZExpSigmaHi[fZExpCurrIdx] = sighi; }
   void SetZExpPath    (string p){ fZExpPath    = p;   }

 private:

   void   Init              (void);
   double CalcWeightNorm    (const EventRecord & event);
   double CalcWeightMaShape (const EventRecord & event);
   double CalcWeightMa      (const EventRecord & event);

   double CalcWeightZExp    (const EventRecord & event);

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked dipole model
   Registry *       fXSecModelConfig; ///< config in tweaked model
   string fFFModel;

   int    fMode;         ///< 0: Ma, 1: Norm and MaShape
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   string fMaPath;       ///< M_{A} path in config Registry
   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<

   int     fZExpCurrIdx;
   int     fZExpMaxCoef;
   string  fZExpPath;
   double* fZExpTwkDial;  
   double* fZExpDef;  
   double* fZExpCurr;  
   double* fZExpSigmaLo;
   double* fZExpSigmaHi;

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew   namespace
} // genie namespace

#endif
