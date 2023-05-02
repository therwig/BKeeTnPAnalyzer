#include "../interface/Def.h"
bool CutBasedLooseID(float full5x5_sigmaIetaIeta,float dEtaSeed,float dPhiIn, float HoverE,float relIso,float Ep,int ExpMissInnerHits,
                   bool PassConversionVeto, float E_SC, float rho, float pT, float eta_SC)
{
  if(std::abs(eta_SC)<=1.479)
  {
    if(full5x5_sigmaIetaIeta<0.0112 && std::abs(dEtaSeed)<0.00377 && std::abs(dPhiIn)<0.0884 && HoverE<0.05+1.16/E_SC+0.0324*rho/E_SC
       && std::abs(Ep)<0.193 && ExpMissInnerHits<=1 && PassConversionVeto)
       return 1;
    else
    return 0;
  }
  else
  {
    if(full5x5_sigmaIetaIeta<0.0425 && std::abs(dEtaSeed)<0.00674 && std::abs(dPhiIn)<0.169 && HoverE<0.0441+2.54/E_SC+0.183*rho/E_SC
       && std::abs(Ep)<0.111 && ExpMissInnerHits<=1 && PassConversionVeto)
       return 1;
    else
    return 0;
  }
}
bool CutBasedMediumID(float full5x5_sigmaIetaIeta,float dEtaSeed,float dPhiIn, float HoverE,float relIso,float Ep, int ExpMissInnerHits,
                   bool PassConversionVeto, float E_SC, float rho, float pT, float eta_SC)
{
  if(std::abs(eta_SC)<=1.479)
   {
     if(full5x5_sigmaIetaIeta<0.0106 && std::abs(dEtaSeed)<0.0032 && std::abs(dPhiIn)<0.0547 && HoverE<0.046+1.16/E_SC+0.0324*rho/E_SC,relIso<0.0478+0.506/pT
        && std::abs(Ep)<0.184 && ExpMissInnerHits<=1 && PassConversionVeto)
        return 1;
     else
     return 0;
   }
   else
   {
     if(full5x5_sigmaIetaIeta<0.0387 && std::abs(dEtaSeed)<0.00632 && std::abs(dPhiIn)<0.0394 && HoverE<0.0275+2.52/E_SC+0.183*rho/E_SC,relIso<0.0658+0.963/pT
        && std::abs(Ep)<0.0721 && ExpMissInnerHits<=1 && PassConversionVeto)
        return 1;
     else
     return 0;
   }
}
bool CutBasedTightID(float full5x5_sigmaIetaIeta,float dEtaSeed,float dPhiIn, float HoverE,float relIso,float Ep, int ExpMissInnerHits,
                   bool PassConversionVeto, float E_SC, float rho, float pT, float eta_SC)
{
  if(std::abs(eta_SC)<=1.479)
  {
    if(full5x5_sigmaIetaIeta<0.0104 && std::abs(dEtaSeed)<0.00255 && std::abs(dPhiIn)<0.022 && HoverE<0.026+1.15/E_SC+0.0324*rho/E_SC,relIso<0.0287+0.506/pT
       && std::abs(Ep)<0.159 && ExpMissInnerHits<=1 && PassConversionVeto)
       return 1;
    else
    return 0;
  }
  else
  {
    if(full5x5_sigmaIetaIeta<0.0353 && std::abs(dEtaSeed)<0.00501 && std::abs(dPhiIn)<0.0236 && HoverE<0.0188+2.06/E_SC+0.183*rho/E_SC,relIso<0.0445+0.963/pT
       && std::abs(Ep)<0.0197 && ExpMissInnerHits<=1 && PassConversionVeto)
       return 1;
    else
    return 0;
  }
}
