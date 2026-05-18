//SEIRDas

#include <Rcpp.h>
using namespace Rcpp;
//using namespace std;
#include <array>
// [[Rcpp::export]]
List model(List parscpp) { 

  //natural history parameters
  const double beta( parscpp["beta"]);     //probability of infection of susceptible upon contact x contact renormalisation 
  const double rEI(  parscpp["rEI"]);      //1/latent period - time to infection
  const double rI1I2(parscpp["rI1I2"]);    //1/time to clinical symptoms
  const double rI2R( parscpp["rI2R"]);     //1/infectious period of clinical infections
  const double rUR(  parscpp["rUR"]);      //1/infectious period of sub-clinical infections
  const double rE    = 2*rEI;  
  const double rU    = 2*rUR;  
  const double f(    parscpp["f"]);        //infectiousness of sub-clinical relative to clinical
  //vaccination parameters
  const NumericVector vcov( parscpp["vcov"]);   //vacc coverage by age and SE
  const NumericVector veff( parscpp["veff"]);   //vacc efficacy by age and SE
  const double          rV( parscpp["rV"]);     //vacc rate
  const double        vcln( parscpp["vcln"]);   //vacc reduction in clinical fraction
  
  //integration parameters
  const double dt(   parscpp["dt"]);       //time step for integration
  const int    nt(   parscpp["nt"]);       //number of time steps
  const int    nd(   parscpp["nd"]);       //number of days (start of day, to end of last day)
  const int    na(   parscpp["na"]);       //number of age groups
  const int    ns(   parscpp["nimd"]);     //number of se groups, e.g. imd
  const int    ng = na*ns;                 //number of age x se groups
  
  //it index for start of each day
  IntegerVector iW(nd);
  //vector parameters
  const std::vector<double> y45(parscpp["y45"]); //critical fraction age x SES
  const std::vector<double>  u( parscpp["u"]);   //susceptibility
  const std::vector<double> mI( parscpp["mI"]);  //rate of mortality of clinically infected
  const std::vector<double> rrep( parscpp["rrep"]); //rate of reporting
  
  // Contact matrices
  const std::vector<double> cm( parscpp["cm"]);  //ng x ng matrix - rows=participant, cols=contact
  const int   cmdim1( parscpp["cmdim1"]);  //row dimension of cm
  
  // initial states at time[0] - age-vector
  const std::vector<double> Sg0( parscpp["Sg0"]);  // Susceptible
  const std::vector<double> E1g0(parscpp["E1g0"]); // Exposed
  const std::vector<double> E2g0(parscpp["E1g0"]); // Exposed
  const std::vector<double> U1g0(parscpp["U1g0"]); // Sub-clinical cases
  const std::vector<double> U2g0(parscpp["U2g0"]); // Sub-clinical cases
  const std::vector<double> I1g0(parscpp["I1g0"]); // Pre-clinical cases
  const std::vector<double> I2g0(parscpp["I2g0"]); // Clinical cases
  const std::vector<double> Rg0(parscpp["Rg0"]);   // Recovered
  const std::vector<double> Dg0(parscpp["Dg0"]);   // Dead
  const std::vector<double> oNg(parscpp["oNg"]);   // 1/number per age group - constant over time
  
  // States(age,time) = 0 - need initialise S(na,0) etc
  std::vector<double> S_0(ng),  Sv_0(ng);  //current and next time step
  std::vector<double> E1_0(ng), Ev1_0(ng); //vector requires #include array
  std::vector<double> E2_0(ng), Ev2_0(ng);  
  std::vector<double> I1_0(ng), Iv1_0(ng);  
  std::vector<double> I2_0(ng), Iv2_0(ng);   
  std::vector<double> U1_0(ng), Uv1_0(ng);    
  std::vector<double> U2_0(ng), Uv2_0(ng);    
  std::vector<double>  R_0(ng),  Rv_0(ng);    
  std::vector<double>  D_0(ng),  Dv_0(ng), Vc_0(ng);   
  std::vector<double> Cc_0(ng), Cvc_0(ng);   
  
  // States(time) = 0 - need initialise St[0] etc, time[0]
  NumericVector time(nt);
  std::vector<double> St(nt); //TODO: use NumericVector if returning at the end
  std::vector<double> E1t(nt); 
  std::vector<double> E2t(nt); 
  std::vector<double> I1t(nt); 
  std::vector<double> I2t(nt); 
  std::vector<double> U1t(nt); 
  std::vector<double> U2t(nt); 
  std::vector<double> Rt(nt);  
  std::vector<double> Dt(nt);  
  std::vector<double> Cct(nt); //cumulative clinical cases
  
//LATER: dayly by age
  NumericVector Sw(nd); //, Svpw(nd); // Rcpp::clone(Sd);
  NumericVector Ew(nd); //, Evpw(nd); // if wanted output no. vaccinated
  NumericVector Uw(nd); //, Ivpw(nd);
  NumericVector Iw(nd); //, Uvpw(nd);
  NumericVector Rw(nd); //, Rvpw(nd);
  NumericVector Dw(nd); //,
  NumericVector Ccw(nd); //, Cvcpw(nd);
  //SE-stratified daily incidence: rows=day, cols=SE group (size nd x ns)
  NumericMatrix Ew_s(nd, ns);
  NumericMatrix Iw_s(nd, ns);
  NumericMatrix Uw_s(nd, ns);
  NumericMatrix Rw_s(nd, ns);
  //age-stratified daily incidence: rows=day, cols=age group (size nd x na)
  NumericMatrix Iw_a(nd, na);
  NumericMatrix Uw_a(nd, na);
  
  int ig;
  double Sig, E1ig, E2ig, I1ig, I2ig, U1ig, U2ig, Rig, Dig;
  
  // age group and population states initialised
  // Note: vaccinated states are empty initially (Sv_0[]=0 etc)
  for (int is = 0; is < ns; is++) { //ses
  for (int ia = 0; ia < na; ia++) { //age
    ig = ia + is*na;       
    Sig = Sg0[ig]; E1ig = E1g0[ig]; E2ig = E2g0[ig]; I1ig = I1g0[ig]; I2ig = I2g0[ig]; U1ig = U1g0[ig]; U2ig = U2g0[ig]; Rig = Rg0[ig]; Dig = Dg0[ig];
      //Assuming: each imd using same by-age-IC
      S_0[ig]  = Sig;        St[0]  += Sig;
      E1_0[ig] = E1ig;       E1t[0] += E1ig;
      E2_0[ig] = E2ig;       E2t[0] += E2ig;
      I1_0[ig] = I1ig;       I1t[0] += I1ig;
      I2_0[ig] = I2ig;       I2t[0] += I2ig;
      U1_0[ig] = U1ig;       U1t[0] += U1ig;
      U2_0[ig] = U2ig;       U2t[0] += U2ig;
      R_0[ig]  = Rig;        Rt[0]  += Rig;
      D_0[ig]  = Dig;        Dt[0]  += Dig;
      //TODO: Ew, Iw, Uw,
  }} //is, ia

  // Dynamics of state variables - fine Euler integration
  time[0]    = 0;
  iW[0]      = 0;  
  int  day  = 1;
  int  day0 = 1;
  
  double Spw = 0,   Epw = 0,   Upw = 0,   Ipw = 0,   Rpw=0,   Dpw = 0,   Ccpw = 0;
  //double Svpw = 0, Evpw = 0,  Uvpw = 0,  Ivpw = 0,  Rvpw = 0,           Cvcpw = 0; //vaccination
  
  NumericVector Epw_s(ns),   Upw_s(ns),   Ipw_s(ns),   Rpw_s(ns);
  NumericVector Upw_a(na),   Ipw_a(na);


  //temp parameters
  int    icm, ig2;                   
  double yas, ua, mIa, rrepa, cmi;
  double rVas, veas; //, Ng;  //vaccination
  //temp variables - naive and vacc populations
  double Sat,  E1at,  E2at,  U1at,  U2at,  I1at,  I2at,  Rat,  Dat;
  double Svat, Ev1at, Ev2at, Uv1at, Uv2at, Iv1at, Iv2at, Rvat; //, Dvat; //,  Nvg;  //vaccination
  //changes
  double dS,  dE1,  dE2,  dI1,  dI2,  dU1,  dU2,  dR,  dD,  dCc,  FOI,  FOIS; 
  double dSv, dEv1, dEv2, dIv1, dIv2, dUv1, dUv2, dRv, dDv, dCvc,       FOIvS, dVc; 
  //rates and inputs
  double rE1,  rE2,  rU1,  rU2,  rI1,  rI2,  dEin,  dIin,  dUin;
  double rEv1, rEv2, rUv1, rUv2, rIv1, rIv2; //, dEvin, dIvin, dUvin;    //vaccination
  
  
  for (int it = 0; it < (nt-1); it++) {	//Crucial: nt-1 ////////////////////////
    day0 = day;
    day  = 1 + (int) time[it];

    for (int is = 0; is < ns; is++) { //////////////////////////////////////////
    for (int ia = 0; ia < na; ia++) { //////////////////////////////////////////
      ig = is*na + ia; 
      // vector parameters
      yas  = y45[ig]*(1-vcln);
      ua   = u[ia];
      mIa  = mI[ia];
      rrepa= rrep[ia];
      rVas = vcov[ig]*rV;
      veas = veff[ig];
      // current matrix cells
      Sat  = S_0[ig];
      E1at = E1_0[ig];
      E2at = E2_0[ig];
      I1at = I1_0[ig];
      I2at = I2_0[ig];
      U1at = U1_0[ig];
      U2at = U2_0[ig];
      Rat  =  R_0[ig];
      Dat  =  D_0[ig];
      //Ng   = 1/oNg[ig];
      // current matrix cells - vaccinated
      Svat  = Sv_0[ig];
      Ev1at = Ev1_0[ig];
      Ev2at = Ev2_0[ig];
      Iv1at = Iv1_0[ig];
      Iv2at = Iv2_0[ig];
      Uv1at = Uv1_0[ig];
      Uv2at = Uv2_0[ig];
      Rvat  =  Rv_0[ig];
      //Dvat  =  Dv_0[ig];
      //Nvg   =  Vc_0[ig];
      // force of infection on group ia
      FOI = 0;
      //double beta_ua  = beta_infectivity*u[ia];
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
           ig2 = is2*na + ia2; 
           icm = ig2*cmdim1 + ig;
           cmi = cm[icm];        
           FOI       += beta*ua*cmi*( I1_0[ig2]  +  I2_0[ig2] +  f*U1_0[ig2] +  f*U2_0[ig2] + 
                                      Iv1_0[ig2] + Iv2_0[ig2] + f*Uv1_0[ig2] + f*Uv2_0[ig2] )*oNg[ig2];
      }} //ir, ib
      //Infection update for next timestep
             FOIS = FOI*Sat, rE1 =    (rE+rVas)*E1at,   rE2 =   (rE+rVas)*E2at, 
                             rU1 =    (rU+rVas)*U1at,   rU2 =   (rU+rVas)*U2at, 
                             rI1 = (rI1I2+rVas)*I1at,   rI2 = (rI2R+rVas)*I2at,
                             dEin = dt*FOIS, 
                             dIin = dt*yas*rE2, 
                             dUin = dt*(1-yas)*rE2;
          //vaccination
          FOIvS = FOI*Svat,  rEv1 =    rE*Ev1at, rEv2 =   rE*Ev2at, 
                             rUv1 =    rU*Uv1at, rUv2 =   rU*Uv2at, 
                             rIv1 = rI1I2*Iv1at, rIv2 = rI2R*Iv2at; //,
                             //dEvin = dt*FOIvS, 
                             //dIvin = dt*yas*rEv2, 
                             //dUvin = dt*(1-yas)*rEv2;

      dS  = dt*( -FOIS          - rVas*Sat);
      dE1 = dt*(  FOIS          - rE1);
      dE2 = dt*(  rE1           - rE2);
      dI1 = dt*(  rE2*yas       - rI1);
      dI2 = dt*(  rI1           - rI2);
      dU1 = dt*(  rE2*(1-yas)   - rU1);
      dU2 = dt*(  rU1           - rU2);
      dR  = dt*(  rI2*(1-mIa)   + rU2);
      dD  = dt*(  rI2*mIa           );
      dCc = dt*(  yas*rE2*rrepa     );
      //vaccination
      // TODO: vaccination count - rescale rVas * (Ng -Nvg)/Nvg = no.unvacc/no.vacc
      dSv  = dt*( -FOIvS               +  rVas*Sat*(1-veas));
      dEv1 = dt*(  FOIvS        - rEv1 +  rVas*E1at);
      dEv2 = dt*(  rEv1         - rEv2 +  rVas*E2at);
      dIv1 = dt*(  yas*rEv2     - rIv1 +  rVas*I1at);
      dIv2 = dt*(  rIv1         - rIv2 +  rVas*I2at);
      dUv1 = dt*(  rEv2*(1-yas) - rUv1 +  rVas*U1at);
      dUv2 = dt*(  rUv1         - rUv2 +  rVas*U2at);
      dRv  = dt*(  rIv2*(1-mIa) + rUv2 +  rVas*(Rat + veas*Sat));
      dDv  = dt*(  rIv2*mIa         );
      dVc  = dt*(  rVas*(Sat + E1at + E2at + I1at + I2at + U1at + U2at + Rat));
      dCvc = dt*(  yas*rEv2*rrepa   );      
      
      S_0[ig]  = Sat  + dS;
      E1_0[ig] = E1at + dE1;
      E2_0[ig] = E2at + dE2;
      I1_0[ig] = I1at + dI1;
      I2_0[ig] = I2at + dI2;
      U1_0[ig] = U1at + dU1;
      U2_0[ig] = U2at + dU2;
      R_0[ig]  = Rat  + dR;
      D_0[ig]  = Dat  + dD;
      Cc_0[ig] += dCc;
      // vacc update
      Sv_0[ig]  = Svat  + dSv;
      Ev1_0[ig] = Ev1at + dEv1;
      Ev2_0[ig] = Ev2at + dEv2;
      Iv1_0[ig] = Iv1at + dIv1;
      Iv2_0[ig] = Iv2at + dIv2;
      Uv1_0[ig] = Uv1at + dUv1;
      Uv2_0[ig] = Uv2at + dUv2;
      Rv_0[ig]  = Rvat  + dRv;
      Dv_0[ig]  += dDv; //Dvat  + dDv;
      Vc_0[ig]  += dVc;  //number vaccinated so far in group ig
      Cvc_0[ig] += dCvc;      
      // time counters
      time[it+1]  = (it+1)*dt;
      day        = 1 + (int) time[it+1];
      // population states evaluated (t>0)
      St[it+1]  += Sat  + dS;
      E1t[it+1] += E1at + dE1;
      E2t[it+1] += E2at + dE2;
      I1t[it+1] += I1at + dI1;
      I2t[it+1] += I2at + dI2;
      U1t[it+1] += U1at + dU1;
      U2t[it+1] += U2at + dU2;
      Rt[it+1]  += Rat  + dR;
      Dt[it+1]  += Dat  + dD;
      Cct[it+1] += Cct[it] + dCc;
      // day incidence (pw) 
      Spw       += dS;
      Epw       += dEin;
      Upw       += dUin;
      Ipw       += dIin;
      Rpw       += dR;
      Dpw       += dD;
      Ccpw      += dCc;
      //vaccination
      //Svpw      += dSv;
      //Evpw      += dEvin;
      //Uvpw      += dUvin;
      //Ivpw      += dIvin;
      //Rvpw      += dRv;
      //Cvcpw     += dCvc;      
      // day incidence - SE stratified - each is=1:5 - add age 1:9
      Epw_s[is]  +=  dEin;
      Ipw_s[is]  +=  dIin;
      Upw_s[is]  +=  dUin;
      Rpw_s[is]  +=  dR;
      // day incidence - age stratified - each ia=1:9 - add is=1:5
      Ipw_a[ia]  +=  dIin;
      Upw_a[ia]  +=  dUin;
      }; //ia ////////////////////////////////////////////////////////////////////
    }; //is ////////////////////////////////////////////////////////////////////
    
    if (day - day0 == 1) { //day incidence - update at the end of each day
      iW[day-1]  = it+1;
      Sw[day-1]  = Spw;  Spw=0;
      Ew[day-1]  = Epw;  Epw=0; 
      Iw[day-1]  = Ipw;  Ipw=0; 
      Uw[day-1]  = Upw;  Upw=0; 
      Rw[day-1]  = Rpw;  Rpw=0; 
      Dw[day-1]  = Dpw;  Dpw=0;
      Ccw[day-1] = Ccpw; Ccpw=0;
      // vaccination //if wanted output no. vaccinated
      //Svw[day-1]  = Svpw;  Svpw=0;
      //Evw[day-1]  = Evpw;  Evpw=0; 
      //Ivw[day-1]  = Ivpw;  Ivpw=0; 
      //Uvw[day-1]  = Uvpw;  Uvpw=0; 
      //Rvw[day-1]  = Rvpw;  Rvpw=0; 
      //Cvcw[day-1] = Cvcpw; Cvcpw=0;
      //SES-stratified daily incidence
      for (int is = 0; is < ns; is++) {
        Ew_s(day-1, is) = Epw_s[is]; Epw_s[is] = 0;
        Iw_s(day-1, is) = Ipw_s[is]; Ipw_s[is] = 0;
        Uw_s(day-1, is) = Upw_s[is]; Upw_s[is] = 0;
        Rw_s(day-1, is) = Rpw_s[is]; Rpw_s[is] = 0;
      }
      //age-stratified daily incidence
      for (int ia = 0; ia < na; ia++) {
        Iw_a(day-1, ia) = Ipw_a[ia]; Ipw_a[ia] = 0;
        Uw_a(day-1, ia) = Upw_a[ia]; Upw_a[ia] = 0;
      }
    }

  }; //it //////////////////////////////////////////////////////////////////////

  
  // byw: overall daily incidence vectors + SE-stratified matrices (nd x ns)
  Rcpp::List byw = Rcpp::List::create(
    Named("iW")   = iW,
    Named("time") = time[iW],
    Named("Sw")   = Sw,
    Named("Ew")   = Ew,
    Named("Iw")   = Iw,
    Named("Uw")   = Uw,
    Named("Dw")   = Dw,
    Named("Iw_s") = Iw_s,
    Named("Uw_s") = Uw_s,
    Named("Ew_s") = Ew_s,
    Named("Rw_s") = Rw_s);
  // byaw: age-stratified daily incidence matrices (nd x na)
  Rcpp::List byaw = Rcpp::List::create(
    Named("Iw_a") = Iw_a,
    Named("Uw_a") = Uw_a);

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
  
  
  
  