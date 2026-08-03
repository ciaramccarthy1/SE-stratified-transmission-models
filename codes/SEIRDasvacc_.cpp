//SEIRDHV - weekly incidence
// See SEIRDasvaccday_.cpp for compartment / equation documentation.

#include <Rcpp.h>
using namespace Rcpp;
#include <array>
#include <cmath>

// [[Rcpp::export]]
List model(List parscpp) {

  const double beta0( parscpp["beta"]);       //baseline transmission (seasonal multiplier applied per day)
  const double b1(   parscpp["b1"]);          //seasonal amplitude (Gaussian pulse, David rsvie)
  const double phi(  parscpp["phi"]);         //seasonal phase (fraction of year at peak)
  const double psi(  parscpp["psi"]);         //seasonal width (fraction of year)
  const double rEI(  parscpp["rEI"]);              //latency E -> I/U (single stage, no Erlang)
  const std::vector<double> rIR( parscpp["rIR"]);  //recovery I -> R/H rate, by age (David exposure-group avg)
  const std::vector<double> rUR( parscpp["rUR"]);  //recovery U -> R   rate, by age
  const double f(    parscpp["f"]);
  const double rH(   parscpp["rH"]);

  const NumericVector vcov(   parscpp["vcov"]);
  const NumericVector VE_inf( parscpp["VE_inf"]);
  const NumericVector VE_sym( parscpp["VE_sym"]);
  const NumericVector VE_hosp(parscpp["VE_hosp"]);
  const NumericVector VE_mort( parscpp["VE_mort"]);
  const double rV(     parscpp["rV"]);
  const double rW(     parscpp["rW"]);
  const double rW_nat( parscpp["rW_nat"]);

  const double dt(   parscpp["dt"]);
  const int    nt(   parscpp["nt"]);
  const int    nw(   parscpp["nw"]);
  const int    na(   parscpp["na"]);
  const int    ns(   parscpp["nimd"]);
  const int    ng = na*ns;

  IntegerVector iW(nw);

  const std::vector<double> y45(parscpp["y45"]);
  const std::vector<double>  u( parscpp["u"]);
  const std::vector<double>  h( parscpp["h"]);
  const std::vector<double> mH( parscpp["mH"]);
  const std::vector<double> rrep(parscpp["rrep"]);

  const std::vector<double> cm( parscpp["cm"]);
  const int   cmdim1( parscpp["cmdim1"]);

  const std::vector<double> Sg0( parscpp["Sg0"]);
  const std::vector<double> Eg0( parscpp["Eg0"]);
  const std::vector<double> Ug0( parscpp["Ug0"]);
  const std::vector<double> Ig0( parscpp["Ig0"]);
  const std::vector<double> Rg0( parscpp["Rg0"]);
  const std::vector<double> Dg0( parscpp["Dg0"]);
  const std::vector<double> oNg( parscpp["oNg"]);

  std::vector<double> S_0(ng), E_0(ng);
  std::vector<double> I_0(ng), H_0(ng);
  std::vector<double> U_0(ng);
  std::vector<double> R_0(ng),  D_0(ng);
  std::vector<double> V_0(ng);
  std::vector<double> Ev_0(ng);
  std::vector<double> Iv_0(ng), Hv_0(ng);
  std::vector<double> Uv_0(ng);
  std::vector<double> Rv_0(ng),  Dv_0(ng);
  std::vector<double> Vc_0(ng), Cc_0(ng), Ccv_0(ng);

  NumericVector Sw(nw), Ew(nw), Iw(nw), Uw(nw), Hw(nw), Rw(nw), Dw(nw);
  NumericVector Ccw(nw);
  NumericVector Ewv(nw), Iwv(nw), Uwv(nw), Hwv(nw), Dwv(nw), Ccwv(nw);
  NumericVector Vw(nw);
  NumericVector Vlev(nw);
  NumericMatrix Ew_s(nw, ns), Iw_s(nw, ns), Uw_s(nw, ns), Hw_s(nw, ns), Rw_s(nw, ns);
  //(age x IMD)-stratified hospitalisation incidence (rows=week, cols=ig=is*na+ia)
  NumericMatrix Hw_ag(nw, ng), Hwv_ag(nw, ng);
  NumericMatrix Iw_a(nw, na), Uw_a(nw, na), Hw_a(nw, na);
  NumericMatrix Vlev_a(nw, na), Vlev_s(nw, ns);

  NumericVector time(nt);
  int ig;

  for (int is = 0; is < ns; is++) {
  for (int ia = 0; ia < na; ia++) {
    ig = ia + is*na;
    S_0[ig]  = Sg0[ig];
    E_0[ig]  = Eg0[ig];
    I_0[ig]  = Ig0[ig];
    U_0[ig]  = Ug0[ig];
    H_0[ig]  = 0.0;
    R_0[ig]  = Rg0[ig];
    D_0[ig]  = Dg0[ig];
    V_0[ig] = 0.0;
    Ev_0[ig] = 0.0;
    Iv_0[ig] = 0.0;
    Uv_0[ig] = 0.0;
    Hv_0[ig]  = 0.0;
    Rv_0[ig]  = 0.0; Dv_0[ig] = 0.0;
    Vc_0[ig]  = 0.0; Cc_0[ig] = 0.0; Ccv_0[ig] = 0.0;
  }}

  time[0]    = 0;
  iW[0]      = 0;
  int week   = 1;
  int week0  = 1;

  double Spw=0, Epw=0, Ipw=0, Upw=0, Hpw=0, Rpw=0, Dpw=0;
  double Ccpw=0;
  double Epw_v=0, Ipw_v=0, Upw_v=0, Hpw_v=0, Dpw_v=0, Ccpw_v=0;
  double Vpw=0;
  NumericVector Epw_s(ns), Ipw_s(ns), Upw_s(ns), Hpw_s(ns), Rpw_s(ns);
  NumericVector Ipw_a(na), Upw_a(na), Hpw_a(na);
  NumericVector Hpw_ag(ng), Hpw_vag(ng);  //weekly H accumulators per (age x IMD)

  double yas, ua, ha, mHa, rrepa, cmi;
  double rVas, ve_i, ve_y, ve_h, ve_m;
  double Sat, Eat, Uat, Iat, Hat, Rat, Dat;
  double Vat, Evat, Uvat, Ivat, Hvat, Rvat, Dvat;
  double FOI, FOIS, FOIvV;
  double rEo, rUo, rIo, rHo;
  double rEvo, rUvo, rIvo, rHvo;
  double dVin, dEin, dIin, dUin, dHin;
  double dEvin, dIvin, dUvin, dHvin;
  double dS, dE, dI, dH, dU, dR, dD, dCc;
  double dV, dEv, dIv, dHv, dUv, dRv, dDv, dCcv, dVc;
  int    icm, ig2;

  for (int it = 0; it < (nt-1); it++) {
    week0 = week;
    week  = 1 + (int) time[it]/7;

    //seasonal forcing (Gaussian pulse). Reproduces rsvie EXACTLY (Hodgson,
    //RunInterventions.h:968). Two properties are INTENTIONAL (not bugs) - keep
    //them to preserve the David calibration:
    // (1) inner (1.0+...) gives a multiplier floor of 1+b1 (~3), so beta0 is NOT
    //     a baseline-R0 beta (see R0_.r warning). Dropping it changes the seasonal
    //     amplitude and de-calibrates unless beta0/b1 are re-derived.
    // (2) phase (t1/365-phi) is NOT wrapped, so beta is discontinuous at each
    //     365-day boundary - but that boundary is the seasonal trough, and rsvie
    //     does the same (no wrap). Do not "fix" without re-fitting to match David.
    double t1   = std::fmod(time[it], 365.0);
    double beta = beta0*(1.0 + b1*(1.0 + std::exp(-(t1/365.0-phi)*(t1/365.0-phi)/(2.0*psi*psi))));

    // Snapshot the infectious compartments at time t so every group's FOI is
    // computed from the SAME (time-t) state - a Jacobi forward-Euler step. The
    // per-group loop below writes I_0/U_0/Iv_0/Uv_0 in place, so without this
    // snapshot the FOI would mix updated (ig2<ig) and un-updated (ig2>=ig)
    // values, making the result depend on is/ia iteration order.
    std::vector<double> I_s = I_0, U_s = U_0, Iv_s = Iv_0, Uv_s = Uv_0;

    for (int is = 0; is < ns; is++) {
    for (int ia = 0; ia < na; ia++) {
      ig    = is*na + ia;
      yas   = y45[ig];
      ua    = u[ia];
      ha    = h[ia];
      mHa   = mH[ia];
      rrepa = rrep[ia];
      rVas  = vcov[ig]*rV;
      ve_i  = VE_inf[ig];
      ve_y  = VE_sym[ig];
      ve_h  = VE_hosp[ig];
      ve_m  = VE_mort[ig];

      Sat  =  S_0[ig];
      Eat  =  E_0[ig];
      Iat  =  I_0[ig];
      Uat  =  U_0[ig];
      Hat  =  H_0[ig];
      Rat  =  R_0[ig]; Dat  =  D_0[ig];

      Vat   =  V_0[ig];
      Evat  =  Ev_0[ig];
      Ivat  =  Iv_0[ig];
      Uvat  =  Uv_0[ig];
      Hvat  =  Hv_0[ig];
      Rvat  =  Rv_0[ig]; Dvat = Dv_0[ig];

      FOI = 0;
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
        ig2 = is2*na + ia2;
        icm = ig2*cmdim1 + ig;
        cmi = cm[icm];
        FOI += beta*ua*cmi*(
          I_s[ig2]  + f*U_s[ig2] +
          Iv_s[ig2] + f*Uv_s[ig2]
        ) * oNg[ig2];
      }}

      FOIS   = FOI*Sat;
      FOIvV  = FOI*Vat*(1.0 - ve_i);

      rEo  = rEI*Eat;
      rUo  = rUR[ia]*Uat;
      rIo  = rIR[ia]*Iat;
      rHo  = rH*Hat;
      rEvo = rEI*Evat;
      rUvo = rUR[ia]*Uvat;
      rIvo = rIR[ia]*Ivat;
      rHvo = rH*Hvat;

      dEin  = dt*FOIS;
      dIin  = dt*yas*rEo;
      dUin  = dt*(1-yas)*rEo;
      dHin  = dt*ha*rIo;
      dVin  = dt*rVas*Sat;
      dEvin = dt*FOIvV;
      dIvin = dt*yas*(1.0 - ve_y)*rEvo;
      dUvin = dt*(1.0 - yas*(1.0 - ve_y))*rEvo;
      dHvin = dt*ha*(1.0 - ve_h)*rIvo;

      dS  = dt*( -FOIS - rVas*Sat + rW*Vat + rW_nat*Rat + rW_nat*Rvat );
      dE  = dt*(  FOIS - rEo );
      dI  = dt*(  yas*rEo - rIo );
      dU  = dt*(  (1-yas)*rEo - rUo );
      dH  = dt*(  ha*rIo - rHo );
      dR  = dt*(  (1-ha)*rIo + rUo + (1.0-mHa)*rHo - rW_nat*Rat );
      dD  = dt*(  mHa*rHo );
      dCc = dt*(  yas*rEo*rrepa );

      dV   = dt*( rVas*Sat - FOIvV - rW*Vat );
      dEv  = dt*(  FOIvV - rEvo );
      dIv  = dt*(  yas*(1.0 - ve_y)*rEvo - rIvo );
      dUv  = dt*(  (1.0 - yas*(1.0 - ve_y))*rEvo - rUvo );
      dHv  = dt*(  ha*(1.0 - ve_h)*rIvo - rHvo );
      dRv  = dt*(  (1.0 - ha*(1.0 - ve_h))*rIvo + rUvo + (1.0 - mHa*(1.0 - ve_m))*rHvo - rW_nat*Rvat );
      dDv  = dt*(  mHa*(1.0 - ve_m)*rHvo );
      dCcv = dt*(  yas*(1.0 - ve_y)*rEvo*rrepa );
      dVc  = dt*(  rVas*Sat );

      S_0[ig]  = Sat  + dS;
      E_0[ig]  = Eat  + dE;
      I_0[ig]  = Iat  + dI;
      U_0[ig]  = Uat  + dU;
      H_0[ig]  = Hat  + dH;
      R_0[ig]  = Rat  + dR;
      D_0[ig]  = Dat  + dD;
      Cc_0[ig] = Cc_0[ig] + dCc;

      V_0[ig]   = Vat   + dV;
      Ev_0[ig]  = Evat  + dEv;
      Iv_0[ig]  = Ivat  + dIv;
      Uv_0[ig]  = Uvat  + dUv;
      Hv_0[ig]  = Hvat  + dHv;
      Rv_0[ig]  = Rvat  + dRv;
      Dv_0[ig]  = Dvat  + dDv;
      Vc_0[ig]  = Vc_0[ig]  + dVc;
      Ccv_0[ig] = Ccv_0[ig] + dCcv;

      time[it+1] = (it+1)*dt;
      week       = 1 + (int) time[it+1]/7;

      Spw  += dS;
      Epw  += dEin;   Epw_v += dEvin;
      Ipw  += dIin;   Ipw_v += dIvin;
      Upw  += dUin;   Upw_v += dUvin;
      Hpw  += dHin;   Hpw_v += dHvin;
      Dpw  += dD;     Dpw_v += dDv;
      Rpw  += dR + dRv;
      Ccpw += dCc;
      Ccpw_v += dCcv;
      Vpw  += dVin;
      Epw_s[is] += dEin + dEvin;
      Ipw_s[is] += dIin + dIvin;
      Upw_s[is] += dUin + dUvin;
      Hpw_s[is] += dHin + dHvin;
      Hpw_ag[ig]  += dHin;
      Hpw_vag[ig] += dHvin;
      Rpw_s[is] += dR + dRv;
      Ipw_a[ia] += dIin + dIvin;
      Upw_a[ia] += dUin + dUvin;
      Hpw_a[ia] += dHin + dHvin;
    }
    }

    if (week - week0 == 1) {
      iW[week-1]   = it+1;
      Sw[week-1]   = Spw;   Spw  = 0;
      Ew[week-1]   = Epw;   Epw  = 0;
      Iw[week-1]   = Ipw;   Ipw  = 0;
      Uw[week-1]   = Upw;   Upw  = 0;
      Hw[week-1]   = Hpw;   Hpw  = 0;
      Rw[week-1]   = Rpw;   Rpw  = 0;
      Dw[week-1]   = Dpw;   Dpw  = 0;
      Ccw[week-1]  = Ccpw;  Ccpw = 0;
      Ewv[week-1]  = Epw_v;  Epw_v = 0;
      Iwv[week-1]  = Ipw_v;  Ipw_v = 0;
      Uwv[week-1]  = Upw_v;  Upw_v = 0;
      Hwv[week-1]  = Hpw_v;  Hpw_v = 0;
      Dwv[week-1]  = Dpw_v;  Dpw_v = 0;
      Ccwv[week-1] = Ccpw_v; Ccpw_v = 0;
      Vw[week-1]   = Vpw;    Vpw   = 0;

      double Vtot = 0;
      for (int is = 0; is < ns; is++) {
        double Vs = 0;
        for (int ia = 0; ia < na; ia++) Vs += V_0[is*na + ia];
        Vlev_s(week-1, is) = Vs;
        Vtot += Vs;
      }
      Vlev[week-1] = Vtot;
      for (int ia = 0; ia < na; ia++) {
        double Va = 0;
        for (int is = 0; is < ns; is++) Va += V_0[is*na + ia];
        Vlev_a(week-1, ia) = Va;
      }

      for (int is = 0; is < ns; is++) {
        Ew_s(week-1, is) = Epw_s[is]; Epw_s[is] = 0;
        Iw_s(week-1, is) = Ipw_s[is]; Ipw_s[is] = 0;
        Uw_s(week-1, is) = Upw_s[is]; Upw_s[is] = 0;
        Hw_s(week-1, is) = Hpw_s[is]; Hpw_s[is] = 0;
        Rw_s(week-1, is) = Rpw_s[is]; Rpw_s[is] = 0;
      }
      for (int ia = 0; ia < na; ia++) {
        Iw_a(week-1, ia) = Ipw_a[ia]; Ipw_a[ia] = 0;
        Uw_a(week-1, ia) = Upw_a[ia]; Upw_a[ia] = 0;
        Hw_a(week-1, ia) = Hpw_a[ia]; Hpw_a[ia] = 0;
      }
      for (int ig2 = 0; ig2 < ng; ig2++) {
        Hw_ag(week-1, ig2)  = Hpw_ag[ig2];  Hpw_ag[ig2]  = 0;
        Hwv_ag(week-1, ig2) = Hpw_vag[ig2]; Hpw_vag[ig2] = 0;
      }
    }
  }

  Rcpp::List byw = Rcpp::List::create(
    Named("iW")    = iW,
    Named("time")  = time[iW],
    Named("Sw")    = Sw,
    Named("Ew")    = Ew,
    Named("Iw")    = Iw,
    Named("Uw")    = Uw,
    Named("Hw")    = Hw,
    Named("Rw")    = Rw,
    Named("Dw")    = Dw,
    Named("Ccw")   = Ccw,
    Named("Ewv")   = Ewv,
    Named("Iwv")   = Iwv,
    Named("Uwv")   = Uwv,
    Named("Hwv")   = Hwv,
    Named("Dwv")   = Dwv,
    Named("Ccwv")  = Ccwv,
    Named("Vw")    = Vw,
    Named("Vlev")  = Vlev,
    Named("Iw_s")  = Iw_s,
    Named("Uw_s")  = Uw_s,
    Named("Hw_s")  = Hw_s,
    Named("Ew_s")  = Ew_s,
    Named("Rw_s")  = Rw_s,
    Named("Vlev_s")= Vlev_s);
  Rcpp::List byaw = Rcpp::List::create(
    Named("Iw_a")  = Iw_a,
    Named("Uw_a")  = Uw_a,
    Named("Hw_a")  = Hw_a,
    Named("Hw_ag") = Hw_ag,    //weekly H by (age x IMD), unvacc chain
    Named("Hwv_ag")= Hwv_ag,   //weekly H by (age x IMD), vacc breakthrough chain
    Named("Vlev_a")= Vlev_a);

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
