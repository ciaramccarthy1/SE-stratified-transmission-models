//SEIRDHV - daily incidence
// Leaky vacc compartment V with breakthrough shadow chain (Ev/Iv/Uv/Hv/Rv/Dv)
// and 4 vaccine-efficacy knobs (per (age x IMD)):
//   VE_inf  - reduces susceptibility of V (leaky)
//   VE_sym  - reduces clinical fraction (E2_v -> I1_v vs U1_v)
//   VE_hosp - reduces hospitalisation fraction (I2_v -> H_v vs R_v)
//   VE_sev  - reduces mortality given hospitalised (H_v -> D vs R_v)
//
// Compartments per (age x IMD):
//   Unvacc: S, E1, E2, I1, I2, H, U1, U2, R, D
//   Vacc  : V, Ev1, Ev2, Iv1, Iv2, Hv, Uv1, Uv2, Rv, Dv
//   Counters: Vc (cumulative doses), Cc (reported clinical, unvacc),
//             Ccv (reported clinical, vacc breakthroughs)
//
// Vaccination: S -> V at rate rV*vcov (no efficacy split at vaccination - leaky)
// Waning: V -> S at rate rW; R/Rv -> S at rate rW_nat
// H and Hv assumed isolated (no FOI contribution)
//
// Force of infection on group ig sums contributions from BOTH chains:
//   I1, I2 (unvacc), Iv1, Iv2 (vacc)
//   U1, U2 (unvacc), Uv1, Uv2 (vacc) - subclinical, scaled by f
// Vaccinated infectious individuals transmit at the SAME rate as unvaccinated
// (no VE_transmission knob in this version).

#include <Rcpp.h>
using namespace Rcpp;
#include <array>

// [[Rcpp::export]]
List model(List parscpp) {

  //natural history
  const double beta( parscpp["beta"]);
  const double rEI(  parscpp["rEI"]);
  const double rI1I2(parscpp["rI1I2"]);
  const double rI2R( parscpp["rI2R"]);
  const double rUR(  parscpp["rUR"]);
  const double rE    = 2*rEI;
  const double rU    = 2*rUR;
  const double f(    parscpp["f"]);
  const double rH(   parscpp["rH"]);

  //vacc + waning
  const NumericVector vcov(   parscpp["vcov"]);    //coverage by age x IMD, length ng
  const NumericVector VE_inf( parscpp["VE_inf"]);  //efficacy against infection, length ng
  const NumericVector VE_sym( parscpp["VE_sym"]);  //efficacy against symptoms, length ng
  const NumericVector VE_hosp(parscpp["VE_hosp"]); //efficacy against hospitalisation, length ng
  const NumericVector VE_sev( parscpp["VE_sev"]);  //efficacy against mortality given hosp, length ng
  const double rV(     parscpp["rV"]);             //vacc rate (per day)
  const double rW(     parscpp["rW"]);             //vaccine waning rate V -> S
  const double rW_nat( parscpp["rW_nat"]);         //natural waning rate R/Rv -> S (0 disables)

  //integration
  const double dt(   parscpp["dt"]);
  const int    nt(   parscpp["nt"]);
  const int    nd(   parscpp["nd"]);
  const int    na(   parscpp["na"]);
  const int    ns(   parscpp["nimd"]);
  const int    ng = na*ns;

  IntegerVector iW(nd);

  const std::vector<double> y45(parscpp["y45"]);
  const std::vector<double>  u( parscpp["u"]);
  const std::vector<double>  h( parscpp["h"]);
  const std::vector<double> mH( parscpp["mH"]);
  const std::vector<double> rrep(parscpp["rrep"]);

  const std::vector<double> cm( parscpp["cm"]);
  const int   cmdim1( parscpp["cmdim1"]);

  const std::vector<double> Sg0( parscpp["Sg0"]);
  const std::vector<double> E1g0(parscpp["E1g0"]);
  const std::vector<double> E2g0(parscpp["E1g0"]);
  const std::vector<double> U1g0(parscpp["U1g0"]);
  const std::vector<double> U2g0(parscpp["U2g0"]);
  const std::vector<double> I1g0(parscpp["I1g0"]);
  const std::vector<double> I2g0(parscpp["I2g0"]);
  const std::vector<double> Rg0( parscpp["Rg0"]);
  const std::vector<double> Dg0( parscpp["Dg0"]);
  const std::vector<double> oNg( parscpp["oNg"]);

  //unvacc states
  std::vector<double> S_0(ng), E1_0(ng), E2_0(ng);
  std::vector<double> I1_0(ng), I2_0(ng), H_0(ng);
  std::vector<double> U1_0(ng), U2_0(ng);
  std::vector<double> R_0(ng),  D_0(ng);
  //vacc states
  std::vector<double> V_0(ng);
  std::vector<double> Ev1_0(ng), Ev2_0(ng);
  std::vector<double> Iv1_0(ng), Iv2_0(ng), Hv_0(ng);
  std::vector<double> Uv1_0(ng), Uv2_0(ng);
  std::vector<double> Rv_0(ng),  Dv_0(ng);
  //counters
  std::vector<double> Vc_0(ng), Cc_0(ng), Ccv_0(ng);

  //daily outputs - totals (vacc + unvacc unless otherwise noted)
  NumericVector Sw(nd), Ew(nd), Iw(nd), Uw(nd), Hw(nd), Rw(nd), Dw(nd);
  NumericVector Ccw(nd);     //reported clinical, unvacc
  NumericVector Ewv(nd), Iwv(nd), Uwv(nd), Hwv(nd), Dwv(nd), Ccwv(nd); //vacc-only flows
  NumericVector Vw(nd);      //daily new vaccinations (= rVas*S)
  NumericVector Vlev(nd);    //V compartment level at end of day (stock)
  //stratified matrices (totals = unvacc + vacc)
  NumericMatrix Ew_s(nd, ns), Iw_s(nd, ns), Uw_s(nd, ns), Hw_s(nd, ns), Rw_s(nd, ns);
  NumericMatrix Iw_a(nd, na), Uw_a(nd, na), Hw_a(nd, na);
  //V-specific stratification (stock)
  NumericMatrix Vlev_a(nd, na), Vlev_s(nd, ns);

  NumericVector time(nt);
  int ig;

  for (int is = 0; is < ns; is++) {
  for (int ia = 0; ia < na; ia++) {
    ig = ia + is*na;
    //unvacc
    S_0[ig]  = Sg0[ig];
    E1_0[ig] = E1g0[ig]; E2_0[ig] = E2g0[ig];
    I1_0[ig] = I1g0[ig]; I2_0[ig] = I2g0[ig];
    U1_0[ig] = U1g0[ig]; U2_0[ig] = U2g0[ig];
    H_0[ig]  = 0.0;
    R_0[ig]  = Rg0[ig];
    D_0[ig]  = Dg0[ig];
    //vacc (all start empty)
    V_0[ig] = 0.0;
    Ev1_0[ig] = 0.0; Ev2_0[ig] = 0.0;
    Iv1_0[ig] = 0.0; Iv2_0[ig] = 0.0;
    Uv1_0[ig] = 0.0; Uv2_0[ig] = 0.0;
    Hv_0[ig]  = 0.0;
    Rv_0[ig]  = 0.0; Dv_0[ig] = 0.0;
    Vc_0[ig]  = 0.0; Cc_0[ig] = 0.0; Ccv_0[ig] = 0.0;
  }}

  time[0]   = 0;
  iW[0]     = 0;
  int day   = 1;
  int day0  = 1;

  //per-day accumulators (combine unvacc + vacc for totals; track vacc separately too)
  double Spw=0, Epw=0, Ipw=0, Upw=0, Hpw=0, Rpw=0, Dpw=0;
  double Ccpw=0;
  double Epw_v=0, Ipw_v=0, Upw_v=0, Hpw_v=0, Dpw_v=0, Ccpw_v=0;
  double Vpw=0;  //daily new vaccinations
  NumericVector Epw_s(ns), Ipw_s(ns), Upw_s(ns), Hpw_s(ns), Rpw_s(ns);
  NumericVector Ipw_a(na), Upw_a(na), Hpw_a(na);

  //temp working vars
  double yas, ua, ha, mHa, rrepa, cmi;
  double rVas, ve_i, ve_y, ve_h, ve_m;
  //unvacc temps
  double Sat, E1at, E2at, U1at, U2at, I1at, I2at, Hat, Rat, Dat;
  //vacc temps
  double Vat, Ev1at, Ev2at, Uv1at, Uv2at, Iv1at, Iv2at, Hvat, Rvat, Dvat;
  //flows
  double FOI, FOIS, FOIvV;
  double rE1, rE2, rU1, rU2, rI1, rI2, rHo;
  double rEv1, rEv2, rUv1, rUv2, rIv1, rIv2, rHvo;
  double dVin, dEin, dIin, dUin, dHin;
  double dEvin, dIvin, dUvin, dHvin;
  double dS, dE1, dE2, dI1, dI2, dH, dU1, dU2, dR, dD, dCc;
  double dV, dEv1, dEv2, dIv1, dIv2, dHv, dUv1, dUv2, dRv, dDv, dCcv, dVc;
  int    icm, ig2;

  for (int it = 0; it < (nt-1); it++) {
    day0 = day;
    day  = 1 + (int) time[it];

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
      ve_m  = VE_sev[ig];

      Sat  =  S_0[ig];
      E1at = E1_0[ig]; E2at = E2_0[ig];
      I1at = I1_0[ig]; I2at = I2_0[ig];
      U1at = U1_0[ig]; U2at = U2_0[ig];
      Hat  =  H_0[ig];
      Rat  =  R_0[ig]; Dat  =  D_0[ig];

      Vat   =  V_0[ig];
      Ev1at = Ev1_0[ig]; Ev2at = Ev2_0[ig];
      Iv1at = Iv1_0[ig]; Iv2at = Iv2_0[ig];
      Uv1at = Uv1_0[ig]; Uv2at = Uv2_0[ig];
      Hvat  =  Hv_0[ig];
      Rvat  =  Rv_0[ig]; Dvat = Dv_0[ig];

      //FOI: sum across unvacc + vacc infectious (H/Hv excluded - isolated)
      FOI = 0;
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
        ig2 = is2*na + ia2;
        icm = ig2*cmdim1 + ig;
        cmi = cm[icm];
        FOI += beta*ua*cmi*(
          I1_0[ig2]  + I2_0[ig2]  + f*U1_0[ig2]  + f*U2_0[ig2] +
          Iv1_0[ig2] + Iv2_0[ig2] + f*Uv1_0[ig2] + f*Uv2_0[ig2]
        ) * oNg[ig2];
      }}

      FOIS   = FOI*Sat;
      FOIvV  = FOI*Vat*(1.0 - ve_i);

      rE1  = rE*E1at;     rE2  = rE*E2at;
      rU1  = rU*U1at;     rU2  = rU*U2at;
      rI1  = rI1I2*I1at;  rI2  = rI2R*I2at;
      rHo  = rH*Hat;
      rEv1 = rE*Ev1at;    rEv2 = rE*Ev2at;
      rUv1 = rU*Uv1at;    rUv2 = rU*Uv2at;
      rIv1 = rI1I2*Iv1at; rIv2 = rI2R*Iv2at;
      rHvo = rH*Hvat;

      //daily incidence flows (for reporting)
      dEin  = dt*FOIS;
      dIin  = dt*yas*rE2;
      dUin  = dt*(1-yas)*rE2;
      dHin  = dt*ha*rI2;
      dVin  = dt*rVas*Sat;                              //new vaccinations
      dEvin = dt*FOIvV;
      dIvin = dt*yas*(1.0 - ve_y)*rEv2;
      dUvin = dt*(1.0 - yas*(1.0 - ve_y))*rEv2;
      dHvin = dt*ha*(1.0 - ve_h)*rIv2;

      //--- unvacc chain ---
      dS  = dt*( -FOIS - rVas*Sat + rW*Vat + rW_nat*Rat + rW_nat*Rvat );
      dE1 = dt*(  FOIS - rE1 );
      dE2 = dt*(  rE1  - rE2 );
      dI1 = dt*(  yas*rE2 - rI1 );
      dI2 = dt*(  rI1  - rI2 );
      dU1 = dt*(  (1-yas)*rE2 - rU1 );
      dU2 = dt*(  rU1  - rU2 );
      dH  = dt*(  ha*rI2 - rHo );
      dR  = dt*(  (1-ha)*rI2 + rU2 + (1.0-mHa)*rHo - rW_nat*Rat );
      dD  = dt*(  mHa*rHo );
      dCc = dt*(  yas*rE2*rrepa );

      //--- vacc chain ---
      dV   = dt*( rVas*Sat - FOIvV - rW*Vat );
      dEv1 = dt*(  FOIvV - rEv1 );
      dEv2 = dt*(  rEv1  - rEv2 );
      dIv1 = dt*(  yas*(1.0 - ve_y)*rEv2 - rIv1 );
      dIv2 = dt*(  rIv1  - rIv2 );
      dUv1 = dt*(  (1.0 - yas*(1.0 - ve_y))*rEv2 - rUv1 );
      dUv2 = dt*(  rUv1  - rUv2 );
      dHv  = dt*(  ha*(1.0 - ve_h)*rIv2 - rHvo );
      dRv  = dt*(  (1.0 - ha*(1.0 - ve_h))*rIv2 + rUv2 + (1.0 - mHa*(1.0 - ve_m))*rHvo - rW_nat*Rvat );
      dDv  = dt*(  mHa*(1.0 - ve_m)*rHvo );
      dCcv = dt*(  yas*(1.0 - ve_y)*rEv2*rrepa );
      dVc  = dt*(  rVas*Sat );  //cumulative doses

      //--- update ---
      S_0[ig]  = Sat  + dS;
      E1_0[ig] = E1at + dE1; E2_0[ig] = E2at + dE2;
      I1_0[ig] = I1at + dI1; I2_0[ig] = I2at + dI2;
      U1_0[ig] = U1at + dU1; U2_0[ig] = U2at + dU2;
      H_0[ig]  = Hat  + dH;
      R_0[ig]  = Rat  + dR;
      D_0[ig]  = Dat  + dD;
      Cc_0[ig] = Cc_0[ig] + dCc;

      V_0[ig]   = Vat   + dV;
      Ev1_0[ig] = Ev1at + dEv1; Ev2_0[ig] = Ev2at + dEv2;
      Iv1_0[ig] = Iv1at + dIv1; Iv2_0[ig] = Iv2at + dIv2;
      Uv1_0[ig] = Uv1at + dUv1; Uv2_0[ig] = Uv2at + dUv2;
      Hv_0[ig]  = Hvat  + dHv;
      Rv_0[ig]  = Rvat  + dRv;
      Dv_0[ig]  = Dvat  + dDv;
      Vc_0[ig]  = Vc_0[ig]  + dVc;
      Ccv_0[ig] = Ccv_0[ig] + dCcv;

      //time
      time[it+1] = (it+1)*dt;
      day        = 1 + (int) time[it+1];

      //accumulators
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
      //totals stratification (unvacc + vacc)
      Epw_s[is] += dEin + dEvin;
      Ipw_s[is] += dIin + dIvin;
      Upw_s[is] += dUin + dUvin;
      Hpw_s[is] += dHin + dHvin;
      Rpw_s[is] += dR + dRv;
      Ipw_a[ia] += dIin + dIvin;
      Upw_a[ia] += dUin + dUvin;
      Hpw_a[ia] += dHin + dHvin;
    }
    }

    if (day - day0 == 1) {
      iW[day-1]   = it+1;
      Sw[day-1]   = Spw;   Spw  = 0;
      Ew[day-1]   = Epw;   Epw  = 0;
      Iw[day-1]   = Ipw;   Ipw  = 0;
      Uw[day-1]   = Upw;   Upw  = 0;
      Hw[day-1]   = Hpw;   Hpw  = 0;
      Rw[day-1]   = Rpw;   Rpw  = 0;
      Dw[day-1]   = Dpw;   Dpw  = 0;
      Ccw[day-1]  = Ccpw;  Ccpw = 0;
      Ewv[day-1]  = Epw_v;  Epw_v = 0;
      Iwv[day-1]  = Ipw_v;  Ipw_v = 0;
      Uwv[day-1]  = Upw_v;  Upw_v = 0;
      Hwv[day-1]  = Hpw_v;  Hpw_v = 0;
      Dwv[day-1]  = Dpw_v;  Dpw_v = 0;
      Ccwv[day-1] = Ccpw_v; Ccpw_v = 0;
      Vw[day-1]   = Vpw;    Vpw   = 0;

      //V level stock
      double Vtot = 0;
      for (int is = 0; is < ns; is++) {
        double Vs = 0;
        for (int ia = 0; ia < na; ia++) Vs += V_0[is*na + ia];
        Vlev_s(day-1, is) = Vs;
        Vtot += Vs;
      }
      Vlev[day-1] = Vtot;
      for (int ia = 0; ia < na; ia++) {
        double Va = 0;
        for (int is = 0; is < ns; is++) Va += V_0[is*na + ia];
        Vlev_a(day-1, ia) = Va;
      }

      for (int is = 0; is < ns; is++) {
        Ew_s(day-1, is) = Epw_s[is]; Epw_s[is] = 0;
        Iw_s(day-1, is) = Ipw_s[is]; Ipw_s[is] = 0;
        Uw_s(day-1, is) = Upw_s[is]; Upw_s[is] = 0;
        Hw_s(day-1, is) = Hpw_s[is]; Hpw_s[is] = 0;
        Rw_s(day-1, is) = Rpw_s[is]; Rpw_s[is] = 0;
      }
      for (int ia = 0; ia < na; ia++) {
        Iw_a(day-1, ia) = Ipw_a[ia]; Ipw_a[ia] = 0;
        Uw_a(day-1, ia) = Upw_a[ia]; Upw_a[ia] = 0;
        Hw_a(day-1, ia) = Hpw_a[ia]; Hpw_a[ia] = 0;
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
    Named("Vlev_a")= Vlev_a);

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
