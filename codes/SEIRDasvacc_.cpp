//SEIRDHV - weekly incidence
// See SEIRDasvaccday_.cpp for compartment / equation documentation.

#include <Rcpp.h>
using namespace Rcpp;
#include <array>

// [[Rcpp::export]]
List model(List parscpp) {

  const double beta( parscpp["beta"]);
  const double rEI(  parscpp["rEI"]);
  const double rI1I2(parscpp["rI1I2"]);
  const double rI2R( parscpp["rI2R"]);
  const double rUR(  parscpp["rUR"]);
  const double rE    = 2*rEI;
  const double rU    = 2*rUR;
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
  const std::vector<double> E1g0(parscpp["E1g0"]);
  const std::vector<double> E2g0(parscpp["E2g0"]);
  const std::vector<double> U1g0(parscpp["U1g0"]);
  const std::vector<double> U2g0(parscpp["U2g0"]);
  const std::vector<double> I1g0(parscpp["I1g0"]);
  const std::vector<double> I2g0(parscpp["I2g0"]);
  const std::vector<double> Rg0( parscpp["Rg0"]);
  const std::vector<double> Dg0( parscpp["Dg0"]);
  const std::vector<double> oNg( parscpp["oNg"]);

  std::vector<double> S_0(ng), E1_0(ng), E2_0(ng);
  std::vector<double> I1_0(ng), I2_0(ng), H_0(ng);
  std::vector<double> U1_0(ng), U2_0(ng);
  std::vector<double> R_0(ng),  D_0(ng);
  std::vector<double> V_0(ng);
  std::vector<double> Ev1_0(ng), Ev2_0(ng);
  std::vector<double> Iv1_0(ng), Iv2_0(ng), Hv_0(ng);
  std::vector<double> Uv1_0(ng), Uv2_0(ng);
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
    E1_0[ig] = E1g0[ig]; E2_0[ig] = E2g0[ig];
    I1_0[ig] = I1g0[ig]; I2_0[ig] = I2g0[ig];
    U1_0[ig] = U1g0[ig]; U2_0[ig] = U2g0[ig];
    H_0[ig]  = 0.0;
    R_0[ig]  = Rg0[ig];
    D_0[ig]  = Dg0[ig];
    V_0[ig] = 0.0;
    Ev1_0[ig] = 0.0; Ev2_0[ig] = 0.0;
    Iv1_0[ig] = 0.0; Iv2_0[ig] = 0.0;
    Uv1_0[ig] = 0.0; Uv2_0[ig] = 0.0;
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
  double Sat, E1at, E2at, U1at, U2at, I1at, I2at, Hat, Rat, Dat;
  double Vat, Ev1at, Ev2at, Uv1at, Uv2at, Iv1at, Iv2at, Hvat, Rvat, Dvat;
  double FOI, FOIS, FOIvV;
  double rE1, rE2, rU1, rU2, rI1, rI2, rHo;
  double rEv1, rEv2, rUv1, rUv2, rIv1, rIv2, rHvo;
  double dVin, dEin, dIin, dUin, dHin;
  double dEvin, dIvin, dUvin, dHvin;
  double dS, dE1, dE2, dI1, dI2, dH, dU1, dU2, dR, dD, dCc;
  double dV, dEv1, dEv2, dIv1, dIv2, dHv, dUv1, dUv2, dRv, dDv, dCcv, dVc;
  int    icm, ig2;

  for (int it = 0; it < (nt-1); it++) {
    week0 = week;
    week  = 1 + (int) time[it]/7;

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

      dEin  = dt*FOIS;
      dIin  = dt*yas*rE2;
      dUin  = dt*(1-yas)*rE2;
      dHin  = dt*ha*rI2;
      dVin  = dt*rVas*Sat;
      dEvin = dt*FOIvV;
      dIvin = dt*yas*(1.0 - ve_y)*rEv2;
      dUvin = dt*(1.0 - yas*(1.0 - ve_y))*rEv2;
      dHvin = dt*ha*(1.0 - ve_h)*rIv2;

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
      dVc  = dt*(  rVas*Sat );

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
