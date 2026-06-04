//SEIRDH (with hospitalisation H compartment) - weekly incidence
// See SEIRDasday_.cpp for compartment / equation documentation.

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
  const double rW_nat(parscpp["rW_nat"]);       //natural waning rate R -> S (per day); 0 disables

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
  std::vector<double> I1_0(ng), I2_0(ng);
  std::vector<double> U1_0(ng), U2_0(ng);
  std::vector<double>  H_0(ng);
  std::vector<double>  R_0(ng),  D_0(ng);
  std::vector<double> Cc_0(ng);

  NumericVector Sw(nw), Ew(nw), Uw(nw), Iw(nw), Hw(nw), Rw(nw), Dw(nw), Ccw(nw);
  NumericMatrix Ew_s(nw, ns), Iw_s(nw, ns), Uw_s(nw, ns), Hw_s(nw, ns), Rw_s(nw, ns);
  //(age x IMD)-stratified hospitalisation incidence (rows=week, cols=ig=is*na+ia)
  NumericMatrix Hw_ag(nw, ng);
  NumericMatrix Iw_a(nw, na), Uw_a(nw, na), Hw_a(nw, na);

  NumericVector time(nt);
  int ig;

  for (int is = 0; is < ns; is++) {
  for (int ia = 0; ia < na; ia++) {
    ig = ia + is*na;
    S_0[ig]  = Sg0[ig];
    E1_0[ig] = E1g0[ig];
    E2_0[ig] = E2g0[ig];
    I1_0[ig] = I1g0[ig];
    I2_0[ig] = I2g0[ig];
    U1_0[ig] = U1g0[ig];
    U2_0[ig] = U2g0[ig];
    H_0[ig]  = 0.0;
    R_0[ig]  = Rg0[ig];
    D_0[ig]  = Dg0[ig];
    Cc_0[ig] = 0.0;
  }}

  time[0]    = 0;
  iW[0]      = 0;
  int week   = 1;
  int week0  = 1;

  double Spw=0, Epw=0, Upw=0, Ipw=0, Hpw=0, Rpw=0, Dpw=0, Ccpw=0;
  NumericVector Epw_s(ns), Upw_s(ns), Ipw_s(ns), Hpw_s(ns), Rpw_s(ns);
  NumericVector Upw_a(na), Ipw_a(na), Hpw_a(na);
  NumericVector Hpw_ag(ng);  //weekly H accumulator per (age x IMD)

  double yas, ua, ha, mHa, rrepa, cmi;
  double Sat, E1at, E2at, U1at, U2at, I1at, I2at, Hat, Rat, Dat;
  double FOI, FOIS;
  double rE1, rE2, rU1, rU2, rI1, rI2, rHo;
  double dEin, dIin, dUin, dHin;
  double dS, dE1, dE2, dI1, dI2, dU1, dU2, dH, dR, dD, dCc;
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

      Sat  =  S_0[ig];
      E1at = E1_0[ig];
      E2at = E2_0[ig];
      I1at = I1_0[ig];
      I2at = I2_0[ig];
      U1at = U1_0[ig];
      U2at = U2_0[ig];
      Hat  =  H_0[ig];
      Rat  =  R_0[ig];
      Dat  =  D_0[ig];

      FOI = 0;
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
        ig2 = is2*na + ia2;
        icm = ig2*cmdim1 + ig;
        cmi = cm[icm];
        FOI += beta*ua*cmi*( I1_0[ig2] + I2_0[ig2] + f*U1_0[ig2] + f*U2_0[ig2] )*oNg[ig2];
      }}

      FOIS = FOI*Sat;
      rE1  = rE*E1at;    rE2 = rE*E2at;
      rU1  = rU*U1at;    rU2 = rU*U2at;
      rI1  = rI1I2*I1at; rI2 = rI2R*I2at;
      rHo  = rH*Hat;
      dEin = dt*FOIS;
      dIin = dt*yas*rE2;
      dUin = dt*(1-yas)*rE2;
      dHin = dt*ha*rI2;

      dS  = dt*( -FOIS + rW_nat*Rat     );
      dE1 = dt*(  FOIS - rE1            );
      dE2 = dt*(  rE1  - rE2            );
      dI1 = dt*(  yas*rE2 - rI1         );
      dI2 = dt*(  rI1  - rI2            );
      dU1 = dt*(  (1-yas)*rE2 - rU1     );
      dU2 = dt*(  rU1  - rU2            );
      dH  = dt*(  ha*rI2 - rHo          );
      dR  = dt*(  (1-ha)*rI2 + rU2 + (1-mHa)*rHo - rW_nat*Rat );
      dD  = dt*(  mHa*rHo               );
      dCc = dt*(  yas*rE2*rrepa         );

      S_0[ig]  = Sat  + dS;
      E1_0[ig] = E1at + dE1;
      E2_0[ig] = E2at + dE2;
      I1_0[ig] = I1at + dI1;
      I2_0[ig] = I2at + dI2;
      U1_0[ig] = U1at + dU1;
      U2_0[ig] = U2at + dU2;
      H_0[ig]  = Hat  + dH;
      R_0[ig]  = Rat  + dR;
      D_0[ig]  = Dat  + dD;
      Cc_0[ig] = Cc_0[ig] + dCc;

      time[it+1] = (it+1)*dt;
      week       = 1 + (int) time[it+1]/7;

      Spw  += dS;
      Epw  += dEin;
      Ipw  += dIin;
      Upw  += dUin;
      Hpw  += dHin;
      Rpw  += dR;
      Dpw  += dD;
      Ccpw += dCc;
      Epw_s[is] += dEin;
      Ipw_s[is] += dIin;
      Upw_s[is] += dUin;
      Hpw_s[is] += dHin;
      Hpw_ag[ig] += dHin;
      Rpw_s[is] += dR;
      Ipw_a[ia] += dIin;
      Upw_a[ia] += dUin;
      Hpw_a[ia] += dHin;
    }
    }

    if (week - week0 == 1) {
      iW[week-1]  = it+1;
      Sw[week-1]  = Spw;  Spw  = 0;
      Ew[week-1]  = Epw;  Epw  = 0;
      Iw[week-1]  = Ipw;  Ipw  = 0;
      Uw[week-1]  = Upw;  Upw  = 0;
      Hw[week-1]  = Hpw;  Hpw  = 0;
      Rw[week-1]  = Rpw;  Rpw  = 0;
      Dw[week-1]  = Dpw;  Dpw  = 0;
      Ccw[week-1] = Ccpw; Ccpw = 0;
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
        Hw_ag(week-1, ig2) = Hpw_ag[ig2]; Hpw_ag[ig2] = 0;
      }
    }
  }

  Rcpp::List byw = Rcpp::List::create(
    Named("iW")   = iW,
    Named("time") = time[iW],
    Named("Sw")   = Sw,
    Named("Ew")   = Ew,
    Named("Iw")   = Iw,
    Named("Uw")   = Uw,
    Named("Hw")   = Hw,
    Named("Dw")   = Dw,
    Named("Ccw")  = Ccw,
    Named("Iw_s") = Iw_s,
    Named("Uw_s") = Uw_s,
    Named("Hw_s") = Hw_s,
    Named("Ew_s") = Ew_s,
    Named("Rw_s") = Rw_s);
  Rcpp::List byaw = Rcpp::List::create(
    Named("Iw_a") = Iw_a,
    Named("Uw_a") = Uw_a,
    Named("Hw_a") = Hw_a,
    Named("Hw_ag") = Hw_ag);  //weekly H by (age x IMD)

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
