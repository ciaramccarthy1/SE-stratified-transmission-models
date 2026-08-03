//SEIRDH (with hospitalisation H compartment) - daily incidence
//
// Compartments per (age x IMD):  S, E, I, H, U, R, D  (single-stage, no Erlang)
//   - I -> H at rate h[ia]*rIR (hospitalisation among clinical)
//   - I -> R at rate (1-h[ia])*rIR (clinical no-hospital recovery)
//   - H  -> D at rate mH[ia]*rH (mortality among hospitalised)
//   - H  -> R at rate (1-mH[ia])*rH (recovery from hospital)
//
// Note: H is assumed to NOT contribute to onward transmission (isolated).

#include <Rcpp.h>
using namespace Rcpp;
#include <array>
#include <cmath>

// [[Rcpp::export]]
List model(List parscpp) {

  //natural history parameters
  const double beta0( parscpp["beta"]);       //baseline transmission (seasonal multiplier applied per day)
  const double b1(   parscpp["b1"]);          //seasonal amplitude (Gaussian pulse, David rsvie)
  const double phi(  parscpp["phi"]);         //seasonal phase (fraction of year at peak)
  const double psi(  parscpp["psi"]);         //seasonal width (fraction of year)
  const double rEI(  parscpp["rEI"]);              //latency E -> I/U (single stage, no Erlang)
  const std::vector<double> rIR( parscpp["rIR"]);  //recovery I -> R/H rate, by age (David exposure-group avg)
  const std::vector<double> rUR( parscpp["rUR"]);  //recovery U -> R   rate, by age
  const double f(    parscpp["f"]);
  const double rH(   parscpp["rH"]);            //1/hospital stay length (per day)
  const double rW_nat(parscpp["rW_nat"]);       //natural waning rate R -> S (per day); 0 disables

  //integration parameters
  const double dt(   parscpp["dt"]);
  const int    nt(   parscpp["nt"]);
  const int    nd(   parscpp["nd"]);
  const int    na(   parscpp["na"]);
  const int    ns(   parscpp["nimd"]);
  const int    ng = na*ns;

  IntegerVector iW(nd);

  //vector parameters (per age, length na)
  const std::vector<double> y45(parscpp["y45"]);   //clinical fraction, length ng
  const std::vector<double>  u( parscpp["u"]);     //susceptibility, length na
  const std::vector<double>  h( parscpp["h"]);     //hospitalisation fraction (among clinical), length na
  const std::vector<double> mH( parscpp["mH"]);    //mortality fraction (among hospitalised), length na
  const std::vector<double> rrep(parscpp["rrep"]); //reporting rate, length na

  //demographic ageing (continuous): eta = ageing-out rate per band (eta[na-1] = death
  //rate from the top band); births = daily births into each IMD's youngest S.
  //All-zero eta & births disable ageing (default).
  const std::vector<double> eta(   parscpp["eta"]);     //length na
  const std::vector<double> births(parscpp["births"]);  //length ns

  //contact matrix
  const std::vector<double> cm( parscpp["cm"]);
  const int   cmdim1( parscpp["cmdim1"]);

  //initial states
  const std::vector<double> Sg0( parscpp["Sg0"]);
  const std::vector<double> Eg0( parscpp["Eg0"]);
  const std::vector<double> Ug0( parscpp["Ug0"]);
  const std::vector<double> Ig0( parscpp["Ig0"]);
  const std::vector<double> Rg0( parscpp["Rg0"]);
  const std::vector<double> Dg0( parscpp["Dg0"]);
  const std::vector<double> oNg( parscpp["oNg"]);

  //current-step states (size ng)
  std::vector<double> S_0(ng), E_0(ng);
  std::vector<double> I_0(ng);
  std::vector<double> U_0(ng);
  std::vector<double>  H_0(ng);
  std::vector<double>  R_0(ng),  D_0(ng);
  std::vector<double> Cc_0(ng);

  //daily outputs - totals
  NumericVector Sw(nd), Ew(nd), Uw(nd), Iw(nd), Hw(nd), Rw(nd), Dw(nd), Ccw(nd);
  //daily SE-stratified incidence (nd x ns)
  NumericMatrix Ew_s(nd, ns), Iw_s(nd, ns), Uw_s(nd, ns), Hw_s(nd, ns), Rw_s(nd, ns);
  //(age x IMD)-stratified hospitalisation incidence (rows=day, cols=ig=is*na+ia)
  NumericMatrix Hw_ag(nd, ng);
  //daily age-stratified incidence (nd x na)
  NumericMatrix Iw_a(nd, na), Uw_a(nd, na), Hw_a(nd, na);
  NumericMatrix Nw_a(nd, na);  //living population (S+E+I+U+H+R) by age band, per day

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
    Cc_0[ig] = 0.0;
  }}

  time[0]   = 0;
  iW[0]     = 0;
  int day   = 1;
  int day0  = 1;

  double Spw=0, Epw=0, Upw=0, Ipw=0, Hpw=0, Rpw=0, Dpw=0, Ccpw=0;
  NumericVector Epw_s(ns), Upw_s(ns), Ipw_s(ns), Hpw_s(ns), Rpw_s(ns);
  NumericVector Upw_a(na), Ipw_a(na), Hpw_a(na);
  NumericVector Hpw_ag(ng);  //daily H accumulator per (age x IMD)

  double yas, ua, ha, mHa, rrepa, cmi;
  double Sat, Eat, Uat, Iat, Hat, Rat, Dat;
  double FOI, FOIS;
  double rEo, rUo, rIo, rHo;
  double dEin, dIin, dUin, dHin;
  double dS, dE, dI, dU, dH, dR, dD, dCc;
  int    icm, ig2;

  for (int it = 0; it < (nt-1); it++) {
    day0 = day;
    day  = 1 + (int) time[it];

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

    //snapshot of current state for ageing (begin-of-step values -> exactly conserves people)
    std::vector<double> Sp=S_0, Ep=E_0, Ip=I_0, Up=U_0, Hp=H_0, Rp=R_0;

    for (int is = 0; is < ns; is++) {
    for (int ia = 0; ia < na; ia++) {
      ig    = is*na + ia;
      yas   = y45[ig];
      ua    = u[ia];
      ha    = h[ia];
      mHa   = mH[ia];
      rrepa = rrep[ia];

      Sat  =  S_0[ig];
      Eat  =  E_0[ig];
      Iat  =  I_0[ig];
      Uat  =  U_0[ig];
      Hat  =  H_0[ig];
      Rat  =  R_0[ig];
      Dat  =  D_0[ig];

      FOI = 0;
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
        ig2 = is2*na + ia2;
        icm = ig2*cmdim1 + ig;
        cmi = cm[icm];
        //H is not infectious (isolated). Read begin-of-step snapshot Ip/Up so
        //FOI is a Jacobi step, independent of is/ia order (I_0/U_0 updated in place).
        FOI += beta*ua*cmi*( Ip[ig2] + f*Up[ig2] )*oNg[ig2];
      }}

      FOIS = FOI*Sat;
      rEo  = rEI*Eat;
      rUo  = rUR[ia]*Uat;
      rIo  = rIR[ia]*Iat;
      rHo  = rH*Hat;
      dEin = dt*FOIS;
      dIin = dt*yas*rEo;
      dUin = dt*(1-yas)*rEo;
      dHin = dt*ha*rIo;

      dS  = dt*( -FOIS + rW_nat*Rat     );
      dE  = dt*(  FOIS - rEo            );
      dI  = dt*(  yas*rEo - rIo         );
      dU  = dt*(  (1-yas)*rEo - rUo     );
      dH  = dt*(  ha*rIo - rHo          );
      dR  = dt*(  (1-ha)*rIo + rUo + (1-mHa)*rHo - rW_nat*Rat );
      dD  = dt*(  mHa*rHo               );
      dCc = dt*(  yas*rEo*rrepa         );

      // --- demographic ageing: eta[ia] out (to ia+1; death from the top band), eta[ia-1] in ---
      dS -= dt*eta[ia]*Sat; dE -= dt*eta[ia]*Eat; dI -= dt*eta[ia]*Iat;
      dU -= dt*eta[ia]*Uat; dH -= dt*eta[ia]*Hat; dR -= dt*eta[ia]*Rat;
      if (ia > 0) {                 // ageing in from band below (same IMD); begin-of-step values
        double r = eta[ia-1];
        dS += dt*r*Sp[ig-1]; dE += dt*r*Ep[ig-1]; dI += dt*r*Ip[ig-1];
        dU += dt*r*Up[ig-1]; dH += dt*r*Hp[ig-1]; dR += dt*r*Rp[ig-1];
      } else {                      // youngest band: births enter susceptible
        dS += dt*births[is];
      }

      S_0[ig]  = Sat  + dS;
      E_0[ig]  = Eat  + dE;
      I_0[ig]  = Iat  + dI;
      U_0[ig]  = Uat  + dU;
      H_0[ig]  = Hat  + dH;
      R_0[ig]  = Rat  + dR;
      D_0[ig]  = Dat  + dD;
      Cc_0[ig] = Cc_0[ig] + dCc;

      time[it+1] = (it+1)*dt;
      day        = 1 + (int) time[it+1];

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

    if (day - day0 == 1) {
      iW[day-1]  = it+1;
      Sw[day-1]  = Spw;  Spw  = 0;
      Ew[day-1]  = Epw;  Epw  = 0;
      Iw[day-1]  = Ipw;  Ipw  = 0;
      Uw[day-1]  = Upw;  Upw  = 0;
      Hw[day-1]  = Hpw;  Hpw  = 0;
      Rw[day-1]  = Rpw;  Rpw  = 0;
      Dw[day-1]  = Dpw;  Dpw  = 0;
      Ccw[day-1] = Ccpw; Ccpw = 0;
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
        double n = 0;                                        // living pop by age band
        for (int is = 0; is < ns; is++) { int g = is*na+ia; n += S_0[g]+E_0[g]+I_0[g]+U_0[g]+H_0[g]+R_0[g]; }
        Nw_a(day-1, ia) = n;
      }
      for (int ig2 = 0; ig2 < ng; ig2++) {
        Hw_ag(day-1, ig2) = Hpw_ag[ig2]; Hpw_ag[ig2] = 0;
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
    Named("Nw_a") = Nw_a,     //living population by age band, per day
    Named("Hw_ag") = Hw_ag);  //daily H by (age x IMD)

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
