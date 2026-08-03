//SEIRDHV - daily incidence  (single-stage compartments, no Erlang)
// Leaky vacc compartment V with breakthrough shadow chain (Ev/Iv/Uv/Hv/Rv/Dv)
// and 4 vaccine-efficacy knobs (per (age x IMD)):
//   VE_inf  - reduces susceptibility of V (leaky)
//   VE_sym  - reduces clinical fraction (E_v -> I_v vs U_v)
//   VE_hosp - reduces hospitalisation fraction (I_v -> H_v vs R_v)
//   VE_mort  - reduces mortality given hospitalised (H_v -> D vs R_v)
//
// Compartments per (age x IMD):
//   Unvacc: S, E, I, H, U, R, D
//   Vacc  : V, Ev, Iv, Hv, Uv, Rv, Dv
//   Counters: Vc (cumulative doses), Cc (reported clinical, unvacc),
//             Ccv (reported clinical, vacc breakthroughs)
//
// Vaccination: S -> V at rate rV*vcov (no efficacy split at vaccination - leaky)
// Waning: V -> S at rate rW; R/Rv -> S at rate rW_nat
// H and Hv assumed isolated (no FOI contribution)
//
// Force of infection on group ig sums contributions from BOTH chains:
//   I (unvacc), Iv (vacc)
//   U (unvacc), Uv (vacc) - subclinical, scaled by f
// Vaccinated infectious individuals transmit at the SAME rate as unvaccinated
// (no VE_transmission knob in this version).

#include <Rcpp.h>
using namespace Rcpp;
#include <array>
#include <cmath>

// [[Rcpp::export]]
List model(List parscpp) {

  //natural history
  const double beta0( parscpp["beta"]);       //baseline transmission (seasonal multiplier applied per day)
  const double b1(   parscpp["b1"]);          //seasonal amplitude (Gaussian pulse, David rsvie)
  const double phi(  parscpp["phi"]);         //seasonal phase (fraction of year at peak)
  const double psi(  parscpp["psi"]);         //seasonal width (fraction of year)
  const double rEI(  parscpp["rEI"]);              //latency E -> I/U (single stage, no Erlang)
  const std::vector<double> rIR( parscpp["rIR"]);  //recovery I -> R/H rate, by age (David exposure-group avg)
  const std::vector<double> rUR( parscpp["rUR"]);  //recovery U -> R   rate, by age
  const double f(    parscpp["f"]);
  const double rH(   parscpp["rH"]);

  //vacc + waning
  const NumericVector vcov(   parscpp["vcov"]);    //coverage by age x IMD, length ng
  const NumericVector VE_inf( parscpp["VE_inf"]);  //efficacy against infection, length ng
  const NumericVector VE_sym( parscpp["VE_sym"]);  //efficacy against symptoms, length ng
  const NumericVector VE_hosp(parscpp["VE_hosp"]); //efficacy against hospitalisation, length ng
  const NumericVector VE_mort( parscpp["VE_mort"]);  //efficacy against mortality given hosp, length ng
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

  //demographic ageing (continuous): eta = ageing-out rate per band (eta[na-1] = death
  //rate from the top band); births = daily births into each IMD's youngest unvacc S.
  //All-zero eta & births disable ageing (default).
  const std::vector<double> eta(   parscpp["eta"]);     //length na
  const std::vector<double> births(parscpp["births"]);  //length ns

  const std::vector<double> cm( parscpp["cm"]);
  const int   cmdim1( parscpp["cmdim1"]);

  const std::vector<double> Sg0( parscpp["Sg0"]);
  const std::vector<double> Eg0( parscpp["Eg0"]);
  const std::vector<double> Ug0( parscpp["Ug0"]);
  const std::vector<double> Ig0( parscpp["Ig0"]);
  const std::vector<double> Rg0( parscpp["Rg0"]);
  const std::vector<double> Dg0( parscpp["Dg0"]);
  const std::vector<double> oNg( parscpp["oNg"]);

  //unvacc states
  std::vector<double> S_0(ng), E_0(ng);
  std::vector<double> I_0(ng), H_0(ng);
  std::vector<double> U_0(ng);
  std::vector<double> R_0(ng),  D_0(ng);
  //vacc states
  std::vector<double> V_0(ng);
  std::vector<double> Ev_0(ng);
  std::vector<double> Iv_0(ng), Hv_0(ng);
  std::vector<double> Uv_0(ng);
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
  NumericMatrix Nw_a(nd, na);  //living population by age band, per day (unvacc + vacc)
  //(age x IMD)-stratified hospitalisation incidence (rows=day, cols=ig=is*na+ia)
  NumericMatrix Hw_ag(nd, ng), Hwv_ag(nd, ng);
  //V-specific stratification (stock)
  NumericMatrix Vlev_a(nd, na), Vlev_s(nd, ns);

  NumericVector time(nt);
  int ig;

  for (int is = 0; is < ns; is++) {
  for (int ia = 0; ia < na; ia++) {
    ig = ia + is*na;
    //unvacc
    S_0[ig]  = Sg0[ig];
    E_0[ig]  = Eg0[ig];
    I_0[ig]  = Ig0[ig];
    U_0[ig]  = Ug0[ig];
    H_0[ig]  = 0.0;
    R_0[ig]  = Rg0[ig];
    D_0[ig]  = Dg0[ig];
    //vacc (all start empty)
    V_0[ig] = 0.0;
    Ev_0[ig] = 0.0;
    Iv_0[ig] = 0.0;
    Uv_0[ig] = 0.0;
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
  NumericVector Hpw_ag(ng), Hpw_vag(ng);  //daily H accumulators per (age x IMD)

  //temp working vars
  double yas, ua, ha, mHa, rrepa, cmi;
  double rVas, ve_i, ve_y, ve_h, ve_m;
  //unvacc temps
  double Sat, Eat, Uat, Iat, Hat, Rat, Dat;
  //vacc temps
  double Vat, Evat, Uvat, Ivat, Hvat, Rvat, Dvat;
  //flows
  double FOI, FOIS, FOIvV;
  double rEo, rUo, rIo, rHo;
  double rEvo, rUvo, rIvo, rHvo;
  double dVin, dEin, dIin, dUin, dHin;
  double dEvin, dIvin, dUvin, dHvin;
  double dS, dE, dI, dH, dU, dR, dD, dCc;
  double dV, dEv, dIv, dHv, dUv, dRv, dDv, dCcv, dVc;
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

    //snapshot for ageing (begin-of-step values -> exactly conserves people)
    std::vector<double> Sp=S_0, Ep=E_0, Ip=I_0, Up=U_0, Hp=H_0, Rp=R_0;
    std::vector<double> Vp=V_0, Evp=Ev_0, Ivp=Iv_0, Uvp=Uv_0, Hvp=Hv_0, Rvp=Rv_0;

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

      //FOI: sum across unvacc + vacc infectious (H/Hv excluded - isolated).
      //Reads the begin-of-step snapshot (Ip/Up/Ivp/Uvp) so the FOI is a Jacobi
      //step, independent of is/ia iteration order (I_0/U_0/... are updated in
      //place below).
      FOI = 0;
      for (int is2 = 0; is2 < ns; is2++) {
      for (int ia2 = 0; ia2 < na; ia2++) {
        ig2 = is2*na + ia2;
        icm = ig2*cmdim1 + ig;
        cmi = cm[icm];
        FOI += beta*ua*cmi*(
          Ip[ig2]  + f*Up[ig2] +
          Ivp[ig2] + f*Uvp[ig2]
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

      //daily incidence flows (for reporting)
      dEin  = dt*FOIS;
      dIin  = dt*yas*rEo;
      dUin  = dt*(1-yas)*rEo;
      dHin  = dt*ha*rIo;
      dVin  = dt*rVas*Sat;                              //new vaccinations
      dEvin = dt*FOIvV;
      dIvin = dt*yas*(1.0 - ve_y)*rEvo;
      dUvin = dt*(1.0 - yas*(1.0 - ve_y))*rEvo;
      dHvin = dt*ha*(1.0 - ve_h)*rIvo;

      //--- unvacc chain ---
      dS  = dt*( -FOIS - rVas*Sat + rW*Vat + rW_nat*Rat + rW_nat*Rvat );
      dE  = dt*(  FOIS - rEo );
      dI  = dt*(  yas*rEo - rIo );
      dU  = dt*(  (1-yas)*rEo - rUo );
      dH  = dt*(  ha*rIo - rHo );
      dR  = dt*(  (1-ha)*rIo + rUo + (1.0-mHa)*rHo - rW_nat*Rat );
      dD  = dt*(  mHa*rHo );
      dCc = dt*(  yas*rEo*rrepa );

      //--- vacc chain ---
      dV   = dt*( rVas*Sat - FOIvV - rW*Vat );
      dEv  = dt*(  FOIvV - rEvo );
      dIv  = dt*(  yas*(1.0 - ve_y)*rEvo - rIvo );
      dUv  = dt*(  (1.0 - yas*(1.0 - ve_y))*rEvo - rUvo );
      dHv  = dt*(  ha*(1.0 - ve_h)*rIvo - rHvo );
      dRv  = dt*(  (1.0 - ha*(1.0 - ve_h))*rIvo + rUvo + (1.0 - mHa*(1.0 - ve_m))*rHvo - rW_nat*Rvat );
      dDv  = dt*(  mHa*(1.0 - ve_m)*rHvo );
      dCcv = dt*(  yas*(1.0 - ve_y)*rEvo*rrepa );
      dVc  = dt*(  rVas*Sat );  //cumulative doses

      // --- demographic ageing: eta[ia] out (to ia+1; death from the top band), eta[ia-1] in ---
      dS -= dt*eta[ia]*Sat;  dE -= dt*eta[ia]*Eat;  dI -= dt*eta[ia]*Iat;
      dU -= dt*eta[ia]*Uat;  dH -= dt*eta[ia]*Hat;  dR -= dt*eta[ia]*Rat;
      dV -= dt*eta[ia]*Vat;  dEv-= dt*eta[ia]*Evat; dIv-= dt*eta[ia]*Ivat;
      dUv-= dt*eta[ia]*Uvat; dHv-= dt*eta[ia]*Hvat; dRv-= dt*eta[ia]*Rvat;
      if (ia > 0) {                 // ageing in from band below (same IMD); begin-of-step values
        double r = eta[ia-1];
        dS += dt*r*Sp[ig-1];  dE += dt*r*Ep[ig-1];  dI += dt*r*Ip[ig-1];
        dU += dt*r*Up[ig-1];  dH += dt*r*Hp[ig-1];  dR += dt*r*Rp[ig-1];
        dV += dt*r*Vp[ig-1];  dEv+= dt*r*Evp[ig-1]; dIv+= dt*r*Ivp[ig-1];
        dUv+= dt*r*Uvp[ig-1]; dHv+= dt*r*Hvp[ig-1]; dRv+= dt*r*Rvp[ig-1];
      } else {                      // youngest band: births enter unvaccinated susceptible
        dS += dt*births[is];
      }

      //--- update ---
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
      Hpw_ag[ig]  += dHin;
      Hpw_vag[ig] += dHvin;
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
        double n = 0;                                        // living pop (unvacc + vacc) by band
        for (int is = 0; is < ns; is++) { int g = is*na+ia;
          n += S_0[g]+E_0[g]+I_0[g]+U_0[g]+H_0[g]+R_0[g]
             + V_0[g]+Ev_0[g]+Iv_0[g]+Uv_0[g]+Hv_0[g]+Rv_0[g]; }
        Nw_a(day-1, ia) = n;
      }
      for (int ig2 = 0; ig2 < ng; ig2++) {
        Hw_ag(day-1, ig2)  = Hpw_ag[ig2];  Hpw_ag[ig2]  = 0;
        Hwv_ag(day-1, ig2) = Hpw_vag[ig2]; Hpw_vag[ig2] = 0;
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
    Named("Nw_a")  = Nw_a,     //living population by age band, per day
    Named("Hw_ag") = Hw_ag,    //daily H by (age x IMD), unvacc chain
    Named("Hwv_ag")= Hwv_ag,   //daily H by (age x IMD), vacc breakthrough chain
    Named("Vlev_a")= Vlev_a);

  return Rcpp::List::create(Rcpp::Named("byw") = byw, Rcpp::Named("byaw") = byaw);
}
