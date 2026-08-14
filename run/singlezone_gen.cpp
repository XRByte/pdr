// Single-zone Av=0 equilibrium generator, GOW_test.md sec 6/7 IC.
// Reproduces reference/Solution_chiXX_corrected.txt.
// args: chi kconst(0=Tdep default, 1=const 3e-17)
#include "gow17.h"
#include "cvodeDense.h"
#include "thermo.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv) {
  const double chi   = (argc>1)? atof(argv[1]) : 0.1;

  gow17 ode;
  const int dim = ode.Dimen();
  double y0[dim], abstol[dim];
  const double reltol = 1.0e-4, abstol0 = 1.0e-12;
  const int mxsteps = 5000000;
  const double temp = 100.0;
  const double xi_CR = 2e-16;
  const double tolfac = 10.0;
  const double tmax = 2000.0*1e6*3.16e7;
  const double tmin = 0.1*1e6*3.16e7;

  const double xC = 1.6e-4, xO = 3.2e-4, xHe = 0.1, xSi = 1.7e-6;

  const double nH_list[] = {10, 12.5, 15, 17.5, 20, 25, 30, 40, 50, 60, 70,
                            80, 90, 100, 125, 150, 175, 200, 250, 300, 400,
                            500, 600, 700, 800, 900};
  const int n_nH = sizeof(nH_list)/sizeof(nH_list[0]);

  for (int i=0;i<dim;i++){ abstol[i]=abstol0; y0[i]=0.0; }
  y0[ode.id("He+")] = 1.450654e-08;
  y0[ode.id("H3+")] = 2.681411e-07;
  y0[ode.id("C+")]  = 1.6e-4;
  y0[ode.id("H+")]  = 1.0e-4;
  y0[ode.id("S+")]  = 0.0;
  y0[ode.id("Si+")] = 0.0;
  y0[ode.id("CO")]  = 1.0e-7;
  y0[ode.id("H2")]  = 1.0e-10;
  y0[ode.id("E")]   = ode.GetE(temp, 1.0e-10, 3.0e-4);
  abstol[ode.id("He+")]=1e-15; abstol[ode.id("OHx")]=1e-15; abstol[ode.id("CHx")]=1e-15;
  abstol[ode.id("CO")]=1e-15; abstol[ode.id("C+")]=1e-15; abstol[ode.id("HCO+")]=1e-30;
  abstol[ode.id("H2")]=1e-8; abstol[ode.id("H+")]=1e-15; abstol[ode.id("H3+")]=1e-15;
  abstol[ode.id("H2+")]=1e-15; abstol[ode.id("E")]=ode.GetE(1.,0.1,0);

  CvodeDense cv(ode, reltol, abstol);
  cv.SetMxsteps(mxsteps);
  cv.SetMaxOrd(3);

  ode.SetGradv(3*3e-14);
  ode.SetNCOeffGlobal(false);   // Av=0 single zone: local gradv (CO cooling negligible)
  ode.SetbCOL(true);
  ode.SetIonRate(xi_CR);
  ode.SetZg(1.0); ode.SetZd(1.0);
  const double fgr = 1.0;
  ode.SetfHplusgr(fgr); ode.SetfHeplusgr(fgr); ode.SetfCplusgr(fgr);
  ode.SetfSplusgr(fgr); ode.SetfSiplusgr(fgr);
  ode.SetxCtot(xC); ode.SetxOtot(xO);
  ode.SetkgrH2(true);   // constant 3e-17 (match jfunc/athena)

  static double GPE[1], GISRF[1], Gph[7];
  GPE[0]=chi; GISRF[0]=chi; for(int i=0;i<7;i++) Gph[i]=chi;   // Av=0: no attenuation
  ode.SetRadField(GPE, Gph, GISRF);

  char fn[256];
  snprintf(fn, sizeof(fn), "../out_singlezone_gen/Solution_chi%g.txt", chi);
  FILE *pf = fopen(fn, "w");
  if (pf == NULL) {
    fprintf(stderr, "ERROR: cannot open %s for writing\n", fn);
    return 1;
  }
  fprintf(pf, "# nH T CO C C+ HCO+ OHx CHx O O+ 2H2 H+ H2+ H3+ He+ Si+ e\n");
  double y[dim];
  for (int inH=0; inH<n_nH; inH++) {
    const double nH = nH_list[inH];
    ode.SetnH(nH);
    ode.SetInit(0.0, y0);
    cv.ReInit();
    cv.SolveEq(tolfac, tmax, false, tmin);
    ode.CopyAbd(y);
    const double He_p=y[ode.id("He+")], OHx=y[ode.id("OHx")], CHx=y[ode.id("CHx")];
    const double CO=y[ode.id("CO")], Cp=y[ode.id("C+")], HCOp=y[ode.id("HCO+")];
    const double H2=y[ode.id("H2")], Hp=y[ode.id("H+")], H3p=y[ode.id("H3+")];
    const double H2p=y[ode.id("H2+")], Sp=y[ode.id("S+")], Sip=y[ode.id("Si+")];
    const double Op=y[ode.id("O+")], E=y[ode.id("E")];
    // ghost species by conservation
    const double C0 = xC - Cp - CO - HCOp - CHx;
    const double O0 = xO - Op - CO - HCOp - OHx;
    const double e  = Hp + Cp + He_p + Sip + Sp + H2p + H3p + HCOp + Op;
    const double T  = E / Thermo::CvCold(H2, xHe, e);
    fprintf(pf, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e\n",
      nH, T, CO, C0, Cp, HCOp, OHx, CHx, O0, Op, 2*H2, Hp, H2p, H3p, He_p, Sip, e);
    fflush(pf);
  }
  fclose(pf);
  printf("Wrote %s\n", fn);
  return 0;
}
