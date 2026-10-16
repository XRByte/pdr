// Single-density Av-sweep PDR slab (copy of av_sweep.cpp) with Andrea's
// initial abundances / abstols. Same grid, physics and CLI as av_sweep:
//   ./andrea_sweep [chi=1] [nH=1000] [field_geo=1]
// Solved to equilibrium
// per depth zone with self-consistent dust + H2/CO/C shielding via Slab.
// Av <-> NH: AV = NH*Zd/1.87e21  (radfield.cpp:41)
#include "gow17.h"
#include "cvodeDense.h"
#include "slab.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv) {
  gow17 ode;
  const double t0 = 0.0;
  const int dim = ode.Dimen();
  double y0[dim], abstol[dim];
  const double reltol = 1.0e-7, abstol0 = 1.0e-15;
  const int mxsteps = 5000000;
  const double temp = 100.0;
  const double xi_CR = 2e-16;
  const double tolfac = 1.5;
  const double tmax = 2000000.0*1e6*3.16e7;
  const double tmin = 0.1*1e6*3.16e7;

  const double xC = 1.6e-4, xO = 3.2e-4;

  // gas density (arg2, default 1000)
  const double nH = (argc > 2) ? atof(argv[2]) : 1000.0;
  // radiation field strength (one-sided slab -> G0*2 passed to Slab)
  const double chi = (argc > 1) ? atof(argv[1]) : 1.0;
  // Zdg = dust/gas metallicity
  const double Zdg = 1.0;

  // depth grid: log-spaced Av from 0.01 to 10 -> NH = 1.87e21 * Av (Zd=1)
  const long int ngrid = 1000;
  const double NH_min   = 1.0e17 / Zdg;          // Gong+2017 Fig16: N_min = 1e17/Z
  const double NH_total = 10.0 * 1.87e21 / Zdg;  // Av = 10   -> 1.87e22
  const bool islogNH = true;
  // radiation geometry (arg3, default 1): 0 = beamed (normal incidence),
  // 1 = Gong+2017 figs: isotropic approx, 60deg incidence (IsotropicApprox: 2x columns)
  const int field_geo = (argc > 3) ? atoi(argv[3]) : 1;
  const bool isdust = true, isfsH2 = true, isfsCO = true, isfsC = true;

  const bool verbose = false;

  // per-(nH,chi) output dir so a grid sweep does not clobber itself
  char dir[256];
  // beamed runs go to a separate tree so they do not overwrite the isotropic reference
  snprintf(dir, sizeof(dir), "../out_andrea_sweep%s/nH%g_chi%g/",
           field_geo == 0 ? "_beamed" : "", nH, chi);

  // initial abundances (Andrea's setup)
  for (int i=0;i<dim;i++){ abstol[i]=abstol0; y0[i]=0.0; }
  y0[ode.id("He+")] = 1.450654e-08;
  y0[ode.id("H3+")] = 2.681411e-07;
  y0[ode.id("C+")]  = 1.0e-4;
  y0[ode.id("S+")]  = 0.0;
  y0[ode.id("Si+")] = 0.0;
  y0[ode.id("CO")]  = 1.0e-7;
  y0[ode.id("H2")]  = 0.1;
  y0[ode.id("E")]   = ode.GetE(temp, 0.1, 0.0);
  abstol[ode.id("He+")]=1e-15; abstol[ode.id("OHx")]=1e-15; abstol[ode.id("CHx")]=1e-15;
  abstol[ode.id("CO")]=1e-15; abstol[ode.id("C+")]=1e-15; abstol[ode.id("HCO+")]=1e-30;
  abstol[ode.id("H2")]=1e-8; abstol[ode.id("H+")]=1e-15; abstol[ode.id("H3+")]=1e-15;
  abstol[ode.id("H2+")]=1e-15; abstol[ode.id("E")]=ode.GetE(1.,0.1,0);

  CvodeDense cv(ode, reltol, abstol);
  cv.SetMxsteps(mxsteps);
  cv.SetMaxOrd(5);

  ode.SetGradv(3*3e-14);
  ode.SetNCOeffGlobal(true);    // depth sweep: use global NCO for CO cooling
  ode.SetbCOL(true);
  ode.SetIonRate(xi_CR);
  ode.SetZg(Zdg); ode.SetZd(Zdg);
  const double fgr = 1.0;
  ode.SetfHplusgr(fgr); ode.SetfHeplusgr(fgr); ode.SetfCplusgr(fgr);
  ode.SetfSplusgr(fgr); ode.SetfSiplusgr(fgr);
  ode.SetxCtot(xC); ode.SetxOtot(xO);
  ode.SetkgrH2(true);           // constant 3e-17 grain H2 formation

  ode.SetInit(t0, y0);
  ode.SetnH(nH);

  // chemnet dump
  char fn[256];
  snprintf(fn, sizeof(fn), "%schemnet.dat", dir);
  FILE *pf_chemnet = fopen(fn, "w+");
  if (pf_chemnet) { ode.PrintChemNet(pf_chemnet); fclose(pf_chemnet); }

  Slab slab(ode, cv, ngrid, NH_total, chi*2.0, Zdg, islogNH, NH_min);
  slab.SetFieldGeo(field_geo);
  slab.IsH2MolSheilding(isfsH2);
  slab.IsCOMolSheilding(isfsCO);
  slab.IsCselfSheilding(isfsC);
  slab.IsDustSheilding(isdust);

  snprintf(fn, sizeof(fn), "%srates.dat", dir);
  FILE *pf_rates = fopen(fn, "w+");
  slab.SolveEq(tolfac, tmin, tmax, verbose, pf_rates);
  if (pf_rates) fclose(pf_rates);

  // abundances vs depth
  snprintf(fn, sizeof(fn), "%sslab.dat", dir);
  FILE *pf = fopen(fn, "w+");
  slab.WriteAbd(pf);
  fclose(pf);

  // NH column per zone (Av = NH/1.87e21)
  snprintf(fn, sizeof(fn), "%scolH_arr.dat", dir);
  FILE *pf_colH = fopen(fn, "w+");
  slab.WriteNH(pf_colH);
  fclose(pf_colH);

  // thermo rates
  snprintf(fn, sizeof(fn), "%sthermo.dat", dir);
  FILE *pf_thermo = fopen(fn, "w+");
  slab.WriteThermoRates(pf_thermo);
  fclose(pf_thermo);

  // profiling
  snprintf(fn, sizeof(fn), "%sprof.dat", dir);
  FILE *pf_prof = fopen(fn, "w+");
  slab.WriteProf(pf_prof);
  fclose(pf_prof);

  // parameters
  snprintf(fn, sizeof(fn), "%spara.dat", dir);
  FILE *pf_para = fopen(fn, "w+");
  if (pf_para) {
    ode.WriteParam(pf_para);
    fprintf(pf_para, "reltol = %0.4e\n", reltol);
    fprintf(pf_para, "tolfac = %0.4e\n", tolfac);
    fprintf(pf_para, "tmin = %0.4e\n", tmin);
    fprintf(pf_para, "nH = %0.4e\n", nH);
    fprintf(pf_para, "chi = %0.4e\n", chi);
    fprintf(pf_para, "ngrid = %ld\n", ngrid);
    fprintf(pf_para, "NH_min = %0.4e  (Av=%.3g)\n", NH_min, NH_min*Zdg/1.87e21);
    fprintf(pf_para, "NH_total = %0.4e  (Av=%.3g)\n", NH_total, NH_total*Zdg/1.87e21);
    fprintf(pf_para, "islogNH = %d\n", (int)islogNH);
    fclose(pf_para);
  }

  printf("Wrote %s*\n", dir);
  return 0;
}
