#include "legendre_rule.h"
#include "matrix_elements.h"
#include "Spectrum.h"
#include <chrono>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <time.h>
#include <vector>
using namespace std;

struct Settings {
  double me = 0.000511;
  double m2 = 0.001;
  double Gf = 0.0000116637;
  double g = 1.0;
  double m = 0.105660;
  int done = 0;
  string filename = "	";
  double (*pMatrixElement)(double E1, double E2, double E3, double ct1,
                           double ct2, double ct3, double ph2, double ph3,
                           double Gf, double m, double me, double m2, double g);
  bool operator==(const Settings &a) const {
    return (m2 == a.m2 && done == a.done);
  }
} sett;


// Calculate derivative appearing in the phasespace integral
double Deriv(double E1, double E2, double ct1, double ct2, double ct3,
             double ph2, double ph3, double m, double me, double m2) {
  return 2 * (E1 + E2 - m -
              (ct2 * ct3 + cos(ph2 - ph3) * pow(1 - pow(ct2, 2), 0.5) *
                               pow(1 - pow(ct3, 2), 0.5)) *
                  pow(pow(E2, 2) - pow(m2, 2), 0.5) -
              pow(pow(E1, 2) - pow(me, 2), 0.5) *
                  (ct1 * ct3 + pow(1 - pow(ct1, 2), 0.5) *
                                   pow(1 - pow(ct3, 2), 0.5) * sin(ph3)));
}

// Solve kinetic equation for E3
double CalcE3(double E1, double E2, double ct1, double ct2, double ct3,
              double ph2, double ph3, double m, double me, double m2) {
  return -(pow(E1 + E2 - m -
                   (ct2 * ct3 + cos(ph2 - ph3) * pow(1 - pow(ct2, 2), 0.5) *
                                    pow(1 - pow(ct3, 2), 0.5)) *
                       pow(pow(E2, 2) - pow(m2, 2), 0.5) -
                   pow(pow(E1, 2) - pow(me, 2), 0.5) *
                       (ct1 * ct3 + pow(1 - pow(ct1, 2), 0.5) *
                                        pow(1 - pow(ct3, 2), 0.5) * sin(ph3)),
               -1) *
           (2 * E1 * E2 - 2 * E1 * m - 2 * E2 * m + pow(m, 2) + pow(m2, 2) +
            pow(me, 2) -
            2 * pow(pow(E2, 2) - pow(m2, 2), 0.5) *
                pow(pow(E1, 2) - pow(me, 2), 0.5) *
                (ct1 * ct2 + pow(1 - pow(ct1, 2), 0.5) *
                                 pow(1 - pow(ct2, 2), 0.5) * sin(ph2)))) /
         2.;
}

// Implements a test for physicality, thus avoiding complicated region
int isPhysical(double E1, double E2, double E3, double m, double me,
               double m2) {
  if (E1 < me || E2 < 0 || E3 < 0 || E1 + E2 + E3 > m)
    return 0;
  return 1;
}

//Integrand of the spectrum in the form for GSL MonteCarlo integration
double MCIntegrand_Spectrum(double *k, size_t dim, void *params) {
  double res = 0.0;
  double ct2 = k[0];
  double ct3 = k[1];
  double ph2 = k[2];
  double ph3 = k[3];
  double E2 = k[4];
  double *par = (double *)params;
  double ct1 = par[0];
  double E1 = par[1];
  double Gf = par[2];
  double m = par[3];
  double me = par[4];
  double m2 = par[5];
  double g = par[6];

  double E3 = CalcE3(E1, E2, ct1, ct2, ct3, ph2, ph3, m, me, m2);
  if (isPhysical(E1, E2, E3, m, me, m2)) {
    double weight = pow(2.0 * M_PI, -7) / 8.0 * sqrt(E1 * E1 - me * me) *
                    sqrt(E2 * E2 - m2 * m2) * E3 /
                    abs(Deriv(E1, E2, ct1, ct2, ct3, ph2, ph3, m, me, m2));
    res = weight * sett.pMatrixElement(E1, E2, E3, ct1, ct2, ct3, ph2, ph3, Gf,
                                       m, me, m2, g);
  }
  return res;
}

//Integrand of the total width in the form for GSL MonteCarlo integration
double MCIntegrand_TotalWidth(double *k, size_t dim, void *params) {
  double res = 0.0;
  double ct1 = k[0];
  double ct2 = k[1];
  double ct3 = k[2];
  double ph2 = k[3];
  double ph3 = k[4];
  double E1 = k[5];
  double E2 = k[6];
  double *par = (double *)params;
  double Gf = par[0];
  double m = par[1];
  double me = par[2];
  double m2 = par[3];
  double g = par[4];

  double E3 = CalcE3(E1, E2, ct1, ct2, ct3, ph2, ph3, m, me, m2);
  if (isPhysical(E1, E2, E3, m, me, m2)) {
    double weight = pow(2.0 * M_PI, -7) / 8.0 * sqrt(E1 * E1 - me * me) *
                    sqrt(E2 * E2 - m2 * m2) * E3 /
                    abs(Deriv(E1, E2, ct1, ct2, ct3, ph2, ph3, m, me, m2));
    res = weight * sett.pMatrixElement(E1, E2, E3, ct1, ct2, ct3, ph2, ph3, Gf,
                                       m, me, m2, g);
  }
  return res;
}

//Calculate BSM spectrum using GSL-MC-integration at E1 and cost1
double BSM_Spectrum(double E1, double cost1) {
  double res, err;
  double xl[5] = {-1.0, -1.0, 0.0, 0.0, sett.m2};
  double xu[5] = {1.0, 1.0, 2.0 * M_PI, 2.0 * M_PI, (sett.m + sett.m2) / 2.0};

    //Prepare parameters
  const gsl_rng_type *T;
  gsl_rng *r;
  size_t calls = 100000;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_monte_function G;
  G.f = &MCIntegrand_Spectrum;
  G.dim = 5;

  double par[7];
  par[0] = cost1;
  par[1] = E1;
  par[2] = sett.Gf;
  par[3] = sett.m;
  par[4] = sett.me;
  par[5] = sett.m2;
  par[6] = sett.g;
  G.params = par;

  //Setup integrator
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);
  gsl_monte_vegas_integrate(&G, xl, xu, 5, calls / 5, r, s, &res, &err);
  int i = 0;
  int iMax = 10;
  //Run integration at most 10 rounds
  do {
    i++;
    gsl_monte_vegas_integrate(&G, xl, xu, 5, calls / 5, r, s, &res, &err);
  } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 && i <= iMax);
  gsl_monte_vegas_free(s);
  return res / (2.0 * sett.m);
}

//Generate BSM-Spectrum with nodes given by xcost and xE1 for use in quadrature integration
vector<vector<double>> GenerateBSMSpectrum(int orderE1, int ordercost1,
                                           double *xcost1, double *xE1) {
  vector<vector<double>> spectrum(orderE1, vector<double>(ordercost1));
  int N = 0;
  omp_set_dynamic(0);
  omp_set_num_threads(8);
#pragma omp parallel for shared(N, spectrum)
  for (int i = 0; i < orderE1; i++) {
    for (int j = 0; j < ordercost1; j++) {
      spectrum[i][j] = BSM_Spectrum(xE1[i], xcost1[j]);
#pragma omp critical
      {
        printf("%.1f %% %s\n", 100.0 * ((double)(++N) / (orderE1 * ordercost1)),
               sett.filename.c_str());
      }
    }
  }
  printf("Done\n");

  return spectrum;
}

//Convert the given Spectrum to its total width
double GetTotalWidth(int orderE1, int ordercost1, double *xcost1,
                     double *wcost1, double *xE1, double *wE1,
                     vector<vector<double>> spectrum) {
  double res = 0.0;
  for (int i = 0; i < orderE1; i++) {
    for (int j = 0; j < ordercost1; j++) {
      res += spectrum[i][j] * wE1[i] * wcost1[j];
    }
  }
  return res;
}

//Convert raw BSM spectrum to clean object
Spectrum GenerateBSM(int orderE1, int ordercost1) {
  vector<double> xE1;
  vector<double> xcost1;
  vector<vector<double>> spectrum;
  double totalBSMWidth;

  double *_wcost1 = new double[ordercost1];
  double *_xcost1 = new double[ordercost1];
  double *_wE1 = new double[orderE1];
  double *_xE1 = new double[orderE1];
  cgqf(ordercost1, 1, 0, 0, -1.0, 1.0, _xcost1, _wcost1);
  cgqf(orderE1, 1, 0, 0, sett.me, sett.m / 2, _xE1, _wE1);
  spectrum = GenerateBSMSpectrum(orderE1, ordercost1, _xcost1, _xE1);
  totalBSMWidth = GetTotalWidth(orderE1, ordercost1, _xcost1, _wcost1, _xE1,
                                _wE1, spectrum);

  xE1.resize(orderE1);
  xcost1.resize(ordercost1);
  for (int i = 0; i < orderE1; i++) {
    xE1[i] = _xE1[i];
  }
  for (int i = 0; i < ordercost1; i++) {
    xcost1[i] = _xcost1[i];
  }
  delete[] _wcost1;
  delete[] _xcost1;
  delete[] _xE1;
  delete[] _wE1;
  Spectrum spec(orderE1, ordercost1, spectrum, xE1, xcost1, totalBSMWidth);
  return spec;
}

//Gives the known standard model description of the muon decay spectrum
double diffrate(double *x, double *param) {
  /* Param
   * 0 : Scale
   * 1 : rho
   * 2 : eta
   * 3 : xi
   * 4 : delta
   *
   * /*/
  double mmu = sett.m;
  double x0 = sett.me * 2.0 * mmu / (sett.me * sett.me + mmu * mmu);
  double pmu = -1.0; // decay to e-
  return x[0] * x[0] * param[0] *
         (3.0 * (1.0 - x[0]) + 2.0 * param[1] / 3.0 * (4.0 * x[0] - 3.0) +
          3.0 * param[2] * x0 * (1 - x[0]) / x[0] +
          pmu * param[3] * x[1] *
              (1 - x[0] + 2.0 * param[4] / 3.0 * (4.0 * x[0] - 3.0)));
}

//Calculates the squared difference in order to find the best fit to the
//combined spectrum by minimisation
double SquaredDiff(const gsl_vector *v, void *params) {
  Spectrum *spec = static_cast<Spectrum *>(params);
  double SM_Width = sett.Gf * sett.Gf * pow(sett.m, 5) / 192.0 / pow(M_PI, 3);
  double param[] = {spec->totalWidth / sett.m * 4, gsl_vector_get(v, 0), 0.0,
                    gsl_vector_get(v, 1), gsl_vector_get(v, 2)};

  double res = 0.0;
  double x[2];
  for (size_t i = 0; i < spec->orderX; i++) {
    for (size_t j = 0; j < spec->orderY; j++) {
      x[0] = spec->X_values[i];
      x[1] = spec->Y_values[j];
      res += pow(spec->d_spectrum[i][j] - diffrate(x, param), 2.0);
    }
  }
  return res;
}

//Extract best-fit-parameters by minimising squared differences
int ExtraxtParams(Spectrum spec, double params[3]) {
  //Setup simplex-minimiser
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;
  //Set initial fit-parameters to 1
  x = gsl_vector_alloc(3);
  gsl_vector_set(x, 0, 1);
  gsl_vector_set(x, 1, 1);
  gsl_vector_set(x, 2, 1);
  ss = gsl_vector_alloc(3);
  gsl_vector_set_all(ss, 1.0);

  //Defince function to be fittet with parameters
  minex_func.n = 3;
  minex_func.f = SquaredDiff;
  minex_func.params = &spec;

  s = gsl_multimin_fminimizer_alloc(T, 3);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  //Run minimiser
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status)
      break;
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-8);
    if (status == GSL_SUCCESS) {
    }

  } while (status == GSL_CONTINUE && iter < 1000);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  //Extract best-fit-parameters
  params[0] = gsl_vector_get(s->x, 0);
  params[1] = gsl_vector_get(s->x, 1);
  params[2] = gsl_vector_get(s->x, 2);
  gsl_multimin_fminimizer_free(s);

  return status;
}

//Compare best-fit parameters with experimental results up to simgma confidence limit
int IsIncompatible(double parm[3], double sigma) {
  double rho = 0.74979;
  double rho_err = 0.00026;
  double xi = 1.0009;
  double xi_err_m = 0.0007;
  double xi_err_p = 0.0016;
  double delta = 0.75047;
  double delta_err = 0.00034;
  if (parm[0] < rho - sigma * rho_err || parm[0] > rho + sigma * rho_err)
    return 1;
  if (parm[1] < xi - sigma * xi_err_m || parm[1] > xi + sigma * xi_err_p)
    return 2;
  if (parm[2] < delta - sigma * delta_err ||
      parm[2] > delta + sigma * delta_err)
    return 3;
  return 0;
}

//Find coupling constant that saturates experimental results
double FitToSpectrum(Spectrum bsmSpec, double sigma) {
  vector<double> xX(bsmSpec.orderX);
  //Prepare pure standard model spectrum
  double rho = 0.75;
  double delta = 0.75;
  double xi = 1.0;
  double SM_Width = sett.Gf * sett.Gf * pow(sett.m, 5) / 192.0 / pow(M_PI, 3);
  double totalRate;
  double param[] = {SM_Width / sett.m * 4, rho, 0.0, xi, delta};
  vector<vector<double>> CombinedSpectrum(bsmSpec.orderX,
                                          vector<double>(bsmSpec.orderY));

  double par_extracted[3];
  double x[2];
  sett.g = 21;
  double dg = 1.0;
  int comp = 0;
  //Find correct coupling constant
  for (size_t i = 0; i < 1000; i++) {
    if ((int)i == 99) {
      sett.g = 999.0;
      break;
    }
    if (dg < 1e-10)
      break;
    sett.g -= dg;
    for (int i = 0; i < bsmSpec.orderX; i++) {
      for (int j = 0; j < bsmSpec.orderY; j++) {
        x[0] = bsmSpec.X_values[i] * 2.0 * sett.m /
               (sett.me * sett.me + sett.m * sett.m);
        x[1] = bsmSpec.Y_values[j];
        xX[i] = x[0];
        CombinedSpectrum[i][j] =
            bsmSpec.d_spectrum[i][j] * sett.g * sett.g + diffrate(x, param);
      }
    }

    totalRate = SM_Width + bsmSpec.totalWidth * sett.g * sett.g;
    Spectrum spec(bsmSpec.orderX, bsmSpec.orderY, CombinedSpectrum, xX,
                  bsmSpec.Y_values, totalRate);

    ExtraxtParams(spec, par_extracted);
    comp = IsIncompatible(par_extracted, sigma);
    if (comp == 0) {
      sett.g += dg;
      dg /= 10.0;
    }
  }
  return sett.g;
}

//Derive coupling constant for a random mass
void RunRandomScan() {
  double lower_bound = -4.0;
  double upper_bound = log10(sett.m);
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::default_random_engine re(
      std::chrono::system_clock::now().time_since_epoch().count());
  double randmass = unif(re);
  printf("%E\n", pow(10.0, randmass));
  sett.g = 1.0;
  sett.m2 = pow(10.0, randmass);
  Spectrum bsmSpec = GenerateBSM(25, 25);
  //Consider 90% and 95% confidence intervals
  double res90 = FitToSpectrum(bsmSpec, 1.64485);
  double res95 = FitToSpectrum(bsmSpec, 1.95996);
  FILE *pFile;
  pFile = fopen(sett.filename.c_str(), "a");
  fprintf(pFile, "%E,%E,%E\n", sett.m2, res90, res95);
  fclose(pFile);
  printf("%E,%E,%E\n", sett.m2, res90, res95);
}

void GetScalarElectron() {
  sett.filename = "results/ScalarElectron.csv";
  sett.pMatrixElement = M2ScalarElectronCoupling;
  RunRandomScan();
}

void GetScalarMuon() {
  sett.filename = "results/ScalarMuon.csv";
  sett.pMatrixElement = M2ScalarMuonCoupling;
  RunRandomScan();
}

void GetScalarLepton() {
  sett.filename = "results/ScalarLepton.csv";
  sett.pMatrixElement = M2ScalarToLepton;
  RunRandomScan();
}

void GetScalarYukawa() {
  sett.filename = "results/ScalarYukawa.csv";
  sett.pMatrixElement = M2ScalarYukawa;
  RunRandomScan();
}

void GetVectorElectron() {
  sett.filename = "results/VectorElectron.csv";
  sett.pMatrixElement = M2VectorElectronRLCoupling;
  RunRandomScan();
}

void GetVectorMuon() {
  sett.filename = "results/VectorMuon.csv";
  sett.pMatrixElement = M2VectorMuonRLCoupling;
  RunRandomScan();
}

void ScanAll() {
  for (size_t i = 0; i < 1000; i++) {
    GetScalarElectron();
    GetScalarMuon();
    GetScalarLepton();
    GetScalarYukawa();
    GetVectorElectron();
    GetVectorMuon();
  }
}

void testrun() {
  sett.filename = "results/dump.csv";
  sett.pMatrixElement = M2VectorMuonRLCoupling;
  double res;
  double randmass = -3.0;
  printf("%E\n", pow(10.0, randmass));
  sett.g = 1.0;
  sett.m2 = pow(10.0, randmass);
  Spectrum bsmSpec = GenerateBSM(10, 10);
  res = FitToSpectrum(bsmSpec, 1.95996);
  printf("%E,%E,%E\n", sett.m2, FitToSpectrum(bsmSpec, 1.64485),
         FitToSpectrum(bsmSpec, 1.95996));
}

int main(int argc, char *argv[]) {
/*
real	2m18.210s
user	11m10.440s
sys	0m0.487s
  testrun();*/
  ScanAll();
  return 0;
}
