#include <R.h>
#include <stddef.h>
#include <math.h>

static double parameters[10];

#define k1 parameters[0]
#define k2 parameters[1]
#define k3 parameters[2]
#define k4 parameters[3]
#define k5 parameters[4]
#define k10 parameters[5]
#define k11 parameters[6]
#define k17 parameters[7]
#define ks parameters[8]
#define ku parameters[9]

// define these globally so they can be used in `eqns` and the solver function
static const double k_ts = 100.0;
static const double k_phi_a = 10.0;

// params with const values
static const double k6 = 1.064;
static const double k7 = 17.3;
static const double k8 = 0.165;
static const double k9 = 0.063;
static const double k12 = 0.001;
static const double k13 = 1.0;
static const double k14 = 0.0001;
static const double k15 = 0.001;
static const double k16 = 1.0;
static const double k18 = 0.1;
static const double k19 = 10.0;
static const double k20 = 0.001;
static const double k21 = 1.0;
static const double k22 = 0.00545;
static const double rhoa = 0.01;
static const double rhop = 1.0;
static const double rhoc = 1.0;
static const double rhom = 1.0;

// 17, 18, etc in days
// k_phia scaling from nondimensionalisation
static const size_t n_ti = 10;
static double ti[] = {
    (17.0 + k_ts) * k_phi_a,
    (18.0 + k_ts) * k_phi_a,
    (19.0 + k_ts) * k_phi_a,
    (20.0 + k_ts) * k_phi_a,
    (21.0 + k_ts) * k_phi_a,
    (22.0 + k_ts) * k_phi_a,
    (24.0 + k_ts) * k_phi_a,
    (26.0 + k_ts) * k_phi_a,
    (28.0 + k_ts) * k_phi_a,
    (33.0 + k_ts) * k_phi_a
};

void initmod(void (* odeparms)(int *, double *)) {
    int n = 10;
    odeparms(&n, parameters);
}

void derivs(int *neq,
            double *t,
            double *y,
            double *ydot,
            double *yout,
            int *ip) {
    const double time = *t;

    const double a = y[0];
    const double p = y[1];
    const double c = y[2];
    const double m = y[3];

    double k_stim = 0.0;

    for (size_t i = 0; i < n_ti; ++i) {
      k_stim = k_stim +
        (ks / (ku * sqrt(2 * M_PI))) * exp(-0.5 * pow((time - ti[i]) / ku, 2));
    }

    const double da_dt = k1 + (k_stim + k2 * a * c / (k3 + a * c)) * c * m * rhom / rhoa - k4 * ((rhop / rhoc) * p + c) * a - a;
    const double dp_dt = k5 * (1 - p / k6) * (1 + k7 * a * p / (k8 + a * p)) * p + (k9 * rhoc / rhop) * c - k10 * p;
    const double dc_dt = (k10 * rhoc / rhop) * p - k9 * c - (k11 * c - k12 * m / (k13 + m)) * c;
    const double dm_dt = (k14 + k15 * a * c / (k16 + a * c)) * rhoc * c / rhom + (k17 + k18 * a * p / (k19 + a * p)) * rhop * p / rhom + k20 * a / (k21 + a) - k22 * m;

    ydot[0] = da_dt;
    ydot[1] = dp_dt;
    ydot[2] = dc_dt;
    ydot[3] = dm_dt;
}
