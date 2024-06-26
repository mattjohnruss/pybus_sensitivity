#include <R.h>
#include <stddef.h>
#include <math.h>

static double parameters[26];

#define k1 parameters[0]
#define k2 parameters[1]
#define k3 parameters[2]
#define k4 parameters[3]
#define k5 parameters[4]
#define k6 parameters[5]
#define k7 parameters[6]
#define k8 parameters[7]
#define k9 parameters[8]
#define k10 parameters[9]
#define k11 parameters[10]
#define k12 parameters[11]
#define k13 parameters[12]
#define k14 parameters[13]
#define k15 parameters[14]
#define k16 parameters[15]
#define k17 parameters[16]
#define k18 parameters[17]
#define k19 parameters[18]
#define k20 parameters[19]
#define k21 parameters[20]
#define k22 parameters[21]
#define ks parameters[22]
#define ku parameters[23]
#define k_ts parameters[24]
#define k_phi_a parameters[25]

// params with const values
static const double rhoa = 0.01;
static const double rhop = 1.0;
static const double rhoc = 1.0;
static const double rhom = 1.0;

// 17, 18, etc in days
// k_phia scaling from nondimensionalisation
static const size_t n_ti = 10;
static double ti_raw[] = {
    17.0,
    18.0,
    19.0,
    20.0,
    21.0,
    22.0,
    24.0,
    26.0,
    28.0,
    33.0
};

void initmod(void (* odeparms)(int *, double *)) {
    int n = 26;
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
        const double ti = (ti_raw[i] + k_ts) * k_phi_a;
        k_stim +=
            (ks / (ku * sqrt(2 * M_PI))) * exp(-0.5 * pow((time - ti) / ku, 2));
    }

    const double da_dt = k1 + (k_stim + k2 * a * c / (k3 + a * c)) * c * m * rhom / rhoa - k4 * ((rhop / rhoc) * p + c) * a - a;
    const double dp_dt = k5 * (1.0 - p / k6) * (1.0 + k7 * a * p / (k8 + a * p)) * p + (k9 * rhoc / rhop) * c - k10 * p;
    const double dc_dt = (k10 * rhoc / rhop) * p - k9 * c - (k11 * c - k12 * m / (k13 + m)) * c;
    const double dm_dt = (k14 + k15 * a * c / (k16 + a * c)) * rhoc * c / rhom + (k17 + k18 * a * p / (k19 + a * p)) * rhop * p / rhom + k20 * a / (k21 + a) - k22 * m;

    ydot[0] = da_dt;
    ydot[1] = dp_dt;
    ydot[2] = dc_dt;
    ydot[3] = dm_dt;
}
