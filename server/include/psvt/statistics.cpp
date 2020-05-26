//**********************************************************************
//     A version of Uniformization Paleomagnetic Tool
//
//   Andrei V. Khokhlov khokhlov@ipgp.jussieu.fr
//   Copyright 2008 Andrei Khokhlov
//
//   This file is part of PaleoSecularVariation Toolkit.
//
//   PaleoSecularVariation Toolkit is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License,
//   or (at your option) any later version.
//
//   PaleoSecularVariation Toolkit is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//************************************************************************
#include "statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

int dbl_cmp(const void* first, const void* second) {
    if (*(reinterpret_cast<const double*>(first)) > *(reinterpret_cast<const double*>(second)))
        return 1;

    if (*(reinterpret_cast<const double*>(first)) < *(reinterpret_cast<const double*>(second)))
        return -1;

    return 0;
}

//    Anderson-Darling test for uniformity.   Given an ordered set
//              x_1<x_2<...<x_n
// of purported uniform [0,1) variates,  compute
//          a = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
// where z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, then find
//  v=adinf(a) and return  p=v+errfix(v), which should be uniform in [0,1),
//  that is, the p-value associated with the observed x_1<x_2<...<x_n.

// Short, practical version of full ADinf(z), z>0.
double adinf(double zzz) {
    if (zzz < 2.)
        return (std::exp(-1.2337141 / zzz) / std::sqrt(zzz) *
                (2.00012 + (.247105 - (.0649821 - (.0347962 - (.011672 - .00168691 * zzz) * zzz) * zzz) * zzz) * zzz));
    // max |error| < .000002 for zzz<2, (p=.90816...)
    else
        return (std::exp(-std::exp(
            1.0776 - (2.30695 - (.43424 - (.082433 - (.008056 - .0003146 * zzz) * zzz) * zzz) * zzz) * zzz)));
    // max |error|<.0000008 for 4<zzz<infinity
}

// The procedure  errfix(n,x)  corrects the error caused
// by using the asymptotic approximation, x=adinf(zzz).
// Thus x+errfix(n,x) is uniform in [0,1) for practical purposes;
// accuracy may be off at the 5th, rarely at the 4th, digit.
/*
double errfix(int NN, double xx)
{{
    double cf,tmp;
    if(xx>.8)
        return((-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*xx)*xx)*xx)*xx)*xx)/NN);
        cf=.01265+.1757/NN;
        if(xx<cf)
       {
         tmp=xx/cf;
             tmp=sqrt(tmp)*(1.-tmp)*(49*tmp-102);
             return(tmp*(.0037/(NN*NN)+.00078/NN+.00006)/NN);
       }
        tmp=(xx-cf)/(.8-cf);
        tmp=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*tmp)*tmp)*tmp)*tmp)*tmp;
        return(tmp*(.04213+.01365/NN)/NN);
}}*/

// The function AD(n,zzz) returns Prob(A_n<zzz) where
//    A_n = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
//          z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, and
//    x_1<x_2<...<x_n is an ordered set of iid uniform [0,1) variates.
//

double AD(int order, double z_arg) {
    double w_tmp[3];

    w_tmp[2] = adinf(z_arg);
    /* now w_tmp[2]=adinf(z_arg). Next, get w_tmp[1]=errfix(order,w_tmp[2]) and return w_tmp[2]+w_tmp[1]; */
    if (w_tmp[2] > .8) {
        w_tmp[1] =
            (-130.2137 +
             (745.2337 - (1705.091 - (1950.646 - (1116.360 - 255.7844 * w_tmp[2]) * w_tmp[2]) * w_tmp[2]) *
                             w_tmp[2]) *
                 w_tmp[2]) /
            order;
        return (w_tmp[2] + w_tmp[1]);
    }
    w_tmp[0] = .01265 + .1757 / order;
    if (w_tmp[2] < w_tmp[0]) {
        w_tmp[1] = w_tmp[2] / w_tmp[0];
        w_tmp[1] = std::sqrt(w_tmp[1]) * (1. - w_tmp[1]) * (49 * w_tmp[1] - 102);
        return (w_tmp[2] + w_tmp[1] * (.0037 / (order * order) + .00078 / order + .00006) / order);
    }
    w_tmp[1] = (w_tmp[2] - w_tmp[0]) / (.8 - w_tmp[0]);
    w_tmp[1] =
        -.00022633 +
        (6.54034 - (14.6538 - (14.458 - (8.259 - 1.91864 * w_tmp[1]) * w_tmp[1]) * w_tmp[1]) * w_tmp[1]) *
            w_tmp[1];
    return (w_tmp[2] + w_tmp[1] * (.04213 + .01365 / order) / order);
}

/// Marsaglia's C program  K(n, d) = Probability (Dn less  than d)

void mMultiply(double* A, double* B, double* C, int m) {
    int i, j, k;
    double s;
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++) {
            s = 0.;
            for (k = 0; k < m; k++)
                s += A[i * m + k] * B[k * m + j];
            C[i * m + j] = s;
        }
}

void mPower(double* A, int eA, double* V, int* eV, int m, int n) {
    double* B;
    int eB, i;
    if (n == 1) {
        for (i = 0; i < m * m; i++)
            V[i] = A[i];
        *eV = eA;
        return;
    }
    mPower(A, eA, V, eV, m, n / 2);
    // B=(double*)malloc((m*m)*sizeof(double));
    B = new double[m * m];
    mMultiply(V, V, B, m);
    eB = 2 * (*eV);
    if (n % 2 == 0) {
        for (i = 0; i < m * m; i++)
            V[i] = B[i];
        *eV = eB;
    } else {
        mMultiply(A, B, V, m);
        *eV = eA + eB;
    }
    if (V[(m / 2) * m + (m / 2)] > 1e140) {
        for (i = 0; i < m * m; i++)
            V[i] = V[i] * 1e-140;
        *eV += 140;
    }
    // free(B);
    delete[] B;
    return;
}
double KS_Estim(double d, int n) // K(int n,double d)
{
    int k, m, i, j, g, eH, eQ;
    double h, s, *H, *Q;

#ifndef ACCURACY_MORE_THAN_SEVEN_DIGITS
    /*OMIT IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL*/
    s = d * d * n;
    if (s > 7.24 || (s > 3.76 && n > 99))
        return (1 - 2 * std::exp(-(2.000071 + .331 / std::sqrt(1.0 * n) + 1.409 / n) * s));
#endif
    k = (int)(n * d) + 1;
    m = 2 * k - 1;
    h = k - n * d;
    // H=(double*)malloc((m*m)*sizeof(double));
    // Q=(double*)malloc((m*m)*sizeof(double));
    H = new double[2 * m * m];
    Q = &(H[m * m]);

    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
            if (i - j + 1 < 0)
                H[i * m + j] = 0;
            else
                H[i * m + j] = 1;
    for (i = 0; i < m; i++) {
        H[i * m] -= std::pow(h, i + 1);
        H[(m - 1) * m + i] -= std::pow(h, (m - i));
    }
    H[(m - 1) * m] += (2 * h - 1 > 0 ? std::pow(2 * h - 1, m) : 0);
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
            if (i - j + 1 > 0)
                for (g = 1; g <= i - j + 1; g++)
                    H[i * m + j] /= g;
    eH = 0;
    mPower(H, eH, Q, &eQ, m, n);
    s = Q[(k - 1) * m + k - 1];
    for (i = 1; i <= n; i++) {
        s = s * i / n;
        if (s < 1e-140) {
            s *= 1e140;
            eQ -= 140;
        }
    }
    s *= std::pow(10., eQ);
    // free(H);
    // free(Q);
    delete[] H;
    return (1.0 - s);
}
/// end Marsaglia's  C program

#define EPSILON1 0.001
#define EPSILON2 1.0e-8
#define RATHER_BIG 101

double Kolmogorov_D(double value) {
    int jj;
    double atwo, term, fac = 2., sum = 0., termbuf = 0.;

    atwo = -2. * value * value;
    for (jj = 1; jj < RATHER_BIG; jj++) {
        term = fac * std::exp(atwo * jj * jj);
        sum += term;
        if (std::fabs(term) <= EPSILON1 * termbuf || std::fabs(term) <= EPSILON2 * sum)
            return (sum);
        fac     = -fac;
        termbuf = std::fabs(term);
    }
    return (1.0);
}

double Kuiper_Q(double value) {
    if (value < 0.4)
        return (1.0);

    double atwo, sum = 0., termbuf = 0.;

    atwo = 2. * value * value;
    for (int jj = 1; jj < RATHER_BIG; jj++) {
        double fac  = (2.0 * atwo * jj * jj - 1);
        double term = fac * std::exp(-atwo * jj * jj);
        sum += term;
        if (std::fabs(term) <= EPSILON1 * termbuf || std::fabs(term) <= EPSILON2 * sum)
            return (2 * sum);

        termbuf = std::fabs(term);
    }
    return (2 * sum);
}

#undef EPSILON1
#undef EPSILON2
#undef RATHER_BIG

#define _ACCURACY 1.E-13

//  Fasano-Franceschini two-dimensional test
unsigned StTest2D(double& KS2D_measure, double& KS2D_estim, double& PearsonCorr, double* unifP, double* unifQ, unsigned NSample) {
    unsigned wtmp[5], maxloc(0);
    double MeanP(0.0), MeanQ(0.0), cov_sum[3], tmpval; // 0 for max,1,2,3,4 -- numbers in corresp quadrants
    const double InvN(1.0 / NSample), TINY(1.E-20);

    KS2D_measure = 0.0;
    memset(cov_sum, 0, 3 * sizeof(double));

    for (unsigned ii = 0; ii < NSample; ii++) {
        memset(wtmp, 0, 5 * sizeof(unsigned));

        for (unsigned jj = 0; jj < NSample; jj++) {
            if (ii < 1) {
                MeanP += unifP[jj];
                MeanQ += unifQ[jj];
            }

            if (unifQ[jj] > unifQ[ii]) {
                unifP[jj] > unifP[ii] ? ++wtmp[1] : ++wtmp[2];
            } else {
                unifP[jj] > unifP[ii] ? ++wtmp[4] : ++wtmp[3];
            }
        }
        for (unsigned kk = 1; kk < 5; kk++) {
            tmpval = fabs(wtmp[kk] * InvN - std::fabs(((((kk + 1) * (kk + 2)) / 2) % 2 - unifP[ii]) *
                                                      (((kk * (kk + 1)) / 2) % 2 - unifQ[ii])));

            if (tmpval > KS2D_measure) {
                maxloc       = ii;
                KS2D_measure = tmpval;
            }
            KS2D_measure = tmpval > KS2D_measure ? tmpval : KS2D_measure;
        }
        if (ii < 1) {
            MeanP *= InvN;
            MeanQ *= InvN;
        }
        cov_sum[0] += (unifP[ii] - MeanP) * (unifQ[ii] - MeanQ);
        cov_sum[1] += (unifP[ii] - MeanP) * (unifP[ii] - MeanP);
        cov_sum[2] += (unifQ[ii] - MeanQ) * (unifQ[ii] - MeanQ);
    }
    PearsonCorr = cov_sum[0] / (TINY + std::sqrt(cov_sum[1] * cov_sum[2]));
    KS2D_estim  = Kolmogorov_D(KS2D_measure * std::sqrt(static_cast<double>(NSample)) /
                              (1.0 + std::sqrt(1.0 - PearsonCorr * PearsonCorr) * (0.25 - 0.75 * std::sqrt(InvN))));
    return (maxloc);
}

void StTests1D(double* OrdUnif,
               double& max_dist,
               double& ks_estim,
               double& AD_measure,
               double& AD_estim, // reslt
               double* UnifSample,
               unsigned NSample) // args
{
    memcpy(OrdUnif, UnifSample, NSample * sizeof(double));
    qsort(OrdUnif, NSample, sizeof(double), dbl_cmp);
    // Data_Save(&TotDistPlot);
    AD_measure = 0.;

    // to be accurate with 0,1 margins for AD-test:
    if (OrdUnif[0] <= 0.0 || OrdUnif[NSample - 1] >= 1.0) // margin values are present-> shift or correction is needed!
    {
        double SmallShift = OrdUnif[0];
        SmallShift = 1.0 - OrdUnif[NSample - 1] > SmallShift ? 0.5 * (1.0 - OrdUnif[NSample - 1])
                                                             : -0.5 * SmallShift;
        if (SmallShift * SmallShift > 0.0) {
            for (unsigned ii = 0; ii < NSample; ii++)
                AD_measure += (2 * ii + 1) * std::log((OrdUnif[ii] + SmallShift) *
                                                 (1. - OrdUnif[NSample - 1 - ii] - SmallShift));
        } else {
            for (unsigned ii = 0; ii < NSample; ii++) {
                /* Float conversion precision might be unsufficient */
                double prot_val[2];
                prot_val[0] = (1. - OrdUnif[NSample - 1 - ii] > 0.) ? 1. - OrdUnif[NSample - 1 - ii] : _ACCURACY;
                prot_val[1] = OrdUnif[ii] > 0. ? OrdUnif[ii] : _ACCURACY;
                AD_measure += (2 * ii + 1) * std::log(prot_val[1] * prot_val[0]);
            }
        }
    } else {
        for (unsigned ii = 0; ii < NSample; ii++)
            AD_measure += (2 * ii + 1) * std::log((OrdUnif[ii]) * (1. - OrdUnif[NSample - 1 - ii]));
    }

    AD_measure = -(AD_measure / NSample + NSample);
    AD_estim   = 1.0 - AD(NSample, AD_measure);
    max_dist   = 0.;
    for (unsigned iii = 0; iii < NSample; iii++) {
        double dist_now = std::fabs(OrdUnif[iii] - iii / (1.0 * NSample));
        max_dist        = max_dist > dist_now ? max_dist : dist_now;
    }
    ks_estim = Kolmogorov_D(max_dist * (std::sqrt(1. * NSample) + 0.12 + 0.11 / std::sqrt(1. * NSample)));

    return;
}

void StTests1D(double* OrdUnif,
               double& max_dist,
               double& ks_estim,
               double& kuiper_v,
               double& kv_estim,
               double& AD_measure,
               double& AD_estim,
               double* UnifSample,
               unsigned NSample) // args
{
    memcpy(OrdUnif, UnifSample, NSample * sizeof(double));
    qsort(OrdUnif, NSample, sizeof(double), dbl_cmp);
    // Data_Save(&TotDistPlot);
    AD_measure = 0.;

    // to be accurate with 0,1 margins for AD-test:
    if (OrdUnif[0] <= 0.0 || OrdUnif[NSample - 1] >= 1.0) // margin values are present-> shift or correction is needed!
    {
        double SmallShift = OrdUnif[0];
        SmallShift = 1.0 - OrdUnif[NSample - 1] > SmallShift ? 0.5 * (1.0 - OrdUnif[NSample - 1])
                                                             : -0.5 * SmallShift;
        if (SmallShift * SmallShift > 0.0) {
            for (unsigned ii = 0; ii < NSample; ii++)
                AD_measure += (2 * ii + 1) * std::log((OrdUnif[ii] + SmallShift) *
                                                 (1. - OrdUnif[NSample - 1 - ii] - SmallShift));
        } else {
            for (unsigned ii = 0; ii < NSample; ii++) {
                /* Float conversion precision might be unsufficient */
                double prot_val[2];
                prot_val[0] = (1. - OrdUnif[NSample - 1 - ii] > 0.) ? 1. - OrdUnif[NSample - 1 - ii] : _ACCURACY;
                prot_val[1] = OrdUnif[ii] > 0. ? OrdUnif[ii] : _ACCURACY;
                AD_measure += (2 * ii + 1) * std::log(prot_val[1] * prot_val[0]);
            }
        }
    } else {
        for (unsigned ii = 0; ii < NSample; ii++)
            AD_measure += (2 * ii + 1) * std::log((OrdUnif[ii]) * (1. - OrdUnif[NSample - 1 - ii]));
    }

    AD_measure = -(AD_measure / NSample + NSample);
    AD_estim   = 1.0 - AD(NSample, AD_measure);
    max_dist   = 0.;
    for (unsigned iii = 0; iii < NSample; iii++) {
        double dist_now = std::fabs(OrdUnif[iii] - iii / (1.0 * NSample));
        max_dist        = max_dist > dist_now ? max_dist : dist_now;
    }
    ks_estim = Kolmogorov_D(max_dist * (std::sqrt(1. * NSample) + 0.12 + 0.11 / std::sqrt(1. * NSample)));

    double max_dist_plus(0), max_dist_minus(0);
    for (unsigned iii = 0; iii < NSample; iii++) {
        double dist_plus  = (OrdUnif[iii] - iii / (1.0 * NSample)),
               dist_minus = (iii / (1.0 * NSample) - OrdUnif[iii]);
        max_dist_plus     = max_dist_plus > dist_plus ? max_dist_plus : dist_plus;
        max_dist_minus    = max_dist_minus > dist_minus ? max_dist_minus : dist_minus;
    }
    kuiper_v = max_dist_minus + max_dist_plus;

    kv_estim = Kuiper_Q(kuiper_v * (std::sqrt(1. * NSample) + 0.155 + 0.24 / std::sqrt(1. * NSample)));

    return;
}

#undef _ACCURACY

const double kAccuracy = 1.E-13;

UnifTestResults UniformityTests(const std::vector<double>& sample) {
    std::vector<double> sample_sorted(sample.size());
    std::partial_sort_copy(sample.begin(), sample.end(), sample_sorted.begin(), sample_sorted.end());

    UnifTestResults results;
    double AD_measure = 0;

    // we must be accurate with 0 and 1 margins for AD-test:
    if (sample_sorted[0] <= 0. || sample_sorted[sample_sorted.size() - 1] >= 1.) {
        // margin values are present, thus shift or correction is needed

        double shift = sample_sorted[0];
        if (1. - sample_sorted[sample_sorted.size() - 1] > shift)
            shift = 0.5 * (1. - sample_sorted[sample_sorted.size() - 1]);
        else
            shift *= -0.5;

        if (shift * shift > 0.) {
            for (int i = 0; i < sample_sorted.size(); ++i)
                AD_measure += (2 * i + 1) * std::log((sample_sorted[i] + shift) * (1. - sample_sorted[sample_sorted.size() - 1 - i] - shift));
        } else {
            for (int i = 0; i < sample_sorted.size(); ++i) {
                // float conversion precision might be insufficient
                double prot_val[2];
                prot_val[0] = (1. - sample_sorted[sample_sorted.size() - 1 - i] > 0.) ? 1. - sample_sorted[sample_sorted.size() - 1 - i] : kAccuracy;
                prot_val[1] = sample_sorted[i] > 0. ? sample_sorted[i] : kAccuracy;
                AD_measure += (2 * i + 1) * std::log(prot_val[1] * prot_val[0]);
            }
        }
    } else {
        for (int i = 0; i < sample_sorted.size(); ++i)
            AD_measure += (2 * i + 1) * std::log(sample_sorted[i] * (1. - sample_sorted[sample_sorted.size() - 1 - i]));
    }
    
    results.AD_measure = -(AD_measure / sample_sorted.size() + sample_sorted.size());
    results.AD_estimate = 1. - AD(sample_sorted.size(), results.AD_measure);

    double max_distance = 0;
    for (int i = 0; i < sample_sorted.size(); ++i) {
        double distance = std::fabs(sample_sorted[i] - i / (1. * sample_sorted.size()));
        max_distance    = max_distance > distance ? max_distance : distance;
    }

    results.KS_measure = max_distance;
    results.KS_estimate = Kolmogorov_D(max_distance * (std::sqrt(1. * sample_sorted.size()) + 0.12 + 0.11 / std::sqrt(1. * sample_sorted.size())));

    return results;
}
