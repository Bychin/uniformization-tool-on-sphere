/**
 * @file FST_semi_memo.c
 * @brief Routines to perform convolutions on the 2-sphere using a combination of seminaive and
 * naive algorithms.
 *
 * Assumes that all precomputed data is already in memory (e.g. tables are generated).\n
 * For descriptions on calling these functions, see the documentation preceding each function.
 */

#include "s2kit/FST_semi_memo.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "s2kit/cospml.h"
#include "s2kit/naive.h"
#include "s2kit/seminaive.h"
#include "s2kit/weights.h"
#include "s2kit/chebyshev_nodes.h"
#include "s2kit/util.h"

/**
 * @brief Computes <b>spherical harmonic transform</b> using the seminaive and naive algorithms.
 *
 * Output ordering of coeffs (in @p rcoeffs and @p icoeffs) <tt>f(m,l)</tt> is:
 * @code
 * f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
 *        f(1,1) f(1,2) ... f(1,bw-1)
 *        etc.
 *            f(bw-2,bw-2), f(bw-2,bw-1)
 *                          f(bw-1,bw-1)
 *                          f(-(bw-1),bw-1)
 *         f(-(bw-2),bw-2), f(-(bw-2),bw-1)
 *        etc.
 *              f(-2,2) ... f(-2,bw-1)
 *      f(-1,1) f(-1,2) ... f(-1,bw-1)
 * @endcode
 *
 * This only requires an array of size @c bw*bw. If zero-padding is used to make the indexing nice,
 * then you need a an <tt>(2bw-1)*bw</tt> array, but that is not done here. Because of the amount
 * of space necessary for doing large transforms, it is important not to use any more than
 * necessary.
 *
 * @note See an example of use in test_s2_semi_memo.c and test_s2_semi_memo_fwd.c
 *
 * @param rdata array of length <tt>4*bw*bw</tt> of the real part of the function samples
 * @param idata array of length <tt>4*bw*bw</tt> of the imaginary part of the function samples
 * @param rcoeffs array of length <tt>bw*bw</tt> which will contain the real part of harmonic
 * coefficients in a linearized form
 * @param icoeffs array of length <tt>bw*bw</tt> which will contain the imaginary part of harmonic
 * coefficients in a linearized form
 * @param bw bandwidth of problem
 * @param seminaive_naive_table pre-generated spharmonic Pml table
 * @param workspace space for computations of size <tt>(8*bw^2)+(7*bw)</tt>
 * @param data_format determines the format of the computed data
 * @param cutoff determines order to switch from seminaive to naive algorithm
 * @param DCT_plan plan for DLT which is used as argument for call DLTSemi()
 * @param FFT_plan plan for FFT
 * @param weights array which is used as argument for call DLTSemi() and DLTNaive()
 *
 * @note @p seminaive_naive_table should be generated by Spharmonic_Pml_Table(). This table can be
 * re-used in the inverse transform and, for example, in series of convolutions.
 * @note See more info about @p DCT_plan and @p weights in DLTSemi().
 */
void FSTSemiMemo(double* rdata, double* idata, double* rcoeffs, double* icoeffs, const int bw,
                 double** seminaive_naive_table, double* workspace, DataFormat data_format,
                 const int cutoff, fftw_plan* DCT_plan, fftw_plan* FFT_plan, double* weights) {
    int size = 2 * bw;

    // total workspace is (8 * bw^2) + (7 * bw)
    double* rres = workspace;               // needs (size * size) = (4 * bw^2)
    double* ires = rres + (size * size);    // needs (size * size) = (4 * bw^2)
    double* fltres = ires + (size * size);  // needs bw
    double* eval_pts = fltres + bw;         // needs (2 * bw)
    double* scratchpad = eval_pts + (size); // needs (4 * bw)

    // do the FFTs along phi
    fftw_execute_split_dft(*FFT_plan, rdata, idata, rres, ires);

    // Normalize
    // The associated Legendres are of norm 1, so because the spherical
    // harmonics are of norm 1 we are using coeff of sqrt(2*pi) here.
    double normed_coeff = sqrt(2. * M_PI) / size;
    for (int i = 0; i < size * size; ++i) {
        rres[i] *= normed_coeff;
        ires[i] *= normed_coeff;
    }

    // point to start of output data buffers
    double* rdataptr = rcoeffs;
    double* idataptr = icoeffs;

    for (int m = 0; m < cutoff; ++m) { // semi-naive part
        // real part
        DLTSemi(rres + (m * size), bw, m, fltres, scratchpad, seminaive_naive_table[m], weights, DCT_plan);
        // load real part of coefficients into output space
        memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
        rdataptr += bw - m;

        // imaginary part
        DLTSemi(ires + (m * size), bw, m, fltres, scratchpad, seminaive_naive_table[m], weights, DCT_plan);
        // load imaginary part of coefficients into output space
        memcpy(idataptr, fltres, sizeof(double) * (bw - m));
        idataptr += bw - m;
    }

    for (int m = cutoff; m < bw; ++m) { // naive part
        // real part
        DLTNaive(rres + (m * size), bw, m, weights, fltres, seminaive_naive_table[m], scratchpad);
        memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
        rdataptr += bw - m;

        // imaginary part
        DLTNaive(ires + (m * size), bw, m, weights, fltres, seminaive_naive_table[m], scratchpad);
        memcpy(idataptr, fltres, sizeof(double) * (bw - m));
        idataptr += bw - m;
    }

    // upper coefficients
    /*
        If the data is real, we don't have to compute the coeffs whose
        order is less than 0, since the data is real, we know that

        f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m)),

        so use that to get the rest of the coefficients
    */
    if (data_format == REAL) {
        double coeff = 1.;
        for (int i = 1; i < bw; ++i) {
            coeff *= -1.;
            for (int j = i; j < bw; ++j) {
                int index0 = IndexOfHarmonicCoeff(i, j, bw);
                int index1 = IndexOfHarmonicCoeff(-i, j, bw);

                rcoeffs[index1] = coeff * rcoeffs[index0];
                icoeffs[index1] = -coeff * icoeffs[index0];
            }
        }

        return;
    }
    
    // complex data
    /* 
        Note that `m` is greater than `bw` here, but this is for
        purposes of indexing the input data arrays. The "true"
        value of `m` as a parameter for Pml is `size-m` 
    */
    for (int m = bw + 1; m <= size - cutoff; ++m) { // naive part
        // real part
        DLTNaive(rres + (m * size), bw, size - m, weights, fltres, seminaive_naive_table[size - m], scratchpad);
        // load real part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                rdataptr[i] = -fltres[i];
        } else {
            memcpy(rdataptr, fltres, sizeof(double) * (m - bw));
        }
        rdataptr += m - bw;

        // imaginary part
        DLTNaive(ires + (m * size), bw, size - m, weights, fltres, seminaive_naive_table[size - m], scratchpad);
        // load imaginary part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                idataptr[i] = -fltres[i];
        } else {
            memcpy(idataptr, fltres, sizeof(double) * (m - bw));
        }
        idataptr += m - bw;
    }

    for (int m = size - cutoff + 1; m < size; ++m) { // semi-naive part
        // real part
        DLTSemi(rres + (m * size), bw, size - m, fltres, scratchpad, seminaive_naive_table[size - m],
                weights, DCT_plan);
        // load real part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                rdataptr[i] = -fltres[i];
        } else {
            memcpy(rdataptr, fltres, sizeof(double) * (m - bw));
        }
        rdataptr += m - bw;

        // imaginary part
        DLTSemi(ires + (m * size), bw, size - m, fltres, scratchpad, seminaive_naive_table[size - m],
                weights, DCT_plan);
        // load imaginary part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                idataptr[i] = -fltres[i];
        } else {
            memcpy(idataptr, fltres, sizeof(double) * (m - bw));
        }
        idataptr += m - bw;
    }
}

/**
 * @brief Computes <b>inverse spherical harmonic transform</b>.
 *
 * @note See an example of use in test_s2_semi_memo.c and test_s2_semi_memo_inv.c
 *
 * @param rcoeffs array of length <tt>bw*bw</tt> of the real part of harmonic coefficients
 * @param icoeffs array of length <tt>bw*bw</tt> of the imaginary part of harmonic coefficients
 * @param rdata array of length <tt>4*bw*bw</tt> which will contain the real part of transformed
 * result
 * @param idata array of length <tt>4*bw*bw</tt> which will contain the imaginary part of
 * transformed result
 * @param bw bandwidth of problem
 * @param transpose_seminaive_naive_table pre-generated transposed spharmonic Pml table
 * @param workspace space for computations of size <tt>(8*bw^2)+(10*bw)</tt>
 * @param data_format determines the format of the computed data
 * @param cutoff determines order to switch from seminaive to naive algorithm
 * @param inv_DCT_plan plan for inverse DLT which is used as argument for call InvDLTSemi()
 * @param inv_FFT_plan plan for inverse FFT
 *
 * @note The real and imaginary part of harmonic coefficients should be stored in @p rcoeffs and
 * @p icoeffs before passing to this function in the same order as it was mentioned in FSTSemiMemo().
 * @note @p transpose_seminaive_naive_table should be generated by Transpose_Spharmonic_Pml_Table().
 * @note See more info about @p inv_DCT_plan in InvDLTSemi().
*/
void InvFSTSemiMemo(double* rcoeffs, double* icoeffs, double* rdata, double* idata, const int bw,
                    double** transpose_seminaive_naive_table, double* workspace, DataFormat data_format,
                    const int cutoff, fftw_plan* inv_DCT_plan, fftw_plan* inv_FFT_plan) {
    // TODO try to use double** -> double*
    int size = 2 * bw;

    // total workspace = (8 * bw^2) + (10 * bw)
    double* rfourdata = workspace;                  // needs (size * size)
    double* ifourdata = rfourdata + (size * size);  // needs (size * size)
    double* rinvfltres = ifourdata + (size * size); // needs (2 * bw)
    double* iminvfltres = rinvfltres + (size);      // needs (2 * bw)
    double* sin_values = iminvfltres + (size);      // needs (2 * bw)
    double* eval_pts = sin_values + (size);         // needs (2 * bw)
    double* scratchpad = eval_pts + (size);         // needs (2 * bw)

    // load up the sin_values array
    AcosOfChebyshevNodes(size, eval_pts);
    for (int i = 0; i < size; ++i)
        sin_values[i] = sin(eval_pts[i]);

    // do all of the inverse Legendre transforms
    double* rdataptr = rcoeffs;
    double* idataptr = icoeffs;

    for (int m = 0; m < cutoff; ++m) { // semi-naive part
        // real part
        InvDLTSemi(rdataptr, bw, m, rinvfltres, transpose_seminaive_naive_table[m], sin_values, scratchpad,
                   inv_DCT_plan);
        // imaginary part
        InvDLTSemi(idataptr, bw, m, iminvfltres, transpose_seminaive_naive_table[m], sin_values, scratchpad,
                   inv_DCT_plan);

        // store normal, then tranpose before doing inverse fft
        memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
        memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

        rdataptr += bw - m;
        idataptr += bw - m;
    }

    for (int m = cutoff; m < bw; ++m) { // naive part
        // real part
        InvDLTNaive(rdataptr, bw, m, rinvfltres, transpose_seminaive_naive_table[m]);
        // imaginary part
        InvDLTNaive(idataptr, bw, m, iminvfltres, transpose_seminaive_naive_table[m]);

        // store normal, then tranpose before doing inverse fft
        memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
        memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

        rdataptr += bw - m;
        idataptr += bw - m;
    }

    // fill in zero values where m = bw (from problem definition)
    memset(rfourdata + (bw * size), 0, sizeof(double) * size);
    memset(ifourdata + (bw * size), 0, sizeof(double) * size);

    /*
        If the data is real, we don't have to compute the coeffs whose order is less than 0,
        i.e. since the data is real, we know that

        invf-hat(l,-m) = conjugate(invf-hat(l,m)),

        so use that to get the rest of the coefficients
    */
    if (data_format == COMPLEX) {
        // do negative m values
        for (int m = bw + 1; m <= size - cutoff; ++m) { // naive part
            InvDLTNaive(rdataptr, bw, size - m, rinvfltres, transpose_seminaive_naive_table[size - m]);
            InvDLTNaive(idataptr, bw, size - m, iminvfltres, transpose_seminaive_naive_table[size - m]);

            // store normal, then tranpose before doing inverse fft
            if (m % 2)
                for (int i = 0; i < size; ++i) {
                    rinvfltres[i] *= -1.0;
                    iminvfltres[i] *= -1.0;
                }

            memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
            memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

            rdataptr += bw - (size - m);
            idataptr += bw - (size - m);
        }

        for (int m = size - cutoff + 1; m < size; ++m) { // semi-naive part
            InvDLTSemi(rdataptr, bw, size - m, rinvfltres, transpose_seminaive_naive_table[size - m],
                       sin_values, scratchpad, inv_DCT_plan);
            InvDLTSemi(idataptr, bw, size - m, iminvfltres, transpose_seminaive_naive_table[size - m],
                       sin_values, scratchpad, inv_DCT_plan);

            // store normal, then tranpose before doing inverse fft
            if (m % 2)
                for (int i = 0; i < size; ++i) {
                    rinvfltres[i] *= -1.0;
                    iminvfltres[i] *= -1.0;
                }

            memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
            memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

            rdataptr += bw - (size - m);
            idataptr += bw - (size - m);
        }
    } else { // real data
        for (int m = bw + 1; m < size; ++m) {
            memcpy(rfourdata + (m * size), rfourdata + ((size - m) * size), sizeof(double) * size);
            memcpy(ifourdata + (m * size), ifourdata + ((size - m) * size), sizeof(double) * size);

            for (int i = 0; i < size; ++i)
                ifourdata[(m * size) + i] *= -1.0;
        }
    }

    // normalize
    double normed_coeff = 1. / (sqrt(2. * M_PI));
    for (int i = 0; i < 4 * bw * bw; i++) {
        rfourdata[i] *= normed_coeff;
        ifourdata[i] *= normed_coeff;
    }

    fftw_execute_split_dft(*inv_FFT_plan, ifourdata, rfourdata, idata, rdata);
}

/**
 * @brief Computes <b>zonal harmonic transform</b> using seminaive algorithm.
 *
 * This transform is used in convolutions. Only computes spherical harmonics for <tt>m=0</tt>.
 * 
 * @param rdata array of length <tt>4*bw*bw</tt> of the real part of the data samples
 * @param idata array of length <tt>4*bw*bw</tt> of the imaginary part of the data samples
 * @param rres array of length <tt>bw</tt> which will contain the real part of harmonic
 * coefficients
 * @param ires array of length <tt>bw</tt> which will contain the imaginary part of harmonic
 * coefficients
 * @param bw bandwidth of problem
 * @param cos_pml_table pre-generated table of Legendre coefficients of P(0,l) functions
 * @param workspace space for computations of size <tt>12*bw</tt>
 * @param data_format determines the format of the computed data
 * @param DCT_plan plan for DLT which is used as argument for call DLTSemi()
 * @param weights array which is used as argument for call DLTSemi()
 *
 * @note @p cos_pml_table should be generated by GenerateCosPmlTable() with <tt>m = 0</tt>.
 * @note See more info about @p DCT_plan and @p weights in DLTSemi().
 */
void FZTSemiMemo(double* rdata, double* idata, double* rres, double* ires, const int bw, double* cos_pml_table,
                 double* workspace, const DataFormat data_format, fftw_plan* DCT_plan, double* weights) {
    int size = 2 * bw;

    double* r0 = workspace;         // needs (2 * bw)
    double* i0 = r0 + size;         // needs (2 * bw)
    double* scratchpad = i0 + size; // needs (4 * bw)

    // total workspace = 13*bw

    double dsize = sqrt(2. * M_PI) / size;
    // compute the m = 0 components
    for (int i = 0; i < size; ++i) {
        double tmpreal = 0.0;
        double tmpimag = 0.0;

        for (int j = 0; j < size; ++j) {
            tmpreal += rdata[(i * size) + j];
            tmpimag += idata[(i * size) + j];
        }

        // normalize
        r0[i] = tmpreal * dsize;
        i0[i] = tmpimag * dsize;
    }

    // real part
    DLTSemi(r0, bw, 0, rres, scratchpad, cos_pml_table, weights, DCT_plan);

    if (data_format == COMPLEX) // if data is complex, do imaginary part
        DLTSemi(i0, bw, 0, ires, scratchpad, cos_pml_table, weights, DCT_plan);
    else // otherwise set coefficients to zero
        memset(ires, 0, sizeof(double) * size);
}

/**
 * @brief <b>Convolves</b> two functions defined on the 2-sphere.
 *
 * Uses seminaive algorithms for spherical harmonic transforms.
 *
 * @note See an example of use in test_conv_semi_memo.c
 *
 * Conv2Sphere requires memory for spharmonic tables, local workspace and workspace for FSTSemiMemo() (reuses
 * in InvFSTSemiMemo()):\n
 * @code
 * legendre_size = Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw,cutoff)
 * 2*legendre_size + (8*(bw*bw)+10*bw) + (4*(bw*bw)+2*bw) = 2*legendre_size + (12*(bw*bw)+12*bw)
 * @endcode
 *
 * @param rdata array of length <tt>4*bw*bw</tt> of the real part of sampled function
 * @param idata array of length <tt>4*bw*bw</tt> of the imaginary part of sampled function
 * @param rfilter array of length <tt>4*bw*bw</tt> of the real part of sampled filter function
 * @param ifilter array of length <tt>4*bw*bw</tt> of the imaginary part of sampled filter function
 * @param rres array of length <tt>4*bw*bw</tt> which will contain the real part of result function
 * @param ires array of length <tt>4*bw*bw</tt> which will contain the imaginary part of result
 * function
 * @param bw bandwidth of problem
 * @param workspace space for computations of size <tt>(2*legendre_size)+(12*(bw*bw)+12*bw)</tt>
 *
 * @note We assume that data is real.
 * @note This function will do seminaive algorithm for all orders.
 * @note If you want to do multiple convolutions, you can freely reuse once allocated workspace and
 * generated spharmonic Pml tables. Same for FSTSemiMemo(), FZTSemiMemo() and InvFSTSemiMemo()
 * functions.
 */
void ConvOn2SphereSemiMemo(double* rdata, double* idata, double* rfilter, double* ifilter, double* rres, double* ires,
                           const int bw, double* workspace) {
    // TODO cutoff to args? (doc days it's legal but needs hardcode for changing param value now)
    int size = 2 * bw;
    int cutoff = bw;
    int legendre_size = Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff);

    double* spharmonic_result_space = workspace; // needs legendre_size
    double* transpose_spharmonic_result_space = spharmonic_result_space + legendre_size; // needs legendre_size

    double* frres = transpose_spharmonic_result_space + legendre_size; // needs (bw * bw)
    double* fires = frres + (bw * bw);                                 // needs (bw * bw)
    double* trres = fires + (bw * bw);                                 // needs (bw * bw)
    double* tires = trres + (bw * bw);                                 // needs (bw * bw)
    double* filtrres = tires + (bw * bw);                              // needs bw
    double* filtires = filtrres + bw;                                  // needs bw
    double* scratchpad = filtires + bw;                                // needs (8 * bw^2) + (10 * bw)

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    GenerateWeightsForDLT(bw, weights);

    // Make DCT plans. Note that I will be using the GURU interface to execute these plans within the routines
    fftw_plan DCT_plan = fftw_plan_r2r_1d(size, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_plan inv_DCT_plan = fftw_plan_r2r_1d(size, weights, rdata, FFTW_REDFT01, FFTW_ESTIMATE);

    // fftw "preamble"
    // Note that FFT plan places the output in a transposed array
    int rank = 1;
    fftw_iodim dims[rank];
    dims[0].n = size;
    dims[0].is = 1;
    dims[0].os = size;

    int howmany_rank = 1;
    fftw_iodim howmany_dims[howmany_rank];
    howmany_dims[0].n = size;
    howmany_dims[0].is = size;
    howmany_dims[0].os = 1;

    fftw_plan FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                  workspace + (4 * bw * bw), FFTW_ESTIMATE);

    // Note that FFT plans assumes that I'm working with a transposed array, e.g. the inputs for a length 2*bw transform
    // are placed every 2*bw apart, the output will be consecutive entries in the array

    rank = 1;
    dims[0].n = size;
    dims[0].is = size;
    dims[0].os = 1;

    howmany_rank = 1;
    howmany_dims[0].n = size;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = size;

    fftw_plan inv_FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                      workspace + (4 * bw * bw), FFTW_ESTIMATE);

    // precompute the associated Legendre fcts
    double** spharmonic_pml_table = Spharmonic_Pml_Table(bw, spharmonic_result_space, scratchpad);
    double** transpose_spharmonic_pml_table =
        Transpose_Spharmonic_Pml_Table(spharmonic_pml_table, bw, transpose_spharmonic_result_space);

    FSTSemiMemo(rdata, idata, frres, fires, bw, spharmonic_pml_table, scratchpad, 1, bw, &DCT_plan, &FFT_plan, weights);
    FZTSemiMemo(rfilter, ifilter, filtrres, filtires, bw, spharmonic_pml_table[0], scratchpad, 1, &DCT_plan, weights);

    TransMult(frres, fires, filtrres, filtires, trres, tires, bw);

    InvFSTSemiMemo(trres, tires, rres, ires, bw, transpose_spharmonic_pml_table, scratchpad, 1, bw, &inv_DCT_plan,
                   &inv_FFT_plan);

    free(weights);
    free(spharmonic_pml_table);
    free(transpose_spharmonic_pml_table);

    fftw_destroy_plan(inv_FFT_plan);
    fftw_destroy_plan(FFT_plan);
    fftw_destroy_plan(inv_DCT_plan);
    fftw_destroy_plan(DCT_plan);
}
