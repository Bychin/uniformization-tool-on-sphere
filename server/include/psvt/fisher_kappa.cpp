#include <cmath>

const int kThetaIndex = 0;
const int kProbIndex = 1;
const int kKappaIndex = 2;

const int kMaxKappa = 40;
const int kMaxEvalKappa = kMaxKappa;

const double kAccuracy = 1.E-13;
const double kActualZero = 1.E-18;
const double kEvalAccuracy = 1.E-7;

int Fisher_Params(double* Par, int field) {
    double v_tmp[6] = {0., 0., 0., 0., 0., 0.};

    switch (field) {
    case kThetaIndex:
        if (Par[kProbIndex] < 1. && Par[kProbIndex] >= 0) {
            if (Par[kKappaIndex] > kActualZero && Par[kKappaIndex] < kMaxEvalKappa) {
                v_tmp[0] = std::exp(Par[kKappaIndex]);
                Par[kThetaIndex] =
                    std::acos(std::log(v_tmp[0] - Par[kProbIndex] * (v_tmp[0] - 1. / v_tmp[0])) / Par[kKappaIndex]);
            } else if (Par[kKappaIndex] < kActualZero)
                Par[kThetaIndex] = std::acos(1. - 2. * Par[kProbIndex]);
            else
                Par[kThetaIndex] = std::acos(1. + std::log(1 - Par[kProbIndex]) / Par[kKappaIndex]);

            return 0;
        };
        break;
    case kProbIndex:
        if (Par[kKappaIndex] >= 0. && Par[kThetaIndex] > kActualZero) {
            if (Par[kKappaIndex] < kMaxKappa) {
                if (Par[kKappaIndex] > kActualZero) {
                    v_tmp[0]        = std::exp(Par[kKappaIndex]);
                    Par[kProbIndex] = 1. / (v_tmp[0] - 1. / v_tmp[0]) *
                                      (v_tmp[0] - std::pow(v_tmp[0], std::cos(Par[kThetaIndex])));
                } else
                    Par[kProbIndex] = (1. - std::cos(Par[kThetaIndex])) / 2.;
            } else {
                v_tmp[0] = Par[kKappaIndex] * (std::cos(Par[kThetaIndex]) - 1.);
                if (std::exp(v_tmp[0]) < kActualZero)
                    Par[kProbIndex] = 1.;
                else
                    Par[kProbIndex] = 1. - std::exp(Par[kKappaIndex] * (std::cos(Par[kThetaIndex]) - 1.));
            }
            return 0;
        }
        break;
    case kKappaIndex:
        if (Par[kThetaIndex] > 0 && Par[kProbIndex] >= (1. - std::cos(Par[kThetaIndex])) / 2.) {
            if (Par[kProbIndex] - (1. - std::cos(Par[kThetaIndex])) / 2. < kAccuracy) {
                Par[kKappaIndex] = 0.;
                return 0;
            };

            v_tmp[2] = std::cos(Par[kThetaIndex]);

            do {
                v_tmp[1] += 1.;
                if (v_tmp[1] > kMaxKappa)
                    v_tmp[5] = (1. - std::exp(v_tmp[1] * (v_tmp[2] - 1)));
                else {
                    v_tmp[3] = std::exp(v_tmp[1]);
                    v_tmp[5] = 1. / (1. - (1. / v_tmp[3]) * (1. / v_tmp[3])) *
                               (1. - std::exp(v_tmp[1] * (v_tmp[2] - 1.)));
                }
            } while (v_tmp[5] < Par[kProbIndex]);

            do {
                v_tmp[4] = (v_tmp[1] + v_tmp[0]) / 2.;
                if (v_tmp[4] > kMaxKappa)
                    v_tmp[5] = (1. - std::exp(v_tmp[4] * (v_tmp[2] - 1)));
                else {
                    v_tmp[3] = std::exp(v_tmp[4]);
                    v_tmp[5] = 1. / (1. - (1. / v_tmp[3]) * (1. / v_tmp[3])) *
                               (1. - std::exp(v_tmp[4] * (v_tmp[2] - 1)));
                }
                if (v_tmp[5] > Par[kProbIndex])
                    v_tmp[1] = v_tmp[4];
                else
                    v_tmp[0] = v_tmp[4];

            } while (v_tmp[1] - v_tmp[0] > kEvalAccuracy * v_tmp[1]);
            Par[kKappaIndex] = 0.5 * (v_tmp[1] + v_tmp[0]);
            return 0;
        }
        break;
    }

    return 1;
}

double P95FisherKappa(double alpha95Rad) {
    double FParams[3];
    FParams[kProbIndex]  = 0.95;
    FParams[kThetaIndex] = std::fabs(alpha95Rad);
    Fisher_Params(FParams, kKappaIndex);
    return FParams[kKappaIndex];
}

double Alpha95InDegreesToFisherKappa(double alpha95) {
    double alpha95InRadians = alpha95 * M_PI / 180.;
    return P95FisherKappa(alpha95InRadians);
}

double P95FisherAlpha95Rad(double FisherKappa) {
    double FParams[3];
    FParams[kProbIndex]  = 0.95;
    FParams[kKappaIndex] = std::fabs(FisherKappa);
    Fisher_Params(FParams, kThetaIndex);
    return FParams[kThetaIndex];
}

void SetP95FisherDefaults(double* FParams) {
    FParams[kKappaIndex] = kMaxEvalKappa;
    FParams[kThetaIndex] = 0.0;
    FParams[kProbIndex]  = 0.95;
}
