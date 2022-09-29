#ifndef BEZIER_CLASS
#define BEZIER_CLASS

#include <Eigen/Dense>

namespace gm{ // gait filter, not to be confused with af, for attitude filter


typedef Eigen::Matrix<float, 64, 1> Vector64f;
typedef Eigen::Matrix<float, 16, 1> Vector16f;
typedef Eigen::Matrix<float, 4, 1> Vector4f;

class BezierClass
{
protected: // member variables, only availiable to instances
	Vector64f bezierCoeffs64;
    Vector4f phaseDelins;

public:
    BezierClass();
    BezierClass(const Vector64f&, const Vector4f&); // constructor
    Vector16f getBezierCoeffs(float);
    // evaluation functions

    float returnBezierEval3P(float, float, float);
    float returnBezier3PDerivEval_dphase(float, float, float);
    float returnBezier3PDerivEval_dsL(float, float, float);
    float returnBezier3PDerivEval_dincline(float, float, float);
    float returnBezier3P_2ndDerivEval_dphase2(float, float, float);
    float returnBezier3P_2ndDerivEval_dphasedsL(float, float, float);
    float returnBezier3P_2ndDerivEval_dphasedincline(float, float, float);
};

} // gf
#endif