#include "Bezier.hpp"
namespace gm{



BezierClass::BezierClass(){
}

BezierClass::BezierClass(const Vector64f& bezierCoeffs_arg, const Vector4f& phaseDelins_arg){
    bezierCoeffs64 = bezierCoeffs_arg;
    phaseDelins = phaseDelins_arg;
}

Vector16f BezierClass::getBezierCoeffs(float phase){
	Vector16f h_selected;
    if (phase < 0){
        phase = 0;
    }
        
    else if (phase > 1){
        phase = 1;
    }
    
    if (phase <= phaseDelins[0]){
        // phaseIdxs = range(0 + 0*16,16 + 0*16)
        h_selected = bezierCoeffs64.segment<16>(0*16);
    }
    else if (phase <= phaseDelins[1]){
        h_selected = bezierCoeffs64.segment<16>(1*16);
    }
    else if (phase <= phaseDelins[2]){
        h_selected = bezierCoeffs64.segment<16>(2*16);
    }
    else{
        h_selected = bezierCoeffs64.segment<16>(3*16);
    }

    

    return h_selected;

}

// evaluation functions

float BezierClass::returnBezierEval3P(float phase, float stepLength, float incline){

	Vector16f bezierEval3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);



    bezierEval3P << \
        (incline)*stepLength* pow((1-phase),3),\
        (incline)*stepLength*3* pow((1-phase),2) *phase,\
        (incline)*stepLength*3*(1-phase)* pow(phase,2),\
        (incline)*stepLength* pow(phase,3),\
        (incline)*(1 - stepLength)*pow((1-phase),3),\
        (incline)*(1 - stepLength)*3*pow((1-phase),2)*phase,\
        (incline)*(1 - stepLength)*3*(1-phase)*pow(phase,2),\
        (incline)*(1 - stepLength)*pow(phase,3),\
        (1 - incline)*stepLength*pow((1-phase),3),\
        (1 - incline)*stepLength*3*pow((1-phase),2)*phase,\
        (1 - incline)*stepLength*3*(1-phase)*pow(phase,2),\
        (1 - incline)*stepLength*pow(phase,3),\
        (1 - incline)*(1 - stepLength)*pow((1-phase),3),\
        (1 - incline)*(1 - stepLength)*3*pow((1-phase),2)*phase,\
        (1 - incline)*(1 - stepLength)*3*(1-phase)*pow(phase,2),\
        (1 - incline)*(1 - stepLength)*pow(phase,3)\
        ;



    return bezierCoeffs.dot(bezierEval3P);
}

float BezierClass::returnBezier3PDerivEval_dphase(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);

    bezierDerivCoeffs3P << \
        (incline)*stepLength*-3* pow(1-phase,2),\
        (incline)*stepLength* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (incline)*stepLength*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (incline)*stepLength*3* pow(phase,2),\
        (incline)*(1 - stepLength)*-3* pow(1-phase,2), \
        (incline)*(1 - stepLength)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (incline)*(1 - stepLength)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (incline)*(1 - stepLength)*3* pow(phase,2),\
        (1 - incline)*stepLength*-3* pow(1-phase,2),\
        (1 - incline)*stepLength* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (1 - incline)*stepLength*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1 - incline)*stepLength*3* pow(phase,2),\
        (1 - incline)*(1 - stepLength)*-3* pow(1-phase,2), \
        (1 - incline)*(1 - stepLength)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (1 - incline)*(1 - stepLength)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1 - incline)*(1 - stepLength)*3* pow(phase,2)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}

float BezierClass::returnBezier3PDerivEval_dsL(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);

    bezierDerivCoeffs3P << \
        (incline)*1* pow((1-phase),3),\
        (incline)*1*3* pow((1-phase),2) *phase,\
        (incline)*1*3*(1-phase)* pow(phase,2),\
        (incline)*1* pow(phase,3),\
        (incline)*(-1)*pow((1-phase),3),\
        (incline)*(-1)*3*pow((1-phase),2)*phase,\
        (incline)*(-1)*3*(1-phase)*pow(phase,2),\
        (incline)*(-1)*pow(phase,3),\
        (1 - incline)*1*pow((1-phase),3),\
        (1 - incline)*1*3*pow((1-phase),2)*phase,\
        (1 - incline)*1*3*(1-phase)*pow(phase,2),\
        (1 - incline)*1*pow(phase,3),\
        (1 - incline)*(-1)*pow((1-phase),3),\
        (1 - incline)*(-1)*3*pow((1-phase),2)*phase,\
        (1 - incline)*(-1)*3*(1-phase)*pow(phase,2),\
        (1 - incline)*(-1)*pow(phase,3)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}


float BezierClass::returnBezier3PDerivEval_dincline(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);


    bezierDerivCoeffs3P << \
        (1)*stepLength* pow((1-phase),3),\
        (1)*stepLength*3* pow((1-phase),2) *phase,\
        (1)*stepLength*3*(1-phase)* pow(phase,2),\
        (1)*stepLength* pow(phase,3),\
        (1)*(1 - stepLength)*pow((1-phase),3),\
        (1)*(1 - stepLength)*3*pow((1-phase),2)*phase,\
        (1)*(1 - stepLength)*3*(1-phase)*pow(phase,2),\
        (1)*(1 - stepLength)*pow(phase,3),\
        (-1)*stepLength*pow((1-phase),3),\
        (-1)*stepLength*3*pow((1-phase),2)*phase,\
        (-1)*stepLength*3*(1-phase)*pow(phase,2),\
        (-1)*stepLength*pow(phase,3),\
        (-1)*(1 - stepLength)*pow((1-phase),3),\
        (-1)*(1 - stepLength)*3*pow((1-phase),2)*phase,\
        (-1)*(1 - stepLength)*3*(1-phase)*pow(phase,2),\
        (-1)*(1 - stepLength)*pow(phase,3)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}


float BezierClass::returnBezier3P_2ndDerivEval_dphase2(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);

    bezierDerivCoeffs3P << \
        (incline)*stepLength*(6 * (1 - phase)), \
        (incline)*stepLength* ( -2*(6*(1 - phase)) +  (6*phase) ), \
        (incline)*stepLength*( 6*(1 - phase) + (-2 * (6 * phase)) ), \
        (incline)*stepLength*(6 * phase),\
        (incline)*(1 - stepLength)*(6 * (1 - phase)), \
        (incline)*(1 - stepLength)*( -2*(6*(1 - phase)) +  (6*phase) ), \
        (incline)*(1 - stepLength)* ( 6*(1 - phase) + (-2 * (6 * phase)) ),\
        (incline)*(1 - stepLength)*(6 * phase),\
        (1 - incline)*stepLength*(6 * (1 - phase)), \
        (1 - incline)*stepLength* ( -2*(6*(1 - phase)) +  (6*phase) ), \
        (1 - incline)*stepLength*( 6*(1 - phase) + (-2 * (6 * phase)) ), \
        (1 - incline)*stepLength*(6 * phase),\
        (1 - incline)*(1 - stepLength)*(6 * (1 - phase)), \
        (1 - incline)*(1 - stepLength)*( -2*(6*(1 - phase)) +  (6*phase) ), \
        (1 - incline)*(1 - stepLength)* ( 6*(1 - phase) + (-2 * (6 * phase)) ),\
        (1 - incline)*(1 - stepLength)*(6 * phase)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}

float BezierClass::returnBezier3P_2ndDerivEval_dphasedsL(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);

    bezierDerivCoeffs3P << \
        (incline)*1*-3* pow(1-phase,2),\
        (incline)*1* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (incline)*1*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (incline)*1*3* pow(phase,2),\
        (incline)*(-1)*-3* pow(1-phase,2), \
        (incline)*(-1)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (incline)*(-1)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (incline)*(-1)*3* pow(phase,2),\
        (1 - incline)*1*-3* pow(1-phase,2),\
        (1 - incline)*1* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (1 - incline)*1*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1 - incline)*1*3* pow(phase,2),\
        (1 - incline)*(-1)*-3* pow(1-phase,2), \
        (1 - incline)*(-1)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (1 - incline)*(-1)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1 - incline)*(-1)*3* pow(phase,2)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}


float BezierClass::returnBezier3P_2ndDerivEval_dphasedincline(float phase, float stepLength, float incline){

    Vector16f bezierDerivCoeffs3P;
    Vector16f bezierCoeffs = getBezierCoeffs(phase);

    bezierDerivCoeffs3P << \
        (1)*stepLength*-3* pow(1-phase,2),\
        (1)*stepLength* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (1)*stepLength*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1)*stepLength*3* pow(phase,2),\
        (1)*(1 - stepLength)*-3* pow(1-phase,2), \
        (1)*(1 - stepLength)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (1)*(1 - stepLength)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (1)*(1 - stepLength)*3* pow(phase,2),\
        (-1)*stepLength*-3* pow(1-phase,2),\
        (-1)*stepLength* (3* pow(1-phase,2) - 6*(1 - phase)*phase),\
        (-1)*stepLength*(6*(1 - phase)*phase - 3* pow(phase,2)),\
        (-1)*stepLength*3* pow(phase,2),\
        (-1)*(1 - stepLength)*-3* pow(1-phase,2), \
        (-1)*(1 - stepLength)*(3* pow(1-phase,2) - 6*(1 - phase)*phase), \
        (-1)*(1 - stepLength)* (6*(1 - phase)*phase - 3* pow(phase,2)),\
        (-1)*(1 - stepLength)*3* pow(phase,2)\
        ;

    return bezierCoeffs.dot(bezierDerivCoeffs3P);
}



}//gf