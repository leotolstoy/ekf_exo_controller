#!/usr/bin/env python
# import CBoostedEKF as cpt
from .CBoostedEKF import *
import numpy as np

def loadBezierCurves(filename):
    # doesn't need to be a member function!
    data = np.loadtxt(filename,delimiter=',')

    best_fit_params_footAngle = data[0,:]
    best_fit_params_shankAngle = data[1,:]

    return (best_fit_params_footAngle,best_fit_params_shankAngle)


class TorqueProfile():
    """ A Boosted Torque Profile Class """
    def __init__(self, model_filepath='TorqueProfile/torqueProfileCoeffs_dataport3P.csv'):
        self.phaseDelins = np.array([0.1,0.5,0.65,1])
        data = np.loadtxt(model_filepath,delimiter=',')
        self.best_fit_params_torque = data[:]
        self.bgmw = BoostedGaitModelWrapper(self.phaseDelins.reshape((4,1)), self.best_fit_params_torque.reshape((64,1)), self.best_fit_params_torque.reshape((64,1)))

    def evalTorqueProfile(self,phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE):
        # return self.profileMesh(phase_estimate_PE, incline_estimate_PE)[0]

        return self.evalBiologicalProfile(phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE)

    def evalBiologicalProfile(self,phase_estimate_PE, stepLength_estimate_PE, incline_estimate_PE):
        torque = 0
        if phase_estimate_PE <= 0.65:
            torque = 1/5.* self.bgmw.getFootAngle(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)

        if torque < 0:
            torque = 0

        return torque
        

class GaitModel():
    """ A Boosted Gait-Model Class """

    def __init__(self, model_filepath='BestFitParams/regressionMatrices_dataport3P.csv'):
        self.phaseDelins = np.array([[0.1,0.5,0.65,1]]).T # Vector4d
        (self.best_fit_params_footAngle, self.best_fit_params_shankAngle) = loadBezierCurves(model_filepath)
        self.bgmw = BoostedGaitModelWrapper(self.phaseDelins.reshape((4,1)), self.best_fit_params_footAngle.reshape((64,1)), self.best_fit_params_shankAngle.reshape((64,1)))

    ## Accounts for my renaming of all these functions
    def returnFootAngle(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getFootAngle(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE, self.best_fit_params_footAngle, self.phaseDelins)

    def returnShankAngle(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getShankAngle(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE, self.best_fit_params_shankAngle, self.phaseDelins);

    def returnFootAngleDeriv_dphase(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdFootAngle_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngleDeriv_dphase(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdShankAngle_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);

    def returnFootAngleDeriv_dsL(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdFootAngle_dsl(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dsL(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngleDeriv_dsL(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdShankAngle_dsl(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dsL(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);

    def returnFootAngleDeriv_dincline(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdFootAngle_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngleDeriv_dincline(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getdShankAngle_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3PDeriv_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);
        

    def returnFootAngle2ndDeriv_dphase2(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddFootAngle_dphase_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphase2(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngle2ndDeriv_dphase2(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddShankAngle_dphase_dphase(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphase2(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);

    def returnFootAngle2ndDeriv_dphasedsL(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddFootAngle_dphase_dsl(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphasedsL(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngle2ndDeriv_dphasedsL(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddShankAngle_dphase_dsl(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphasedsL(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);

    def returnFootAngle2ndDeriv_dphasedincline(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddFootAngle_dphase_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphasedincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_footAngle, self.phaseDelins);

    def returnShankAngle2ndDeriv_dphasedincline(self, phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE):
        return self.bgmw.getddShankAngle_dphase_dincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE)
        # PiecewiseBezier3P_2ndDeriv_dphasedincline(phase_estimate_PE,stepLength_estimate_PE,incline_estimate_PE,self.best_fit_params_shankAngle, self.phaseDelins);


def InstantiateBGMW(model_filepath='BestFitParams/regressionMatrices_dataport3P.csv'):
    # If you want to use my interface, we just need a special constructor for the BGMW
    phaseDelins = [0.1,0.5,0.65,1]
    (best_fit_params_footAngle, best_fit_params_shankAngle) = loadBezierCurves(model_filepath)
    bgmw = BoostedGaitModelWrapper(phaseDelins, best_fit_params_footAngle, best_fit_params_shankAngle)
    return bgmw





def hello():
    print("Hello world 42!\n\t-pyCppWrapper")




