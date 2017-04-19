function [kalman] = Create_New_Kalman(init_state_vec,Q,R)
% This function creats a new predictor for an object.

StateTransitionModel = [
                        1 0 0 0 0 1 0 0 0 0
                        0 1 0 0 0 0 1 0 0 0
                        0 0 1 0 0 0 0 0 0 0 
                        0 0 0 1 0 0 0 0 0 0
                        0 0 0 0 1 0 0 0 0 0
                        0 0 0 0 0 1 0 0 0 0
                        0 0 0 0 0 0 1 0 0 0
                        0 0 0 0 0 0 0 1 0 0
                        0 0 0 0 0 0 0 0 1 0
                        0 0 0 0 0 0 0 0 0 1
                       ];
MeasurementModel = [
                        
                        1 0 0 0 0 0 0 0 0 0
                        0 1 0 0 0 0 0 0 0 0
                        0 0 1 0 0 0 0 0 0 0
                        0 0 0 1 0 0 0 0 0 0
                        0 0 0 0 1 0 0 0 0 0
                        0 0 0 0 0 1 0 0 0 0
                        0 0 0 0 0 0 1 0 0 0
                        0 0 0 0 0 0 0 1 0 0
                        0 0 0 0 0 0 0 0 1 0
                        0 0 0 0 0 0 0 0 0 1
                        ];
               
kalman = vision.KalmanFilter(StateTransitionModel,MeasurementModel,'State',init_state_vec,'ProcessNoise',Q,'MeasurementNoise',R);
