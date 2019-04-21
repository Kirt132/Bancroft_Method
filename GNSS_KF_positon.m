function [est_r_ea_e, est_clock, x_old_old, P_old_old] = GNSS_KF_positon(...
    GNSS_measurements,no_GNSS_meas,x_old_old,P_old_old,F_matrix,ProcessCOV_matrix, R)
% GNSS_Kalman_position - Calculates position, clock offset, 
% using weighted iterated least squares with priori.
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 11/4/2012 by Paul Groves
%
% Inputs:
%   GNSS_measurements     GNSS measurement data:
%     Column 1              Pseudo-range measurements (m)
%     Column 2              Pseudo-range rate measurements (m/s)
%     Columns 3-5           Satellite ECEF position (m)
%     Columns 6-8           Satellite ECEF velocity (m/s)
%     no_GNSS_meas          Number of satellites for which measurements are
%                         supplied
%   predicted_r_ea_e      prior predicted ECEF user position (m)
%
% Outputs:
%   est_r_ea_e            estimated ECEF user position (m)
%   est_clock             estimated receiver clock offset (m) and drift (m/s)
%   GDOP                  error measurement
 
% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details
%
%   P. Closas (2018): version only for position and clock offset LS estimation
% 

% Constants (sone of these could be changed to inputs at a later date)
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

% Begins

% POSITION AND CLOCK OFFSET

% Setup predicted state

test_convergence = 1;

% initialize
pred_meas = zeros(no_GNSS_meas,1);
H_matrix = zeros(no_GNSS_meas,4);


% Repeat until convergence (i.e. two consecutive iterations provide almost identical results)
%while test_convergence>0.0001
    
%   num_iter = num_iter + 1;
    
    % Loop measurements
    for j = 1:no_GNSS_meas

        % Predict approx range 
        delta_r = GNSS_measurements(j,3:5)' - x_old_old(1:3);
        approx_range = sqrt(delta_r' * delta_r);

        % Calculate frame rotation during signal transit time using (8.36)
        C_e_I = [1, omega_ie * approx_range / c, 0;...
                 -omega_ie * approx_range / c, 1, 0;...
                 0, 0, 1];

        % Predict pseudo-range using (9.143)
        delta_r = C_e_I *  GNSS_measurements(j,3:5)' - x_old_old(1:3);
        range = sqrt(delta_r' * delta_r);
        pred_meas(j,1) = range + x_old_old(4);
        
        % Predict line of sight and deploy in measurement matrix, (9.144)
        H_matrix (j,1:3) = - delta_r' / range;
        H_matrix (j,4) = 1;
        
    end % for j
   
    %calculate the new Kalman Values
    x_new_old = F_matrix*x_old_old;
    P_new_old = F_matrix*P_old_old*F_matrix'+ProcessCOV_matrix;
    K = P_new_old*H_matrix'*(H_matrix*P_new_old*H_matrix'+R)^(-1);
    x_est = x_new_old + K*(GNSS_measurements(1:no_GNSS_meas,1) - pred_meas(1:no_GNSS_meas));
    P_old_old = (eye(4)- K*H_matrix)*P_new_old;
    % Set predictions to estimates for next iteration

    
%end % while

% Set outputs to estimates
x_old_old = x_est;
est_r_ea_e(1:3,1) = x_est(1:3);
est_clock(1) = x_est(4);
covariancematrix = 4*inv(H_matrix(1:no_GNSS_meas,:)'*H_matrix(1:no_GNSS_meas,:));%get the values of sigma from this matrix
GDOP = sqrt(covariancematrix(1,1)+covariancematrix(2,2)+covariancematrix(3,3)+covariancematrix(4,4))/2;

