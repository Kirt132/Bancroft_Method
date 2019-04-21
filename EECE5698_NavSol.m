%
%   EECE5698-ST: GNSS signal processing
%       Simple pseudorange simulator to test navigation solvers
%
%   Pau Closas, Spring 2018 (using matlab functions from SW accompanying P. Groves book)

clearvars
close all
clc

rad2deg = 180/pi;
deg2rad = pi/180;


%% Simulation Parameters 

% number of epochs or pseudoranges
GNSS_config.no_epochs = 200;          
% time between consecutive pseudorange measurements (s)
GNSS_config.sampling = 1;          
% number of satellites
GNSS_config.no_sat = 30;              
% Mask angle (deg)
GNSS_config.mask_angle = 10;
% Code tracking error SD at zenith (m)
GNSS_config.code_track_err_SD = 2;
% Range rate tracking error SD at zenith (m/s)
GNSS_config.rate_track_err_SD = 0.02;
% Receiver clock offset at time=0 (m);
GNSS_config.rx_clock_offset = 10000;
% Receiver clock drift at time=0 (m/s);
GNSS_config.rx_clock_drift = 100;
% Initial estimated position (meters, ECEF)
GNSS_config.init_est_p_eb_ecef = [0;0;0];
% GNSS_config.init_est_p_eb_ecef = [1.5302;-4.467;4.2734]*1e6;  
% GNSS_config.init_est_p_eb_ecef = [1;1;1]*1e7;  

% specify location (static/ intitial) - for instance: 42°20'13.4"N 71°05'25.9"W
true_phi_b_dms = [42 20 13.4];        % latitude (degree, minutes, seconds)
true_lambda_b_dms = -[71 05 25.9];    % longitude (degree, minutes, seconds)
true_h_b = 10;                        % height (meters)
% some transformations
true_phi_b_deg = true_phi_b_dms(1) + true_phi_b_dms(2)/60 + true_phi_b_dms(3)/3600;                % latitude (decimal degrees)
true_phi_b_rad = deg2rad*true_phi_b_deg;                                                                   % latitude (radians)

true_lambda_b_deg = true_lambda_b_dms(1) + true_lambda_b_dms(2)/60 + true_lambda_b_dms(3)/3600;    % longitude (decimal degrees)
true_lambda_b_rad = deg2rad*true_lambda_b_deg;                                                             % longitude (radians)

true_p_eb_ecef = pv_NED_to_ECEF(true_phi_b_rad,true_lambda_b_rad,true_h_b,[0;0;0]);                     % Ellipsoidal-to-Cartesian
%set intital position to anitpodes
%GNSS_config.init_est_p_eb_ecef = -true_p_eb_ecef;
%set initial position to a few hundred meters away (200x 100y away)
expected_pos = true_p_eb_ecef + [200;100;0];
%GNSS_config.init_est_p_eb_ecef = expected_pos;
% norm(true_p_eb_ecef)        % should be comparable to the Earh radius (R_0=6378 km < norm(position_ecef) < R_P=6356 km)

% % % nice google map plot
% figure, plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20)
%     xlabel('longitude [^o]'), ylabel('latitude [^o]')
%     plot_google_map

% user's velocity (Cartesian, ECEF)
true_v_eb_ecef = [0 0 0]';          % static receiver (m/s)

% initialize variables
time = 0;
dummy = 0;
pseudorange = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
pseudorangerate = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
elevations = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
est_p_eb_ecef_LS = zeros(3,GNSS_config.no_epochs+1);
est_clock_LS = zeros(1,GNSS_config.no_epochs);
est_phi_b_LS = zeros(1,GNSS_config.no_epochs);
est_lambda_b_LS = zeros(1,GNSS_config.no_epochs);
est_h_b_LS = zeros(1,GNSS_config.no_epochs);
est_p_eb_ecef_WLS = zeros(3,GNSS_config.no_epochs+1);
est_clock_WLS = zeros(1,GNSS_config.no_epochs);
est_phi_b_WLS = zeros(1,GNSS_config.no_epochs);
est_lambda_b_WLS = zeros(1,GNSS_config.no_epochs);
est_h_b_WLS = zeros(1,GNSS_config.no_epochs);
est_p_eb_ecef_WLSA = zeros(3,GNSS_config.no_epochs+1);
est_clock_WLSA = zeros(1,GNSS_config.no_epochs);
est_phi_b_WLSA = zeros(1,GNSS_config.no_epochs);
est_lambda_b_WLSA = zeros(1,GNSS_config.no_epochs);
est_h_b_WLSA = zeros(1,GNSS_config.no_epochs);
est_p_eb_ecef_KF = zeros(3,GNSS_config.no_epochs+1);
est_clock_KF = zeros(1,GNSS_config.no_epochs);
est_phi_b_KF = zeros(1,GNSS_config.no_epochs);
est_lambda_b_KF = zeros(1,GNSS_config.no_epochs);
est_h_b_KF = zeros(1,GNSS_config.no_epochs);
real_phi_b = zeros(1,GNSS_config.no_epochs);
real_lambda_b = zeros(1,GNSS_config.no_epochs);
real_h_b = zeros(1,GNSS_config.no_epochs);
true_p_eb_ecef_vec = true_p_eb_ecef(:,1)*ones(1,GNSS_config.no_epochs);
predicted_X = [1530400, -4466900, 4273954,10000]'; % predicted value for Boston
Q_matrix = diag([4000,4000,120,1000]); %variance expected in Boston
x_old_old = predicted_X;
P_old_old = Q_matrix;
KF_est = zeros(4, GNSS_config.no_epochs);
KF_cov = zeros(GNSS_config.no_epochs,4);
F_Dynamics = eye(8); %define dynamics matrix for the reciever
F_Dynamics(1:3,4:6) = GNSS_config.sampling * eye(3);
velocity = [0,0,0]; %initial Reciever Velocity
cdt = 10000; %intial time offset
cdotdt = 0; % initial time offset velocity
%% main loop
for epoch = 1:GNSS_config.no_epochs    
    
    %% Generate pseudoranges
    Q = diag(diag(randn(8)))*[20,20,3,1,1,.5,150,50]'; %define covariance and multiply by noise
    x_old = [true_p_eb_ecef_vec(:,epoch)',velocity,cdt,cdotdt];
    x_new = F_Dynamics*x_old' + Q;
    true_p_eb_ecef_vec(:,epoch+1) = x_new(1:3);
    velocity = x_new(4:6)';
    cdt = x_new(7);
    cdotdt = x_new(8);
    
    % Determine satellite positions and velocities
    [sat_pos_es_e,sat_vel_es_e] = Satellite_positions_and_velocities(time,GNSS_config);
    % norm(sat_vel_es_e(1,:))        % magnitude of velocity for a satellite (m/s)

    % Generate GNSS measurements
    [GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(...
        time,sat_pos_es_e,sat_vel_es_e,true_p_eb_ecef_vec(:,epoch),true_phi_b_rad,true_lambda_b_rad,true_v_eb_ecef,GNSS_config);

    pseudorange(:,epoch) = GNSS_measurements(:,1);          % for plotting
    pseudorangerate(:,epoch) = GNSS_measurements(:,2);
    elevations(:,epoch) = GNSS_measurements(:,9);
    GNSS_config.rx_clock_offset = cdt;
    %% Navigation solution
    
    %% LS
    [est_p_eb_ecef_LS(:,epoch),est_clock_LS(epoch),num_iter_LS, GDOP_LS(epoch)] = GNSS_LS_position(GNSS_measurements,no_GNSS_meas,GNSS_config.init_est_p_eb_ecef);
    
    %% WLS
    %Generate R Matrix
    %First find the standard error divided by the sin(elevation) which was used to
    %generate noise, then square the elements this matrix to get sigma
    %squared for each satellite.
    R = diag(GNSS_config.code_track_err_SD./sin(GNSS_measurements(:,9))); 
    R= R.^2;
    %Get rid of the rows and columns of zeros
    R([(no_GNSS_meas+1):GNSS_config.no_sat],:) = [];
    R(:,[(no_GNSS_meas+1):GNSS_config.no_sat]) = [];
    %invert R Matrix to find W Matrix
    W=inv(R);
    %calculate WLS solution
    [est_p_eb_ecef_WLS(:,epoch),est_clock_WLS(epoch),num_iter_WLS, GDOP_WLS(epoch)] = GNSS_WLS_position(GNSS_measurements,no_GNSS_meas,GNSS_config.init_est_p_eb_ecef, W);

    %% WLSA
    %this uses the same R and W matrixes as above but uses a new covariance
    %and estimate defined here.
    
    [est_p_eb_ecef_WLSA(:,epoch),est_clock_WLSA(epoch),num_iter_WLSA, GDOP_WLSA(epoch)] = GNSS_WLSA_position(GNSS_measurements,no_GNSS_meas,predicted_X,W,Q_matrix);
    
    
    %% KF
    %This implements the karman filter which doesn't iterate allowing
    %multiple iterations to occur.
    F_matrix = eye(4);
    ProcessCOV_matrix = diag([20,20,3,150]);
    %R = 10*eye(no_GNSS_meas);
    R = R.^.5;
    KF_est(:, epoch) = x_old_old;
    KF_cov(epoch, :) = diag(P_old_old);
    [est_p_eb_ecef_KF(:,epoch), est_clock_KF(epoch), x_old_old, P_old_old] = GNSS_KF_positon(...
    GNSS_measurements,no_GNSS_meas,x_old_old,P_old_old,F_matrix,ProcessCOV_matrix, R);
    

%% transform estimates to latitude, longitude, and height
    [est_phi_b_LS(epoch),est_lambda_b_LS(epoch),est_h_b_LS(epoch),dummy] = pv_ECEF_to_NED(est_p_eb_ecef_LS(:,epoch),[0;0;0]);
    [est_phi_b_WLS( epoch ),est_lambda_b_WLS(epoch), est_h_b_WLS( epoch ) ,dummy] = pv_ECEF_to_NED( est_p_eb_ecef_WLS ( : , epoch ) , [ 0 ; 0 ; 0 ] ) ;
    [est_phi_b_WLSA( epoch ),est_lambda_b_WLSA(epoch), est_h_b_WLSA( epoch ) ,dummy] = pv_ECEF_to_NED( est_p_eb_ecef_WLSA ( : , epoch ) , [ 0 ; 0 ; 0 ] ) ;
    [est_phi_b_KF( epoch ),est_lambda_b_KF(epoch), est_h_b_KF( epoch ) ,dummy] = pv_ECEF_to_NED( est_p_eb_ecef_KF ( : , epoch ) , [ 0 ; 0 ; 0 ] ) ;
    [real_phi_b( epoch ),real_lambda_b(epoch), real_h_b( epoch ) ,dummy] = pv_ECEF_to_NED(true_p_eb_ecef_vec (: , epoch), [0;0;0]);
    time = time + GNSS_config.sampling;
end

%% Compute error statistics (RMSE)
error_ecef_LS = sqrt(sum(mean((est_p_eb_ecef_LS(:,1:100) - true_p_eb_ecef_vec(:,1:100)).^2,2)))
error_ecef_WLS = sqrt(sum(mean((est_p_eb_ecef_WLS(:,1:100) - true_p_eb_ecef_vec(:,1:100)).^2,2)))
error_ecef_WLSA = sqrt(sum(mean((est_p_eb_ecef_WLSA(:,1:100) - true_p_eb_ecef_vec(:,1:100)).^2,2)))
error_ecef_KF = sqrt(sum(mean((est_p_eb_ecef_KF(:,1:100) - true_p_eb_ecef_vec(:,1:100)).^2,2)))
KF_error = est_p_eb_ecef_KF - true_p_eb_ecef_vec;
WLS_error = est_p_eb_ecef_WLS - true_p_eb_ecef_vec;
LS_error = est_p_eb_ecef_LS - true_p_eb_ecef_vec;
t = 1:GNSS_config.no_epochs;
%%Display number of iterations (LS)
num_iter_LS
num_iter_WLS
num_iter_WLSA
%% Plot figures
t_vec = (0:GNSS_config.no_epochs-1)*GNSS_config.sampling;

% all pseudoranges
figure,
plot(t_vec,pseudorange(1:no_GNSS_meas,:).'), grid
    xlabel('time[s]'), ylabel('Pseudorange [m]')

% pick one satellite    
figure,
plot(t_vec,pseudorange(1,:).'), grid
    xlabel('time[s]'), ylabel('Pseudorange [m]')

% pseudorange(1,1) - pseudorange(1,end)    % pseudorange difference in no_epochs 
    
% elevations for all satellites    
figure,
plot(t_vec,rad2deg*elevations(1:no_GNSS_meas,:).'), grid
    xlabel('time[s]'), ylabel('elevations [m]')
    
% elevations for one satellites    
figure,
plot(t_vec,rad2deg*elevations(1,:).'), grid
    xlabel('time[s]'), ylabel('elevations [m]')

% true and estimated locations (LS)    
figure, plot(rad2deg*est_lambda_b_LS, rad2deg*est_phi_b_LS,'x','MarkerSize',5), hold on
    plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20), grid
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    axis([true_lambda_b_deg-0.001 true_lambda_b_deg+0.001 true_phi_b_deg-0.001 true_phi_b_deg+0.001])
%     plot_google_map    
    legend('LS estimates','True locations')
    
    % save KML file for processing in Google Earth
    kmlwritepoint('LS_solution', [true_phi_b_deg rad2deg*est_phi_b_LS], ...
        [true_lambda_b_deg rad2deg*est_lambda_b_LS]);
    hold off;
%plot GDOP
figure, plot(t_vec,GDOP_LS), grid
    xlabel('time[s]'), ylabel('GDOP')

% true and estimated locations (WLS)    
figure, plot(rad2deg*est_lambda_b_WLS, rad2deg*est_phi_b_WLS,'x','MarkerSize',5), hold on
    plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20), grid
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    axis([true_lambda_b_deg-0.001 true_lambda_b_deg+0.001 true_phi_b_deg-0.001 true_phi_b_deg+0.001])
%     plot_google_map    
    legend('WLS estimates','True locations')
    
    % save KML file for processing in Google Earth
    kmlwritepoint('WLS_solution', [true_phi_b_deg rad2deg*est_phi_b_WLS], ...
        [true_lambda_b_deg rad2deg*est_lambda_b_WLS]);
    hold off;
% true and estimated locations (WLSA)    
figure, plot(rad2deg*est_lambda_b_WLSA, rad2deg*est_phi_b_WLSA,'x','MarkerSize',5), hold on
    plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20), grid
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    axis([true_lambda_b_deg-0.001 true_lambda_b_deg+0.001 true_phi_b_deg-0.001 true_phi_b_deg+0.001])
%     plot_google_map    
    legend('WLSA estimates','True locations')
    
    % save KML file for processing in Google Earth
    kmlwritepoint('WLSA_solution', [true_phi_b_deg rad2deg*est_phi_b_WLSA], ...
        [true_lambda_b_deg rad2deg*est_lambda_b_WLSA]);
    hold off;
  % true and estimated locations (KF)    
figure, plot(rad2deg*est_lambda_b_KF, rad2deg*est_phi_b_KF,'x','MarkerSize',5), hold on
    plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20), grid
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    axis([true_lambda_b_deg-0.001 true_lambda_b_deg+0.001 true_phi_b_deg-0.001 true_phi_b_deg+0.001])
%     plot_google_map    
    legend('KF estimates','True locations')
    
    % save KML file for processing in Google Earth
    kmlwritepoint('KF_solution', [true_phi_b_deg rad2deg*est_phi_b_KF], ...
        [true_lambda_b_deg rad2deg*est_lambda_b_KF]);
    hold off;
    figure, plot(t_vec, KF_error(1,1:GNSS_config.no_epochs), 'r'), hold on 
    plot(t_vec, KF_error(2,1:GNSS_config.no_epochs), 'b')
    plot(t_vec, KF_error(3,1:GNSS_config.no_epochs), 'g')
    xlabel('time [s]'), ylabel('error [m]')
    hold off
    %plot the x, y, z covariance values over time
    figure, plot(t_vec, KF_cov(:,1), 'r'), hold on 
    plot(t_vec, KF_cov(:,2), 'b')
    plot(t_vec, KF_cov(:,3), 'g')
    xlabel('time [s]'), ylabel('covariance [m]')
    title('Covariance Values');
    axis([0,30,0,25])
    hold off
 KF_errormag = abs(sqrt(diag(KF_error'*KF_error)));
 WLS_errormag = abs(sqrt(diag(WLS_error'*WLS_error)));
 LS_errormag = abs(sqrt(diag(LS_error'*KF_error)));
  figure, plot(t_vec, KF_errormag(1:GNSS_config.no_epochs), 'r'), hold on 
    plot(t_vec, WLS_errormag(1:GNSS_config.no_epochs), 'b')
    plot(t_vec, LS_errormag(1:GNSS_config.no_epochs), 'g')
    xlabel('time [s]'), ylabel('error [m]')
    title('Error Magnitude Comparison');
    legend('KF', 'WLS', 'LS');
    hold off
    
  figure, plot(rad2deg*real_lambda_b, rad2deg*real_phi_b, 'r'), hold on 
  plot(rad2deg*est_lambda_b_KF, rad2deg*est_phi_b_KF,'x','MarkerSize',5)  
  xlabel('deg'), ylabel('deg')
    title('Actual Position');
    hold off
    
kmlwriteline('dynamicpos.kml', rad2deg*real_phi_b, rad2deg*real_lambda_b, 'Color', 'r','LineWidth',5)    % save KML file for processing in Google Earth
  