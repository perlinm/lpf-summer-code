% this script uses simulated thruster inputs and test mass response
% data to extract information about the space craft

%% clear workspace

clearvars -except set;

%% load data

tic; fprintf('\nloading data (set %i)... ',set);
% list data files
data_dir = '/home/perlinm/lpf/code/data';
data_files = {'2014-06-30_planar.mat',...
              '2014-06-30_planar_destiff.mat',...
              '2014-06-30_planar_destiff_decapact.mat',...
              '2014-06-30_planar_destiff_decapact_tweak.mat',...
              '2014-07-02_3D.mat',...
              '2014-07-02_3D_diag.mat'};
% files which store parameters in "pDefault" instead of "pKey"
p_def_files = [3 5];
% files which only contain 2D data in the x-y plane
two_dim_files = 1:4;

% load data file, making sure that the choice of file is valid
if ~exist('set','var') || set > 6
    error('choose a valid data set');
else
    load(sprintf('%s/%s',data_dir,data_files{set}));
end
if find(p_def_files == set)
    pKey = pDefault;
end

% extract input parameters
for t = 1:length(TC)
    % thruster calibration
    calib_inp(t) = pKey.find(sprintf('CMNT_T%i_CAL',t)).y;
    % location of thruster
    r_th_inp(t,:) = ...
        ao_vec([pKey.find(sprintf('CMNT_T%i_POS_X',t)).y,...
                pKey.find(sprintf('CMNT_T%i_POS_Y',t)).y,...
                pKey.find(sprintf('CMNT_T%i_POS_Z',t)).y],'m');
    % thruster pitch
    alpha_inp(t) = pKey.find(sprintf('CMNT_T%i_ALPHA',t)).y;
    % x-y angle of thruster
    beta_inp(t) = pKey.find(sprintf('CMNT_T%i_BETA',t)).y;
end

% fake 3D data, if necessary
if find(two_dim_files == set)
    tm1z = 0*tm1x;
    tm1theta = 0*tm1phi;
    tm1eta = 0*tm1phi;
    tm2z = 0*tm2x;
    tm2theta = 0*tm2phi;
    tm2eta = 0*tm2phi;
    % correct for thruster pitch if its projection is not in the x-y
    %   plane
    for t = 1:length(TC)
        TC(t) = TC(t)*cos(alpha_inp(t));
    end
end

% data vectors
accel_data = [tm1x,tm1y,tm1z,tm1theta,tm1eta,tm1phi,...
              tm2x,tm2y,tm2z,tm2theta,tm2eta,tm2phi];
input_force = TC; % input command vector
% some derived quantities
num_sets = length(accel_data); % number of response data arrays
num_masses = floor(num_sets/6); % number of masses
num_thrusters = length(input_force); % number of thrusters
time = toc; fprintf('%.0f sec\n',time);

%% define global variables

% location of center of mass
r_com_mechanical = ao_vec([0.005,0.006,0.47],'m');
% locations of test masses
r_tm_mechanical = [ao_vec([0.188,0,0.6093],'m');...
                   ao_vec([-0.188,0,0.6093],'m')];
% locations of test masses relative to center of mass
r_tm = r_tm_mechanical - repmat(r_com_mechanical,num_masses,1);

% mass and moment of inertia of space sraft
m_sc = ao(422.7,plist('yunits','kg'));
I_sc_xx = ao(202.5,plist('yunits','kg m^2'));
I_sc_yy = ao(209.7,plist('yunits','kg m^2'));
I_sc_zz = ao(191.7,plist('yunits','kg m^2'));
radian = ao(1,plist('yunits','rad'));

x_d = @(m) 6*m-5; % for indexing accel_data
y_d = @(m) 6*m-4;
z_d = @(m) 6*m-3;
theta_d = @(m) 6*m-2;
eta_d = @(m) 6*m-1;
phi_d = @(m) 6*m;
x_c = 1; y_c = 2; z_c = 3; % for indexing dimensions

% filter parameters
low_pass_cap = 3e-2; % Hz
low_pass_order = 4;
notch_function = 'BH92'; % window function used for notch filter
notch_order = 2048;
downsample_frequency = 0.1; % Hz

% misc
sample_frequency = input_force(1).fs; % Hz
signal_range = [1e-3 1e-1]; % Hz
input_bandwidth = 1e-4; % Hz; fixme: find bandwidth automatically
ampl_window = 'HFT248D';
histogram_bins = 20;
fig_dir = './figures';

% convenient plists
unit_pl = plist('exceptions','');
pl = plist('legends','off');
time_pl = plist(...
    'legends','off',...
    'xranges',[3 3.1]*1e4);
freq_pl = plist(...
    'legends','off',...
    'xscales',{'all','lin'},...
    'xranges',[6.4,8.2]*1e-3,...
    'complexplottype','absrad');

%% identify thruster signals

tic; fprintf('identifying thruster signals... ');
for t = 1:num_thrusters
    % amplitude spectrum of signal
    input_force_spec(t) = split(psd(input_force(t),...
        plist('win',ampl_window,'scale','AS')),...
        plist('frequencies',signal_range));
    % magnitude of input signal
    input_mag(t) = max(input_force_spec(t))*sqrt(2);
    % time-series index of maximal magnitude point
    [~,max_index] = max(input_force_spec(t).y);
    % frequency of signal
    frequency(t) = max(input_force_spec(t).select(max_index),...
        plist('axis','x'));
end; clear max_index;
time = toc; fprintf('%.0f sec\n',time);

%% make data filters

tic; fprintf('making data filters... ');
% make low pass filter to eliminate high frequency noise
prefilter = miir(plist(...
  'name','lowpass pre-filter',...
  'type','lowpass',...
  'fc',low_pass_cap,...
  'order',low_pass_order,...
  'gain',1,...
  'fs',sample_frequency));
% make band pass filter to isolate signals by frequency
for t = 1:num_thrusters
    notch_filter(t) = mfir(plist(...
        'name','notch filter',...
        'type','bandpass',...
        'fc',(frequency(t).y + ...
              [-input_bandwidth/2,input_bandwidth/2]),...
        'order',notch_order,...
        'win',notch_function,...
        'gain',1,...
        'fs',downsample_frequency));
end
time = toc; fprintf('%.0f sec\n',time);

%% filter data

tic; fprintf('filtering data...\n');
for t = 1:num_thrusters
    fprintf('  thruster %i\n',t);
    for s = 1:num_sets
        % filter data by thruster frequencies
        accel(s,t) = signal_filter(...
            accel_data(s),prefilter,notch_filter(t));
    end
    % filter input signals to match acceleration data
   input_wave(t) = signal_filter(...
        input_force(t),prefilter,notch_filter(t));
end
time = toc; fprintf('    %.0f sec\n',time);

%% characterize filtering effects

tic; fprintf('characterizing filtering effects... ');
for t = 1:num_thrusters
    % construct unfiltered input signal
    input_wave_nofilter(t) = input_force(t)-mean(input_force(t));
    input_wave_nofilter(t).downsample(...
        plist('factor',sample_frequency/downsample_frequency));
    input_wave_nofilter(t).select(...
        notch_order:len(input_wave_nofilter(t)));
    % compute transfer function of filter
    filter_transfer(t) = tfe(input_wave_nofilter(t),input_wave(t));
    % transfer function of filter at signal frequency
    filter_sig_transfer(t) = interp(filter_transfer(t),...
        plist('vertices',frequency(t).y));
    % find actual and expected signal phase shift caused by filter
    actual_shift(t) = angle(filter_sig_transfer(t));
    expected_shift(t) = ...
        angle(prefilter.resp(frequency(t).y)) + ...
        angle(notch_filter(t).resp(frequency(t).y));
end
time = toc; fprintf('%.0f sec\n',time);

%% find thruster time delays

tic; fprintf('finding thruster time delays... ');
for s = 1:num_sets
    for t = 1:num_thrusters
        % find transfer function of input to acceleration
        transfer(s,t) = tfe(input_wave(t),accel(s,t));
        % transfer function at signal frequency
        sig_transfer(s,t) = interp(transfer(s,t),...
            plist('vertices',frequency(t).y));
        % find phase difference between acceleration and input
        phase(s,t) = abs(angle(sig_transfer(s,t)));
        % use phase to find accel_com_dir: the direction of
        %   acceleration of the space craft when the thrusters provide
        %   a positive force
        if phase(s,t) < pi/2
            accel_com_dir(s,t) = 1;
        else
            accel_com_dir(s,t) = -1;
            phase(s,t) = abs(phase(s,t) - pi);
        end
        % find time delay from phase
        delay_all(s,t) = phase(s,t).y/(2*pi*frequency(t).y);
    end
end
for t = 1:num_thrusters
    % find average delay for each thruster
    delay(t) = average(delay_all(:,t));
end
time = toc; fprintf('%.0f sec\n',time);

%% find thruster calibrations and orientations

tic; fprintf('finding thruster calibrations and orientations... ');
for t = 1:num_thrusters
    % calibration (calib), pitch (alpha), and x-y angle (beta) of
    %   thrusters, averaged over results from all test masses
    calib(t) = 0;
    alpha(t) = 0;
    beta(t) = 0;
    for m = 1:num_masses
        % time series acceleration and velocity vectors
        accel_tm(t,m,:) = ...
            [accel(x_d(m),t),accel(y_d(m),t),accel(z_d(m),t)];
        angular_accel_tm(t,m,:) = ...
            [accel(theta_d(m),t),accel(eta_d(m),t),accel(phi_d(m),t)];
        % Euler acceleration
        accel_euler(t,m,:) = ...
            ao_cross(angular_accel_tm(t,m,:),r_tm(m,:));
        for i = 1:3
            % acceleration of center of mass of space craft
            % effects not accounted for:
            %   changing r_tm in time (error unknown)
            %   coriolis force (1e-6 fractional error in 2D data)
            %   centrifugal force (1e-8 fractional error in 2D data)
            accel_com(t,m,i) = accel_tm(t,m,i) - accel_euler(t,m,i);
        end
        % covariance matrix of normalized response
        accel_com_cov(t,m) = cov(accel_com(t,m,x_c),...
                                 accel_com(t,m,y_c),...
                                 accel_com(t,m,z_c));
        % eigenvectors of cov_resp and diagonzliaed cov_resp
        [evecs(t,m,:,:),cov_diag(t,m,:,:)] = ...
            eigs(accel_com_cov(t,m).y);
        % thruster calibration, pitch, and x-y angle
        calib_all(t,m) = ...
            m_sc.y*sqrt(2*cov_diag(t,m,1,1))/input_mag(t).y;
        alpha_all(t,m) = ...
            atan2(-accel_com_dir(z_d(m),t)*abs(evecs(t,m,z_c,1)),...
            sqrt(evecs(t,m,x_c,1)^2 + evecs(t,m,y_c,1)^2));
        beta_all(t,m) = ...
            atan2(-accel_com_dir(y_d(m),t)*abs(evecs(t,m,y_c,1)),...
                  -accel_com_dir(x_d(m),t)*abs(evecs(t,m,x_c,1)));
        calib(t) = calib(t) + calib_all(t,m)/num_masses;
        alpha(t) = alpha(t) + alpha_all(t,m)/num_masses;
        beta(t) = beta(t) + beta_all(t,m)/num_masses;
    end
end
% error in calibration, pitch, and x-y angle
calib_error = calib-calib_inp.y;
alpha_error = alpha-alpha_inp.y;
beta_error = beta-beta_inp.y;
time = toc; fprintf('%.0f sec\n',time);

%% find center of mass responses in the respective eigenbases

tic; fprintf('finding com responses in the respective eigenbases... ');
for t = 1:num_thrusters
    for m = 1:num_masses
        % accelertions of com projected onto respective eigenvectors
        com_resp_proj(t,m,:) = m_sc/input_mag(t)*...
            in_basis(evecs(t,m,:,:),accel_com(t,m,:));
    end
end
time = toc; fprintf('%.0f sec\n',time);

%% save data and calculations

tic; fprintf('saving results... ');
save(sprintf('%s/all-%i.mat',data_dir,set));
save(sprintf('%s/calc-%i.mat',data_dir,set),...
        'num_sets','num_thrusters','num_masses',...
        'delay_all','delay',...
        'calib_inp','calib_all','calib','calib_error',...
        'alpha_inp','alpha_all','alpha','alpha_error',...
        'beta_inp','beta_all','beta','beta_error');
data_tables;
time = toc; fprintf('%.0f sec\n',time);