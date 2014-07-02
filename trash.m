% velocity things

%% compute velocities

tic; fprintf('computing velocities... ');
for s = 1:data_sets
    vel_data(s) = indef_integral(accel_data(s));
end
time = toc; fprintf('%.0f sec\n',time);


%% filter data (exists)

for t = 1:num_thrusters
    for m = 1:num_masses
        vel(s,t) = signal_filter(...
            vel_data(s),prefilter,notch_filter(t));
    end
end


%% find thruster calibrations and orientations (exists)

for s = 1:data_sets
    for t = 1:num_thrusters
        vel_tm(t,m,:) = ...
            [vel(x(m),t),vel(y(m),t),vel(z(m),t)];
        angular_vel_tm(t,m,:) = ...
            [vel(theta(m),t),vel(eta(m),t),vel(phi(m),t)];
        
        accel_coriolis(t,m,:) = ...
            2*cross_ao(angular_vel_tm(t,m,:),vel_tm(t,m,:));
        accel_centrifugal(t,m,:) = ...
            cross_ao(angular_vel_tm(t,m,:),...
                     cross_ao(angular_vel_tm(t,m,:),r_tm(m,:)));
    end
end



