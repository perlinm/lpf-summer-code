clear all;
data_dir = '/home/perlinm/lpf/code/data';
data_files = {'2014-06-30_planar.mat',...
              '2014-06-30_planar_destiff.mat',...
              '2014-06-30_planar_destiff_decapact.mat',...
              '2014-06-30_planar_destiff_decapact_tweak.mat',...
              '2014-07-02_3D.mat',...
              '2014-07-02_3D_diag.mat'};

out_file = fopen('data_tables.tex','w');

fprintf(out_file,'\n\\section*{Calibration}\n');
for set = 1:length(data_files)
    load(sprintf('%s/calc-%i.mat',data_dir,set),...
                 'num_thrusters','num_masses',...
                 'calib_inp','calib_all','calib');

    caption = sprintf('%s\n\\verb|%s|',...
                      ['Calibration measurements for ' ...
                       'simulation ID:\\'],...
                      data_files{set}(1:end-4));

    rows = arrayfun(@(n)num2str(n,'T%i'),1:num_thrusters,'unif',0);
    first_columns{1} = 'Actual';
    first_columns{2} = 'Mean';
    data_columns = arrayfun(@(n)num2str(n,'TM%i'),...
                            1:num_masses,'unif',0);
    columns = [first_columns,data_columns];

    data = [calib_inp',calib',calib_all];
    data = arrayfun(@(n)num2str(n,'%.2f'),data,'unif',0);

    label = sprintf('calib-%i',set);

    write_table(out_file,caption,data,rows,columns,label);
end
fprintf(out_file,'\\newpage');

clearvars -except data_dir data_files out_file

fprintf(out_file,'\n\\section*{Time delay}\n\n');
for set = 1:length(data_files)
    load(sprintf('%s/calc-%i.mat',data_dir,set),...
                 'num_sets','num_thrusters',...
                 'delay_all','delay');

    if strfind(data_files{set},'planar')
        num_sets = 6;
        da = delay_all([1 2 6 7 8 12],:);
        clear delay_all;
        delay_all = da;
        clear da;
    end

    caption = sprintf('%s\n\\verb|%s|',...
                      ['Thruster time delays in seconds for ',...
                       'simulation ID:\\'],...
                      data_files{set}(1:end-4));

    rows = arrayfun(@(n)num2str(n,'T%i'),1:num_thrusters,'unif',0);
    first_columns{1} = 'Actual';
    first_columns{2} = 'Mean';
    data_columns = arrayfun(@(n)num2str(n,'A%i'),...
                            1:num_sets,'unif',0);
    columns = [first_columns,data_columns];

    data = [repmat(0.02,num_thrusters,1),delay',delay_all'];
    data = arrayfun(@(n)num2str(n,'%.2f'),data,'unif',0);

    label = sprintf('delay-%i',set);

    write_table(out_file,caption,data,rows,columns,label);
end
fprintf(out_file,'\\newpage');

clearvars -except data_dir data_files out_file

fprintf(out_file,'\n\\section*{Thruster pitch}\n\n');
for set = 1:length(data_files)
    load(sprintf('%s/calc-%i.mat',data_dir,set),...
                 'num_thrusters','num_masses',...
                 'alpha_inp','alpha_all','alpha');

    if strfind(data_files{set},'planar')
       continue
    end

    caption = sprintf('%s\n\\verb|%s|',...
                      ['Error in mrad in thruster pitch ($\alpha$)'...
                       ' for simulation ID:\\'],...
                      data_files{set}(1:end-4));

    rows = arrayfun(@(n)num2str(n,'T%i'),1:num_thrusters,'unif',0);
    first_columns{1} = 'Actual';
    data_columns = arrayfun(@(n)num2str(n,'TM%i'),...
                            1:num_masses,'unif',0);
    columns = [first_columns,data_columns];

    data = 1e3*[(alpha-alpha_inp)',alpha_all-repmat(alpha_inp',1,2)];
    data = arrayfun(@(n)num2str(n,'%.2f'),data,'unif',0);

    label = sprintf('alpha-%i',set);

    write_table(out_file,caption,data,rows,columns,label);
end
fprintf(out_file,'\\newpage');

clearvars -except data_dir data_files out_file

fprintf(out_file,'\n\\section*{Thruster azimuth}\n\n');
for set = 1:length(data_files)
    load(sprintf('%s/calc-%i.mat',data_dir,set),...
                 'num_thrusters','num_masses',...
                 'beta_inp','beta_all','beta');

    caption = sprintf('%s\n\\verb|%s|',...
                      ['Error in mrad in thruster azimuth ($\beta$)'...
                       ' for simulation ID:\\'],...
                      data_files{set}(1:end-4));

    rows = arrayfun(@(n)num2str(n,'T%i'),1:num_thrusters,'unif',0);
    first_columns{1} = 'Actual';
    data_columns = arrayfun(@(n)num2str(n,'TM%i'),...
                            1:num_masses,'unif',0);
    columns = [first_columns,data_columns];

    data = 1e3*[(beta-beta_inp)',beta_all-repmat(beta_inp',1,2)];
    data = arrayfun(@(n)num2str(n,'%.2f'),data,'unif',0);

    label = sprintf('beta-%i',set);

    write_table(out_file,caption,data,rows,columns,label);
end
