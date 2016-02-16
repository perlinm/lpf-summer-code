%% make latex table of results for postes

clear pest calib_num elev_num azim_num calib_in elev_in azim_in ...
  calib elev azim data rows cols;

key = pest('keys/ModAll3D_set1.mat');
calib_num = y(unpack(out('pca').getObjectAtIndex(1).search('calibration_estimates')));
elev_num = y(unpack(out('pca').getObjectAtIndex(1).search('elevation_estimates')));
azim_num = y(unpack(out('pca').getObjectAtIndex(1).search('azimuth_estimates')));
for i = 1:8
  calib_in{i} = sprintf('%.0f',key.find(['CMNT_T' num2str(i) '_CAL']).y);
  elev_in{i} = sprintf('%.2f',key.find(['CMNT_T' num2str(i) '_ALPHA']).y);
  azim_in{i} = sprintf('%.2f',key.find(['CMNT_T' num2str(i) '_BETA']).y);
  calib{i} = sprintf('%.2f',calib_num(i));
  elev{i} = sprintf('%.2f',elev_num(i));
  azim{i} = sprintf('%.2f',azim_num(i));
  
  data{i,1} = calib{i};
  data{i,2} = elev_in{i};
  data{i,3} = elev{i};
  data{i,4} = azim_in{i};
  data{i,5} = azim{i};
  rows{i} = ['T' num2str(i)];
end

cols = {...
  '$c_{\t{out}}$',...
  '$\alpha_{\t{in}}$','$\alpha_{\t{out}}$',...
  '$\beta_{\t{in}}$','$\beta_{\t{out}}$',...
  };

file = fopen('code/docs/poster/results_table.tex','w');
write_table(file,['Analysis results by thruster; '...
  'angles are given in radians and all input calibrations are 1'],...
  data,rows,cols,'results')
fclose(file);