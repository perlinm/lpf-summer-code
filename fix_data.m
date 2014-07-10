

num_sets = 6;

for set = 1:num_sets
    disp(set)
    
    save_files = {sprintf('data/all-%i.mat',set),...
                  sprintf('data/calc-%i.mat',set)};
    
    for file = 1:length(save_files)
        clearvars -except set save_files file
        
        load(save_files{file});
        
        % stuff here
        
        save(save_files{file});
    end
end