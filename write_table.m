% this function generates a latex table
% requirements: cprotect package (in the latex document)

% columns, rows, and data must all be cell arrays containing
%   column names, row names, and data respectively as strings

% label is optional
function [] = write_table(file,caption,data,rows,columns,label)
    fprintf(file,'\n\\begin{table}[H]\n');
    fprintf(file,'\\centering\n');
    fprintf(file,'\\cprotect\\caption{%s}\n',caption);
    fprintf(file,'\\begin{tabular}{|%s} \\hline\n',...
                             repmat('c|',1,length(columns)+1));
    fprintf(file,'~');
    for c = 1:length(columns)
        fprintf(file,' & %s',columns{c});
    end
    fprintf(file,' \\\\ \\hline');
    for r = 1:length(rows)
        fprintf(file,'\n%s',rows{r});
        for c = 1:length(columns)
            fprintf(file,' & %s',data{r,c});
        end
        fprintf(file,' \\\\');
    end
    fprintf(file,' \\hline\n');
    fprintf(file,'\\end{tabular}\n');
    fprintf(file,'\\label{%s}\n',label);
    fprintf(file,'\\end{table}\n');
end