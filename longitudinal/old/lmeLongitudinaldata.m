function [lme_linear, lme_quad, data_table] = lmeCalc(sid, long_var, score, dummyon)
% Furnction: Calculates linear mixed effects on longitudinal data


s = unique(sid);
% % Center each individual's reading scores
% for ii = 1:length(s)
%    index = find(strcmp(s(ii),sid));
%    total = 0;
%    for jj = 1:length(index)
%        total = plus(total, score(index(jj))); 
%    end
%    avg = total/length(index);
%    
%    for kk = 1:length(index);
%        score_sq_unique(index(kk), 1) = score(index(kk)) - avg;
%    end
% end
% score_adj = score_sq_unique;
if dummyon == 0
    % Centering of time course variable
    for ii = 1:length(s)
        index = find(strcmp(s(ii),sid));
        total = 0;
        for jj = 1:length(index)
            total = plus(total, long_var(index(jj)));
        end
        avg = total/length(index);
        for kk = 1:length(index);
            time_adj(index(kk), 1) = long_var(index(kk)) - avg;
        end
    end
    uncentered = long_var;
    long_var = time_adj;
    % Create squared hours variable to use in quadratic model
    long_var_sq = long_var.^2;
    % Create DataSet
    data_table = dataset(sid, uncentered, long_var, long_var_sq, score);
    % Calculate LME fit
    % Make sid a categorical variable
    data_table.sid = categorical(data_table.sid);
    % Fit the model on the uncentered data as changing linearly
    lme_linear = fitlme(data_table, 'score ~ long_var + (1|sid)');
    % Fit the model on uncentered data as changing quadratically
    lme_quad = fitlme(data_table, 'score ~ long_var + long_var_sq + (1|sid)');
    
    stats(ii).lme_quad = lme_quad;

elseif dummyon == 1
    % make the longitudinal variable a categorical variable
    long_var = categorical(long_var);
    % Create DataSet
    data_table = dataset(sid, long_var, score);
    % Calculate LME fit
    % Make sid a categorical variable
    data_table.sid = categorical(data_table.sid);
    % Fit the model on the centered data as changing linearly
    lme_linear = fitlme(data_table, 'score ~ long_var + (long_var|sid)');
    
end
% Collate data into stats struct

stats(ii).test_name = test_name;  
stats(ii).lme_linear = lme_linear;
stats(ii).data_table = data_table;  



return


