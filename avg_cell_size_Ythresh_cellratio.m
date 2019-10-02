%% average cell size & fraction of cells eliminated due to Ythresh [9/18/2019]
% In response to Vicente's comment (Supplementary Info), include average size of cells and 
% fraction of cells eliminated due to YFP threshold at 50 a.u.

% load 14_aggdata_YFPnorm_t2avg.mat


%% average cell size

cell_size_ypos_avg_expt = [];
for expt = 5:10
     
    % average size of cells in an expt (first mean within an image, then
    % mean across images of all conditions of all time)
    cell_size_ypos_avg_expt(expt-4) = mean(mean(aggdata.cell_size_ypos_avg_xy{expt}));
    
    % first time point only
%    cell_size_ypos_avg_expt(expt-4) = mean(aggdata.cell_size_ypos_avg_xy{expt}(1,:));
    
end

% average cell size (across experiments)
mean(cell_size_ypos_avg_expt)
std(cell_size_ypos_avg_expt)


%% fraction of cells eliminated due to Ythresh
clearvars -except aggdata
for expt = 5:10
    % before YFP threshold
    total_before_ythresh(expt-4) = sum(sum(aggdata.cell_num{expt}));

    % after YFP threshold
    total_after_ythresh(expt-4) = sum(sum(aggdata.cell_num_ypos{expt}));

    cellnum_elim(expt-4) = total_before_ythresh(expt-4) - total_after_ythresh(expt-4);
    cellratio_elim(expt-4) = cellnum_elim(expt-4) / total_before_ythresh(expt-4);

end

min(cellratio_elim)
max(cellratio_elim)




