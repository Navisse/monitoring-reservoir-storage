%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% SECOND DETECTION OF RESERVOIRS AFTER MASK (STEP 3)
% 
% N. Avisse (nicolas.avisse@gmail.com)
% (10/15/2017)
%
% Please cite the following paper:
% Avisse, N., Tilmant, A., Müller, M. F., and Zhang, H.: Monitoring
% small reservoirs storage from satellite remote sensing in inaccessible
% areas, Hydrology and Earth System Sciences Discussions, 2017, 1 – 23,
% doi:/10:5194/hess-2017-373, in review, 2017.
% https://www.hydrol-earth-syst-sci-discuss.net/hess-2017-373/
%
%=========================================================================
%%

clear variables
close all

savedata = 1; % Save data & results?


%% Parameters for computation
aMin = 20; % Minimum quantity of isolated points to consider it is a reservoir
p_prctile = 100; % percentile for thresholds calculation


%% Source files
dir_landsat = '.\Landsat_results\';

load('statistics','Res_tot')
Res_tot = uint8(Res_tot);

% Landsat files
listing = dir([dir_landsat,'\L*']);
L = length(listing);


%% Analysis for each date acquired

fprintf('Detecting reservoirs\n0%%')
progress1 = L/10;
progress2 = L/40;
for l=1:L
    % progression bar
    if l==L
        fprintf('100%%\n')
    elseif l>progress1
        fprintf('%.0f%%',progress1/L*100)
        progress1 = progress1+L/10;
        progress2 = L/40;
    elseif l>progress1-L/10+progress2
        fprintf('.')
        progress2 = progress2+L/40;
    end
    
    load([dir_landsat,listing(l).name,'\data.mat'],'d')
    d.Res_i = [];

    thrshld_NDVI = max(d.NDVI((Res_tot==1).*(d.mask_i==1)==1));
    % if water has not been detected by Fmask
    if isempty(thrshld_NDVI)
        thrshld_NDVI = -0.1;
    end
    W_ndvi = uint8(d.NDVI<=thrshld_NDVI);
    W_ndvi = W_ndvi.*Res_tot.*uint8((d.mask_i~=255));
    % dynamic definition of a MNDWI threshold
    if thrshld_NDVI==-0.1
        thrshld_MNDWI = prctile(d.MNDWI((Res_tot==1).*(W_ndvi==1)==1),p_prctile);
        if isnan(thrshld_MNDWI)
            thrshld_MNDWI = prctile(d.MNDWI((Res_tot==1).*(d.mask_i==4)==1),70);
        end
    else
        thrshld_MNDWI = prctile(d.MNDWI((Res_tot==1).*(d.mask_i==1).*(W_ndvi==1)==1),p_prctile);
    end
    W_mndwi = uint8(d.MNDWI<=thrshld_MNDWI);
    W_mndwi = W_mndwi.*Res_tot.*uint8(d.mask_i~=255);
    [~,Lr,Nr] = bwboundaries(W_mndwi==1,'noholes');
    % removing noise
    for n=1:Nr
        id = find(Lr==n);
        if length(id)<=aMin
            W_mndwi(id) = 0;
        end
    end
    clear n Nr Lr id

    d.Res_i = W_mndwi;
    d.thrshld_MNDWI = thrshld_MNDWI;

    if savedata==1
        save([dir_landsat,listing(l).name,'\data.mat'],'d')
    end
    clear d W_*
end
clear progress*

% ========================================================================