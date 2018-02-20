%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% ASSEMBLING THROUGH NO-DATA VALUES (STEP 5)
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

savedata = 1;
time = '';


%% Parameters for computation
BasicArea = 30*30; % Pixel size [m^2]
qAs = 0.98; % Quantile of dams' hedge DEM pixels used for assembling


%% Source files
dir_region = '.\region_mask\';
% Study area
study_region = imread([dir_region,'mask.tif']);
Iclass = geotiffinfo([dir_region,'mask.tif']);

fprintf('-> Loading Res_tot, h_res, stats_res & r2\n')
load('statistics','Res_tot')
load(['statistics_2',time],'h_res','stats_res','r2')

addpath('.\functions\step_5\')

% Landsat files
dir_landsat = '.\Landsat_results\';
listing = dir([dir_landsat,'\L*']);
L = length(listing);


%% Assembling through nodata values
img2 = struct('stats',cell(1,L),'pct',cell(1,L),'yd',cell(1,L));
fprintf('\n-> Filling reservoirs no-data values\n')
parfor l=1:L
    img2(l) = assemble(img2(l),savedata,time,listing(l).name,dir_landsat, ...
        study_region,Iclass,Res_tot,h_res,BasicArea,qAs)
end
clear l yrb dir_landsat
if savedata==1
    save(['statistics_2',time],'h_res','stats_res','img2','r2')
end

% ========================================================================