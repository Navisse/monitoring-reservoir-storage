%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% CALCULATING STATISTICS ON FLOODED AREAS FOR EACH RESERVOIR (STEP 4)
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
% REQUIRED for this analysis:
% DEM files
%
%=========================================================================
%%

clear variables
close all

% first: method 1, then 2
method = 2; % Method & DEM chosen? 1: original, 2: statistical correction
graph = 1;  % Plot curves?
savedata = 0;% Save data & results?

% possibility to create updated filling curves after step 6
time = ''; % '' or 't2'


%% Parameters for computation
BasicArea = 30*30;% Pixel size [m^2]

% index of reservoirs to associate to specific reservoirs
id_srtmx = [];
id_srtmc = [];

% index of reservoirs to see filling curves
id_plot = [];

% % Degree of polynomial to fit (index of dams Excel file)
n_fit = [0   1   NaN]; % 0 for local pol fit
% First value to use to fit (index of dams Excel file)
x_fit = [0  .35 0]; % NaN for r fit
% Initial H value if constrained fitting
y0_fit = [NaN NaN 730];
% span for rloess fitting
span_fit = [0.3 NaN NaN];


%% Source files
% dams
fprintf('-> Loading data\n')
dir_region = '.\region_mask\';
Iclass = geotiffinfo([dir_region,'mask.tif']);
[~,~,dams] = xlsread([dir_region,'dams_info.xlsx'],'A2:C3');

load('statistics.mat','Res_tot_fmask','Res_tot')
if method==2 || graph==1
    load(['statistics',time],'r')
    % statistics_2 files pre-created considering each DEM only
    r2ASTER = load(['statistics_2ASTER',time],'r2');
    r2SRTMC = load(['statistics_X2SRTM-C',time],'r2');
    r2SRTMX = load(['statistics_X2SRTM-X',time],'r2');
    scrsz = get(groot,'ScreenSize');
end

% DEM files
dir_aster = '.\DEM\ASTER\';
dir_srtm_c = '.\DEM\SRTM_1_Arc-Second_Global\';
dir_srtm_x = '.\DEM\SRTM_X-SAR\';
dir_landsat = '.\Landsat_results\';
addpath('.\functions\step_4\')

[ndd,~] = size(dams);% nb of dams according to studies
[~,Lt,Nt] = bwboundaries(Res_tot,'noholes');

% dam number in the xls file for each reservoir in the image
IMtoDAM = zeros(1,Nt);
% reservoir number on the image for each dam
DAMtoIM = zeros(1,ndd);
for nd=1:ndd
    if ~isnan(dams{nd,2})% if the location of the dam is known
        l = Lt(round((Iclass.BoundingBox(2,2)-dams{nd,3})/Iclass.PixelScale(2)), ...
            round((dams{nd,2}-Iclass.BoundingBox(1,1))/Iclass.PixelScale(1)));
        if l~=0% if dam has been detected by remote sensing
            IMtoDAM(l) = nd;
            DAMtoIM(nd) = l;
        end
        clear l
    end
end
clear nd

% Landsat files
listing = dir([dir_landsat,'\L*']);
L = length(listing);

% converting DEMs in Landsat original CRS
% ASTER GDEM
dos(['gdalwarp -overwrite -s_srs EPSG:4326', ...
    ' -t_srs EPSG:326',num2str(Iclass.Zone),' -r near -te ',...
    num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
    num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
    ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
    ' -of GTiff ',dir_aster,'DEM.tif dem_test.tif']);
DEM_ASTER = imread('dem_test.tif');
DEM_ASTER = double(DEM_ASTER);
delete dem_test.tif
% SRTM Void Filled 1 arc second
dos(['gdalwarp -overwrite -s_srs EPSG:4326', ...
    ' -t_srs EPSG:326',num2str(Iclass.Zone),' -r near -te ',...
    num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
    num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
    ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
    ' -of GTiff ',dir_srtm_c,'SRTM_VF_1arc.tif dem_test.tif']);
DEM_SRTM_C = imread('dem_test.tif');
DEM_SRTM_C = double(DEM_SRTM_C);
delete dem_test.tif
% SRTM X-SAR
dos(['gdalwarp -overwrite -s_srs EPSG:4326', ...
    ' -t_srs EPSG:326',num2str(Iclass.Zone),' -r near -te ',...
    num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
    num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
    ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
    ' -of GTiff ',dir_srtm_x,'srtm_x-sar_dem.tif dem_test.tif']);
DEM_SRTM_X = imread('dem_test.tif');
DEM_SRTM_X = double(DEM_SRTM_X);
delete dem_test.tif
clear dir_aster dir_srtm*


%% Calculating statistics on flooded areas of each reservoir (Step 3)
if method==1
    r = struct('flooded',cell(1,Nt),'theory',cell(1,Nt));
else
    r2 = struct('stats',cell(1,Nt),'theory',cell(1,Nt),'DS',cell(1,Nt),'pct',cell(1,Nt));
    % considering most often flooded pixels as a reference for DEM correcting
    % & nodata filling
    h_res = zeros(Iclass.Height,Iclass.Width);
    stats_res = zeros(Iclass.Height,Iclass.Width);
    % max storage detected
    data_to_write = NaN(ndd,1);
    i_end = NaN(Nt,1);
end

fprintf('\n-> Building filling curves\n')
idx_p = 0;
for nt=1:Nt%DAMtoIM(id_plot) % for each reservoir
    if method==1
        if any(nt==DAMtoIM(id_srtmx))
            DEM = DEM_SRTM_X;
        elseif any(nt==DAMtoIM(id_srtmc))
            DEM = DEM_SRTM_C;
        else
            DEM = DEM_ASTER;
        end
        r(nt) = calculate_theory(r(nt),L,dams,nt,Nt,Lt,IMtoDAM(nt),DEM,BasicArea);
        clear DEM
    elseif nnz([0 0 0 ones(1,L)].*(sum(r(nt).flooded==1)~=0))>=10
        [r2(nt).theory, h_res_nt, stats_res_nt, i_end(nt), idx_p] = calculate_theory2(r(nt), ...
            r2ASTER.r2(nt),r2SRTMC.r2(nt),r2SRTMX.r2(nt),id_srtmc,id_srtmx,idx_p, ...
            savedata,graph,scrsz,L,n_fit,x_fit,y0_fit,span_fit, ...
            dams,nt,Nt,Lt,IMtoDAM(nt),DAMtoIM,Iclass,BasicArea);
        stats_res = stats_res + stats_res_nt;
        h_res = h_res + h_res_nt;
    end
end
if method~=1 && savedata==1
    for nt=1:Nt
        if nnz([0 0 0 ones(1,L)].*(sum(r(nt).flooded==1)~=0))>=10
            data_to_write(IMtoDAM(nt)) = r2(nt).theory.V2f(i_end(nt));
            % saving volume-elevation-area into an Excel file
            xlswrite([dir_region,'dams_info',time,'.xlsx'],[r2(nt).theory.A2/1e6 r2(nt).theory.V2f r2(nt).theory.H2f], ...
                dams{IMtoDAM(nt),1},['A2:C',num2str(i_end(nt)+1)])
        end
    end
end
clear nt i_end

if method==1
    fprintf('\n-> Collecting statistics on flooded pixels\n 0%%')
    progress1 = L/10;
    progress2 = L/40;
    for l=1:L% for each image
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
        
        for nt=1:Nt
            if any(nt==DAMtoIM(id_srtmx))
                DEM = DEM_SRTM_X;
            elseif any(nt==DAMtoIM(id_srtmc))
                DEM = DEM_SRTM_C;
            else
                DEM = DEM_ASTER;
            end
            if isempty(time)
                r(nt).flooded = calculate_flooded(r(nt).flooded,[],l,Lt,nt, ...
                    DEM,d.mask_i,d.Res_i,[]);
            else
                r(nt).flooded = calculate_flooded(r(nt).flooded,time,l,Lt,nt, ...
                    DEM,d.mask_i,d.Res_i,d.Res_f);
            end
            clear DEM
        end
        clear nt d
    end
    clear l progress*
end

if savedata==1
    fprintf('\n-> Saving data\n')
    if method==1
        save(['statistics',time],'Res_tot_fmask','Res_tot','r')
    else
        clear r
        save(['statistics_2',time],'h_res','stats_res','r2')
        xlswrite([dir_region,'dams_info',time,'.xlsx'],data_to_write,'Data','D2:D5')
    end
end

% ========================================================================