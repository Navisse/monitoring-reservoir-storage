%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FILTERING IMAGES (STEP 2)
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

savedata = 1; % Save data & results? (1: Yes)

%% Source files
fprintf('Loading data for region %s\n',region)
dir_landsat = '.\Landsat_results\';

tErrors = 5; % Threshold above which we consider water really appeared

dir_region = '.\region_mask\';
Iclass = geotiffinfo([dir_region,'mask.tif']);
% dams
[~,~,dams] = xlsread([dir_region,'dams_info.xlsx'],'A2:C3');
% study area
[study_region,R] = geotiffread([dir_region,'mask.tif']);
[ndd,~] = size(dams);% nb of dams according to studies

% Landsat files
listing = dir([dir_landsat,'\L*']);
L = length(listing);

Res_tot_fmask = zeros(Iclass.Height,Iclass.Width);
Res_tot = zeros(Iclass.Height,Iclass.Width);
realRes = zeros(Iclass.Height,Iclass.Width);


%% Filtering images
% Superimposing each date acquired
fprintf('Summing images\n0%%')
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
    
    load([dir_landsat,'\',listing(l).name,'\data.mat'],'d')
    
    idMask = find(d.mask_i.*study_region==1); % water body pixels location
    Res_tot_fmask(idMask) = Res_tot_fmask(idMask)+1; % sum

    clear d idMask
end
clear progress* l

% facultative step: if reservoirs are known
fprintf('Selecting real reservoirs\n0%%')
for nd=1:ndd
    % Selection of reservoirs that were detected by Orient
    realRes(round((Iclass.BoundingBox(2,2)-dams{nd,3})/Iclass.PixelScale(2)), ...
        round((dams{nd,2}-Iclass.BoundingBox(1,1))/Iclass.PixelScale(1))) = 1;
end
clear nd

% Selection of reservoirs that appeared more than tErrors times
[~,Res_int,NResint] = bwboundaries(Res_tot_fmask>tErrors,'noholes');
progress1 = NResint/10;
progress2 = NResint/40;
for nri=1:NResint
    % progression bar
    if nri==NResint
        fprintf('100%%\n')
    elseif nri>progress1
        fprintf('%.0f%%',progress1/NResint*100)
        progress1 = progress1+NResint/10;
        progress2 = NResint/40;
    elseif nri>progress1-NResint/10+progress2
        fprintf('.')
        progress2 = progress2+NResint/40;
    end
    id_res_int = find(Res_int==nri);
    % if reservoirs are known
    if nnz(realRes(id_res_int))
        Res_tot(id_res_int) = 1;
    end
end
clear nri realRes Res_int NResint progress*

[~,Lt,Nt] = bwboundaries(Res_tot,'noholes');

% if reservoirs are known
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

Res_tot(study_region==0) = 255;
% map: localisation of dams
figure, imshow(Res_tot+1,FmaskCol,'InitialMagnification',41), hold on
statsnt = regionprops(Lt,'Centroid');
for nt=1:Nt
    if IMtoDAM(nt)~=0
        % dams listed
        text(statsnt(nt).Centroid(1),statsnt(nt).Centroid(2),dams{IMtoDAM(nt),1},'Color','r')
    else
        % dams not listed
        text(statsnt(nt).Centroid(1),statsnt(nt).Centroid(2),num2str(nt),'Color','k')
        if savedata==1
            Res_tot(Lt==nt) = 0;
        end
    end
end
hold off
clear statsnt nt Lt Nt
Res_tot(study_region==0) = 0;

if savedata==1
    fprintf(' Saving data\n')
    save('statistics.mat','Res_tot_fmask','Res_tot')
    Res_tot_fmask = uint16(Res_tot_fmask);
    geotiffwrite('Res_tot_fmask',Res_tot_fmask,R,'CoordRefSysCode','EPSG:32636')
end

% ========================================================================