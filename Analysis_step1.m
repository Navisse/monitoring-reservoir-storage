%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% ANALYSIS FOR EACH DATE ACQUIRED (STEP 1)
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
% gdal_translate (http://www.gdal.org/)
% gdalwarp (http://www.gdal.org/)
% Fmask (https://github.com/prs021/fmask)
%
%=========================================================================
%%

clear variables
close all

savedata = 1; % Save data & results? (1: Yes)

%% Parameters for computation
aMin = 20; % Minimum quantity of isolated points to consider it is a reservoir


%% Source files
dir_fmask = '.\functions\FmaskSentinel\';  % Fmask directory
p = path;
dir_source = '.\Landsat_source\';
dir_landsat = '.\Landsat_results\';

dir_region = '.\region_mask\';
% mask.tif = raster of the study area (resolution from landsat images)
Iclass = geotiffinfo([dir_region,'mask.tif']);
study_region = imread([dir_region,'mask.tif']);
clear dir_region

% Landsat files
listing = dir([dir_source,'\L*']);
L = length(listing);


%% Analysis for each date acquired
% Original data
fprintf('\n-> Detecting reservoirs\n')
for l=1:L
    fprintf(['\nChecking the presence of a cloud mask in ',listing(l).name,' (%.1f %%) \n'],l/L*100)
    
    file_name = listing(l).name;
    cd([dir_landsat,'\',file_name])
    % in case aMin needs to be adjusted...
    test = dir('data.mat');
    if isempty(test)
        cd([dir_source,'\',listing(l).name])
        % capture date of the image (may need to be updated)
        year = str2double(listing(l).name(18:21));
        month = str2double(listing(l).name(22:23));
        day = str2double(listing(l).name(24:25));

        % ----------------------------------------------------------------
        % Fmask algorithm
        % Zhu, Z. and Woodcock, C. E., Improvement and Expansion of the Fmask
        % Algorithm: Cloud, Cloud Shadow, and Snow Detection for Landsats 4-7,
        % 8, and Sentinel 2 images, Remote Sensing of Environment (2014)
        % doi:10.1016/j.rse.2014.12.014
        
        % conditions (see Fmask notice: https://github.com/prs021/fmask)
        if listing(l).name(4)~='8' || year<=2013 || ...
                (year==2014 && month<12) || ...
                (year==2014 && month==12 && day<18)
            path(p, [dir_matlab,'Fmask_3_3'])
        else
            path(p, [dir_matlab,'Fmask_3_3_nothermal'])
        end
        autoFmask;
        path(p)
        
        % all files put at the same resolution
        dos(['gdal_translate ',listing(l).name,'_MTLFmask ',dir_landsat,'\',file_name,'\mask.tif']);
        Imask = geotiffinfo([dir_landsat,'\',file_name,'\mask.tif']);
        % clear land = 0
        % clear water = 1
        % cloud shadow = 2
        % snow = 3
        % cloud = 4
        % outside = 255

        d = struct('mask_i',[],'NDVI',[],'MNDWI',[]);
        
        % NDVI
        dos(['gdalwarp -overwrite -s_srs EPSG:326',num2str(Imask.Zone), ...
            ' -t_srs EPSG:',num2str(Iclass.GeoTIFFCodes.PCS),' -r near -te ',...
            num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
            num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
            ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
            ' -dstnodata 255 -of GTiff NDVI_tmp.tif NDVI.tif']);
        d.NDVI = imread('NDVI.tif');

        % MNDWI
        dos(['gdalwarp -overwrite -s_srs EPSG:326',num2str(Imask.Zone), ...
            ' -t_srs EPSG:',num2str(Iclass.GeoTIFFCodes.PCS),' -r near -te ',...
            num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
            num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
            ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
            ' -dstnodata 255 -of GTiff MNDWI_tmp.tif MNDWI.tif']);
        d.MNDWI = imread('MNDWI.tif');
        
        delete('*_MTLFmask*','NDVI*','MNDWI*')
        cd([dir_landsat,'\',file_name])
    else
        load('data.mat')
        d.mask_i = [];
        Imask = geotiffinfo('mask.tif');
    end
    clear test
    
    % --------------------------------------------------------------------
    % Erasing data outside of the YRB and possible false values

    % converting to the YRB dimensions in the original CRS: EPSG:32636 (WGS 84 / UTM zone 36N)
    dos(['gdalwarp -overwrite -s_srs EPSG:326',num2str(Imask.Zone), ...
        ' -t_srs EPSG:',num2str(Iclass.GeoTIFFCodes.PCS),' -r near -te ',...
        num2str(Iclass.BoundingBox(1,1),'%.12f'),' ',num2str(Iclass.BoundingBox(1,2),'%.12f'),' ',...
        num2str(Iclass.BoundingBox(2,1),'%.12f'),' ',num2str(Iclass.BoundingBox(2,2),'%.12f'),...
        ' -ts ',num2str(Iclass.Width,'%d'),' ',num2str(Iclass.Height,'%d'),...
        ' -dstnodata 255 -of GTiff mask.tif mask_i.tif']);
    mask_i = imread('mask_i.tif');
    delete('mask_i*')

    clear Imask

    % erasing data outside of the YRB
    mask_i(study_region==0) = 255;
    
    % erasing very small reservoirs
    [~,Lr1,Nr1] = bwboundaries(mask_i==1,'noholes');
    fprintf(' Fmask: deleting water areas < %d cells (among %d)\n0%%',aMin,Nr1)
    progress1 = Nr1/10;
    progress2 = Nr1/40;
    for n1=1:Nr1
        % progression bar
        if n1==Nr1
            fprintf('100%%\n')
        elseif n1>progress1
            fprintf('%.0f%%',progress1/Nr1*100)
            progress1 = progress1+Nr1/10;
            progress2 = Nr1/40;
        elseif n1>progress1-Nr1/10+progress2
            fprintf('.')
            progress2 = progress2+Nr1/40;
        end

        id1 = find(Lr1==n1);
        if length(id1)<=aMin
            mask_i(id1) = 0;
        end
    end
    clear n1 Nr1 Lr1 progress* id1
    
    d.mask_i = mask_i;
    clear mask_i

    if savedata==1
        fprintf(' Saving data\n')
        save('data.mat','d')
    end
    clear d
end
clear l aMin p_prctile

%=========================================================================