% ========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FUNCTION FOR STEP 5: assembling through no-data
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

function img2_l = assemble(img2_l,savedata,time,listing_l_name,dir_landsat, ...
    study_region,Iclass,Res_tot,h_res,BasicArea,qAs)

load([dir_landsat,listing_l_name,'\data.mat'],'d')

Res_i = double(d.Res_i);
mask_i = double(d.mask_i);
mask_i = mask_i.*Res_tot;

% filling step
NoData = zeros(Iclass.Height,Iclass.Width);
% NoData = no water nor land
NoData((Res_i~=1).*(mask_i~=0)==1) = 1;

Res = Res_i;
[~,Lr,~] = bwboundaries(Res,'noholes');
% nodata on reservoirs
[~,Lrnd,Nrnd] = bwboundaries(Res_tot.*NoData,'noholes');
error_res = 0*Res_tot;
clear NoData *_i

% 3 conditions for assembling:
% - no-data pixels,
% - in Res_tot,
% - h_res < level detected.
for n=1:Nrnd
    [rn,cn] = find(Lrnd==n);
    K = length(rn);
    k = 0; % progression in nodata areas
    A = [];
    level = [];
    while k<K
        k = k+1;
        error_res(rn(k),cn(k)) = 1; % => nodata & dam at this location
        % detecting reservoirs just next to nodata areas
        a = max(max(Lr(max(rn(k)-1,1):min(rn(k)+1,Iclass.Height), ...
            max(cn(k)-1,1):min(cn(k)+1,Iclass.Width))));
        if a~=0 && ~any(A==a)% dam different than the previous one detected
            level_tmp = h_res(Lr==a);
            % taking qAs quantile for level value
            level_tmp = quantile(level_tmp,qAs);
            % level for each dam's part
            level = [level level_tmp];
            clear level_tmp
            A = [A, a];
        end
        clear a
    end
    % if dam detected
    if ~isempty(level)
        Res((Lrnd==n).*(h_res<max(level))==1) = 1;
    elseif ~isempty(level)
        for idxA = A
            Res(Lr==idxA) = 0;
        end
        clear idxA
    end
    clear A count_A k K level rn cn
end
clear Lr Lrnd Nrnd n
Res = imfill(Res,'holes');
Res = Res.*Res_tot;

% --------------------------------------------------------------------
% Statistics

[~,Res,Nres] = bwboundaries(Res,'noholes');
img2_l.stats = regionprops(Res,'Centroid','Area');
img2_l.pct = struct('Lmask_i',cell(1,Nres),'Wmask_i',cell(1,Nres), ...
    'CSmask_i',cell(1,Nres),'Smask_i',cell(1,Nres), ...
    'Cmask_i',cell(1,Nres),'Omask_i',cell(1,Nres), ...
    'Res_i',cell(1,Nres),'Res_f',cell(1,Nres));

% fprintf(' Calculating image statistics\n')
for n=1:Nres
    Ln = Res;
    Ln(Res~=n) = 0;

    % reservoir has/could have been reconstructed?
    % *4: LT4, *5: LT5, *7: LE7, *8: LC8
    if nnz(Ln.*error_res)
        % sthg was hiding reservoir n
        img2_l.stats(n).trust = 0 + str2double(listing_l_name(3));
    else
        % reservoir n was free of clouds or anythg else
        img2_l.stats(n).trust = 10 + str2double(listing_l_name(3));
    end
    clear Ln

    % vector of the hight of each reservoir's pixels
    hresn = h_res(Res==n);
    hresn = hresn(hresn~=0);
    % bottom's coordinates
    if ~isempty(hresn)
        [Rb,Cb] = find(h_res.*(Res==n)==min(hresn));
    else
        [Rb,Cb] = find(Res==n);
    end
    d2 = (Cb-img2_l.stats(n).Centroid(1)).^2+(Rb-img2_l.stats(n).Centroid(2)).^2;
    id_d2 = find(d2==max(d2));
    img2_l.stats(n).iBottom = max(Rb(id_d2));
    img2_l.stats(n).jBottom = min(Cb(id_d2));
    clear id_d2 d2 Rb Cb hresn

    img2_l.stats(n).xBottom = Iclass.BoundingBox(1,1)+(img2_l.stats(n).jBottom-0.5)*Iclass.PixelScale(1);% [°]
    img2_l.stats(n).yBottom = Iclass.BoundingBox(1,2)+(Iclass.Height-img2_l.stats(n).iBottom+0.5)*Iclass.PixelScale(2);% [°]
    img2_l.stats(n).zBottom = h_res(img2_l.stats(n).iBottom,img2_l.stats(n).jBottom);% [m]

    % pixels that are immersed after assembling
    img2_l.stats(n).Area = img2_l.stats(n).Area*BasicArea;% [m^2]
    % position
    img2_l.stats(n).Immersed = find(Res==n);% [idx]
    
    % for different detections:
    img2_l.pct(n).Lmask_i = nnz((Res==n).*(d.mask_i==0));
    img2_l.pct(n).Wmask_i = nnz((Res==n).*(d.mask_i==1));
    img2_l.pct(n).CSmask_i = nnz((Res==n).*(d.mask_i==2));
    img2_l.pct(n).Smask_i = nnz((Res==n).*(d.mask_i==3));
    img2_l.pct(n).Cmask_i = nnz((Res==n).*(d.mask_i==4));
    img2_l.pct(n).Omask_i = nnz((Res==n).*(d.mask_i==255));
    
    img2_l.pct(n).Res_i = nnz((Res==n).*(d.Res_i==1));
    img2_l.pct(n).Res_f = nnz(Res==n);
end
clear n Nres error_res

% acquiring date
img2_l.yd = [str2double(listing_l_name(10:13)) ,...
    str2double(listing_l_name(14:16))];

Res(Res~=0) = 1;
Res(study_region==0) = 255;
d.(['Res_f',time]) = uint8(Res);
clear Res

if savedata==1
    % data saving:
    % d ~ data for each acquisition date
    save([dir_landsat,listing_l_name,'\data.mat'],'d')
end
clear d

end

% ========================================================================