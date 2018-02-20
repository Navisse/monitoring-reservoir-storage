% ========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FUNCTION FOR STEP 4: calculating "flooded" statistics on each reservoir
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

function r_nt_flooded = calculate_flooded(r_nt_flooded,time,l,Lt,nt,DEM, ...
    d_mask_i,d_Res_i,d_Res_f)

[rt,ct] = find(Lt==nt);
lt = length(rt);

for k=1:lt% for each reservoir's pixel
    if l==1
        r_nt_flooded(k,1) = ct(k);% x
        r_nt_flooded(k,2) = rt(k);% ymax-y
        r_nt_flooded(k,3) = DEM(rt(k),ct(k));% z
    end
    if (isempty(time) && d_Res_i(rt(k),ct(k))==1) || ...
            (~isempty(time) && d_Res_f(rt(k),ct(k))==1)
        % if there is water
        r_nt_flooded(k,l+3) = 1;
    elseif d_mask_i(rt(k),ct(k))~=0
        % if cloud, cloud shadow, snow or no-data
        r_nt_flooded(k,l+3) = -1;
    end
end

end

% ========================================================================