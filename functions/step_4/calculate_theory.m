% ========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FUNCTION FOR STEP 4: calculating theory statistics on each reservoir
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

function r_nt = calculate_theory(r_nt,L,dams,nt,Nt,Lt,IMtoDAM_nt,DEM,BasicArea)

fprintf('Reservoir %d/%d: [%s]\n',nt,Nt,dams{IMtoDAM_nt,1})

% coordinates of each dam's pixels
[rt,~] = find(Lt==nt);
lt = length(rt);

% table of "flooding implications" statistics
r_nt.flooded = zeros(lt,3+L);
% expected results (from the DEM)
r_nt.theory = struct('H',zeros(100,1),'A',zeros(100,1),'V',zeros(100,1));
dem_r = DEM(Lt==nt);
level_pixels = sort(dem_r);%[m]
r_nt.theory.H(1) = level_pixels(1);% [m]
r_nt.theory.A(1) = 0;% [m2]
r_nt.theory.V(1) = 0;% [hm3]
ind_z = 1;
for z=1:length(dem_r)% for each water pixel
    if level_pixels(z)~=r_nt.theory.H(ind_z)
        % ->  for each water level
        ind_z = ind_z+1;
        r_nt.theory.H(ind_z) = level_pixels(z);% [m]
        r_nt.theory.A(ind_z) = (nnz(dem_r<r_nt.theory.H(ind_z))+nnz(dem_r==r_nt.theory.H(ind_z))/2)*BasicArea;% [m2]
        r_nt.theory.V(ind_z) = sum((dem_r<=r_nt.theory.H(ind_z)).*(r_nt.theory.H(ind_z)-dem_r))*BasicArea/1e6;% [hm3]
    end
end
r_nt.theory.H = r_nt.theory.H(1:ind_z);
r_nt.theory.A = r_nt.theory.A(1:ind_z);
r_nt.theory.V = r_nt.theory.V(1:ind_z);
clear z dem_r level_pixels ind_z

end

% ========================================================================