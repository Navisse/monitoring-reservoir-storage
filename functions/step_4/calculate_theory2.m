% ========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FUNCTION FOR STEP 4: calculating statistics on each reservoir
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

function [r2_nt_theory, h_res_nt, stats_res_nt, i_end_nt, idx_p] = calculate_theory2(r_nt, ...
    r2ASTER_nt,r2SRTMC_nt,r2SRTMX_nt,id_srtmc,id_srtmx,idx_p, ...
    savedata,graph,scrsz,L,n_fit,x_fit,y0_fit,span_fit, ...
    dams,nt,Nt,Lt,IMtoDAM_nt,DAMtoIM,Iclass,BasicArea)

fprintf('Reservoir %d/%d: [%s]\n',nt,Nt,dams{IMtoDAM_nt,1})

% immersion statistics
stats_res_nt = zeros(Iclass.Height,Iclass.Width);
% corrected DEM
h_res_nt = zeros(Iclass.Height,Iclass.Width);

% coordinates of each dam's pixels
[rt,ct] = find(Lt==nt);
lt = length(rt);

% H2, A2, V2, S2: observed and basic mean
% H2f, V2f: obtained after fitting
r2_nt_theory = struct('H2',zeros(100,1),'A2',zeros(100,1), ...
    'V2',zeros(100,1),'S2',zeros(100,1), ...
    'H2f',[],'V2f',zeros(100,1));

% DEM for method 2
% table only for images where the dam nt has been detected
im_table = r_nt.flooded(:,[0 0 0 ones(1,L)].*(sum(r_nt.flooded==1)~=0)==1);

stats_h = zeros(lt,1);
for h=1:lt% for each pixel of the reservoir: % of the time it is visible
    stats_h(h) = nnz(im_table(h,:)==1)/nnz(im_table(h,:)~=-1);
end
clear h im_table

[~,I_sh] = sort(stats_h,'descend');% sorting from most often immersed to least often immersed
level_pixels2 = r_nt.flooded(:,3);% hight (or above bottom) [m]
i_sh1 = I_sh(1);
h_tmp = level_pixels2(i_sh1);
lh_tmp = 0;
f_tmp = stats_h(i_sh1);
if savedata==1
    stats_res_nt(rt(i_sh1),ct(i_sh1)) = stats_h(i_sh1);
end
% r2_nt_theory.H2(1) = 0;% [m]
r2_nt_theory.H2(1) = h_tmp;% [m]
r2_nt_theory.A2(1) = 0;% [m2]
r2_nt_theory.V2(1) = 0;% [hm3]
r2_nt_theory.S2(1) = 1;% stats on immersion [0,1]
ind_m = 1;
clear i_sh1
for i_sh = I_sh(2:end)'% for each pixel from potentially deepest to highest
    % conditions to average dem values
    if max(f_tmp)==stats_h(i_sh)
        h_tmp = [h_tmp, level_pixels2(i_sh)];
        f_tmp = [f_tmp, stats_h(i_sh)];% in case larger interval...
    else
        ind_m = ind_m+1;
        r2_nt_theory.H2(ind_m) = mean(h_tmp);% [m]
        r2_nt_theory.A2(ind_m) = r2_nt_theory.A2(ind_m-1) + (length(h_tmp)+lh_tmp)*BasicArea/2;% [m2]
        r2_nt_theory.V2(ind_m) = sum((r2_nt_theory.H2(ind_m)-r2_nt_theory.H2(1:ind_m)).* ...
            (r2_nt_theory.A2(1:ind_m)-[0;r2_nt_theory.A2(1:ind_m-1)]))/1e6;% [hm3]
        r2_nt_theory.S2(ind_m) = max(f_tmp);
        lh_tmp = length(h_tmp);
        h_tmp = level_pixels2(i_sh);
        f_tmp = stats_h(i_sh);
    end
    if savedata==1
        stats_res_nt(rt(i_sh),ct(i_sh)) = stats_h(i_sh);
    end
end
ind_m = ind_m+1;
r2_nt_theory.H2(ind_m) = mean(h_tmp);% [m]
r2_nt_theory.A2(ind_m) = r2_nt_theory.A2(ind_m-1) + (length(h_tmp)+lh_tmp)*BasicArea/2;% [m2]
r2_nt_theory.V2(ind_m) = sum((r2_nt_theory.H2(ind_m)-r2_nt_theory.H2(1:ind_m)).* ...
    (r2_nt_theory.A2(1:ind_m)-[0;r2_nt_theory.A2(1:ind_m-1)]))/1e6;% [hm3]
r2_nt_theory.S2(ind_m) = max(f_tmp);
r2_nt_theory.H2 = r2_nt_theory.H2(1:ind_m);
r2_nt_theory.A2 = r2_nt_theory.A2(1:ind_m);
r2_nt_theory.V2 = r2_nt_theory.V2(1:ind_m);
r2_nt_theory.S2 = r2_nt_theory.S2(1:ind_m);
clear I_sh i_sh h_tmp lh_tmp f_tmp ind_m stats_h level_pixels2

% Fitting
if any(nt==DAMtoIM(id_srtmx))
    name_dem = 'SRTM-X';
elseif any(nt==DAMtoIM(id_srtmc))
    name_dem = 'SRTM-C';
else
    name_dem  = 'ASTER';
end

if any(nt==DAMtoIM([1 28]))
    % Local polynomial fitting
    x = r2_nt_theory.A2;
    y = r2_nt_theory.H2;
    
    % *b = part considered for the basic fitting
    xb = []; yb = [];
    idb = find(x/max(x)>=x_fit(IMtoDAM_nt));
    lx10 = length(x)/10; % decomposition depends on the reservoir
    step_x = ceil(idb(1)/lx10);
    up_lim = prctile(y(round((step_x-1)*lx10)+1:round((step_x)*lx10)),80);
    lo_lim = prctile(y(round((step_x-1)*lx10)+1:round((step_x)*lx10)),20);
    for i_idb = idb'
        if i_idb>round(step_x*lx10)
            up_lim = prctile(y(round(step_x*lx10)+1: ...
                round((step_x+1)*lx10)),90);
            lo_lim = prctile(y(round(step_x*lx10)+1: ...
                round((step_x+1)*lx10)),10);
            step_x = step_x+1;
        end
        if y(i_idb)<up_lim && y(i_idb)>lo_lim
            xb = [xb; x(i_idb)];
            yb = [yb; y(i_idb)];
        end
    end
    clear idb i_idb step_x up_lim lo_lim lx10

    % rloess fit
    idh = find(x/max(x)>=x_fit(IMtoDAM_nt),1);
    H2f = smooth(xb,yb,span_fit(IMtoDAM_nt),'rloess');
    if x_fit(IMtoDAM_nt)~=0
        r2_nt_theory.H2f = [interp1([0;x(idh)], [y0_fit(IMtoDAM_nt);min(H2f)],x(1:idh-1));H2f];
    else
        r2_nt_theory.H2f = H2f;
    end
    r2_nt_theory.H2f = interp1([x(1:idh-1);xb],r2_nt_theory.H2f,x);
    r2_nt_theory.H2f(isnan(r2_nt_theory.H2f).*(x>median(x))==1) = max(r2_nt_theory.H2f);
    r2_nt_theory.H2f(isnan(r2_nt_theory.H2f).*(x<median(x))==1) = min(r2_nt_theory.H2f);
    L2f = length(r2_nt_theory.H2f);
    for i2f = 2:L2f
        if r2_nt_theory.H2f(i2f)<r2_nt_theory.H2f(i2f-1)
            r2_nt_theory.H2f(i2f) = r2_nt_theory.H2f(i2f-1);
        end
    end
    clear L2f i2f

    ybr = H2f;
    R2 = sum((yb-mean(yb)).*(ybr-mean(ybr)))^2/(sum((yb-mean(yb)).^2).*sum((ybr-mean(ybr)).^2));
    fprintf(['R^2 = ',num2str(R2,'%.2f'),'\n'])
    if graph==1
        idx_p = idx_p+1;
        % Plot H
        figure('Position',[scrsz(3)/5 scrsz(4)/4 scrsz(3)*2/5 scrsz(4)/3]), hold on
        
%         % plot uncorrected curve
%         plot(r_nt.theory.A*1e-6,r_nt.theory.H,'--k','DisplayName','H-A w/o sorting')
        % Plot original data
%         if x_fit(IMtoDAM_nt)~=0
%             plot(x(1:idh-1)*1e-6,y(1:idh-1),'.k','DisplayName',name_dem)
%         end
        
        % Original sorted data
        if strcmp(name_dem,'ASTER')
            plot(r2ASTER_nt.theory.A2*1e-6,r2ASTER_nt.theory.H2,'ok','markerfacecolor','b','markersize',5,'DisplayName','ASTER')
        elseif strcmp(name_dem,'SRTM-C')
            plot(r2SRTMC_nt.theory.A2*1e-6,r2SRTMC_nt.theory.H2,'sk','markerfacecolor','g','markersize',5,'DisplayName','SRTM-C')
        else
            plot(r2SRTMX_nt.theory.A2*1e-6,r2SRTMX_nt.theory.H2,'>k','markerfacecolor','y','markersize',5,'DisplayName','SRTM-X')
        end
        
        plot(xb*1e-6,yb,'or','markerfacecolor','r','markersize',2,'DisplayName',[name_dem,' for regression'])
        
        % Plot fitted data
        plot(x*1e-6,r2_nt_theory.H2f,'-k','linewidth',2,'DisplayName',['H_c: LPR (span ',num2str(span_fit(IMtoDAM_nt)),')'])
        lgd = legend('show','Location','SouthEast');
        lgd.FontSize = 11;
        grid on
        xlabel('A [km^2]','FontSize',12)
        ylabel('H [m]','FontSize',12)
        set(gca,'FontSize',12)
        hold off
    end
    clear x y H2f idh H2min lgd
else
    % Polynomial fitting
    x = r2_nt_theory.A2;
    y = r2_nt_theory.H2;
    % *a = part ignored for the fitting
    ida = find(x/max(x)<x_fit(IMtoDAM_nt));
    xa = x(ida); ya = y(ida);
    % *b = part considered for the basic fitting
    xb = []; yb = [];
    idb = find(x/max(x)>=x_fit(IMtoDAM_nt));
    lx10 = length(x)/10;  % decomposition depends on the reservoir
    step_x = ceil(idb(1)/lx10);
    up_lim = prctile(y(round((step_x-1)*lx10)+1:round((step_x)*lx10)),80);
    lo_lim = prctile(y(round((step_x-1)*lx10)+1:round((step_x)*lx10)),20);
    for i_idb = idb'
        if i_idb>round(step_x*lx10)
            up_lim = prctile(y(round(step_x*lx10)+1: ...
                round((step_x+1)*lx10)),90);
            lo_lim = prctile(y(round(step_x*lx10)+1: ...
                round((step_x+1)*lx10)),10);
            step_x = step_x+1;
        end
        if y(i_idb)>up_lim || y(i_idb)<lo_lim
            xa = [xa; x(i_idb)];
            ya = [ya; y(i_idb)];
        else
            xb = [xb; x(i_idb)];
            yb = [yb; y(i_idb)];
        end
    end
    clear ida idb i_idb step_x up_lim lo_lim lx10

    xb = xb(:); %reshape the data into a column vector
    yb = yb(:);
    % 'C' is the Vandermonde matrix for 'x'
    V(:,n_fit(IMtoDAM_nt)+1) = ones(length(xb),1,class(xb));
    for j = n_fit(IMtoDAM_nt):-1:1
        V(:,j) = xb.*V(:,j+1);
    end
    clear j
    C = V;
    % 'd_tv' is the vector of target values, 'y'.
    d_tv = yb;

    % There are no inequality constraints in this case, i.e., 
    A = [];
    b = [];

    options = optimset('Algorithm','active-set','LargeScale','off');
    if isnan(y0_fit(IMtoDAM_nt))
        p = lsqlin( C, d_tv, A, b, [],[],[],[],[],options);
    else
        % We use linear equality constraints to force the curve to hit the required point. In
        % this case, 'Aeq' is the Vandermoonde matrix for 'x0'
        x0 = 0;
        Aeq = x0.^(n_fit(IMtoDAM_nt):-1:0);
        y0 = y0_fit(IMtoDAM_nt);
        beq = y0;
        p = lsqlin( C, d_tv, A, b, Aeq, beq,[],[],[],options);
    end

    % We can then use POLYVAL to evaluate the fitted curve
    r2_nt_theory.H2f = polyval( p, x );

    ybr = polyval(p,xb);
    R2 = sum((yb-mean(yb)).*(ybr-mean(ybr)))^2/(sum((yb-mean(yb)).^2).*sum((ybr-mean(ybr)).^2));
    fprintf(['R^2 = ',num2str(R2,'%.2f'),'\n'])
    if graph==1
        idx_p = idx_p+1;
        % Plot H
        figure('Position',[scrsz(3)/5 scrsz(4)/4 scrsz(3)*2/5 scrsz(4)/3]), hold on
        
%         % plot uncorrected curve
%         plot(r_nt.theory.A*1e-6,r_nt.theory.H,'--k','DisplayName','H-A w/o sorting')
        % Plot original data
%         if x_fit(IMtoDAM_nt)~=0
%         if ~isempty(xa)
%             plot(xa*1e-6,ya,'.k','DisplayName',[name_dem,' ignored in fit'])
%         end

        % Original sorted data
        if strcmp(name_dem,'ASTER')
            plot(r2ASTER_nt.theory.A2*1e-6,r2ASTER_nt.theory.H2,'ok','markerfacecolor','b','markersize',5,'DisplayName','ASTER')
        elseif strcmp(name_dem,'SRTM-C')
            plot(r2SRTMC_nt.theory.A2*1e-6,r2SRTMC_nt.theory.H2,'sk','markerfacecolor','g','markersize',5,'DisplayName','SRTM-C')
        else
            plot(r2SRTMX_nt.theory.A2*1e-6,r2SRTMX_nt.theory.H2,'>k','markerfacecolor','y','markersize',5,'DisplayName','SRTM-X')
        end
        
        plot(xb*1e-6,yb,'or','markerfacecolor','r','markersize',2,'DisplayName',[name_dem,' for regression'])
        
%         % Plot points to go through
%         if ~isnan(y0_fit(IMtoDAM_nt))
%             plot(x0*1e-6,y0,'dr','MarkerFaceColor','k','DisplayName','Constraints') 
%         end
        
        % Plot fitted data
        plot(x*1e-6,r2_nt_theory.H2f,'-k','linewidth',2,'DisplayName',['H_c: PR (deg. ',num2str(n_fit(IMtoDAM_nt)),')'])
        lgd = legend('show','Location','SouthEast');
        lgd.FontSize = 11;
        grid on
        xlabel('A [km^2]','FontSize',12)
        ylabel('H [m]','FontSize',12)
        set(gca,'FontSize',12)
        hold off
    end
    clear V x y xa ya xb yb x0* y0* C d_tv A b Aeq* beq* p options ybr coef_r sse sst R2 R2plot lgd
end
clear name_dem
if savedata==1
    h_res_nt((Lt==nt).*(stats_res_nt==r2_nt_theory.S2(1))==1) = r2_nt_theory.H2f(1);
end
r2_nt_theory.V2f(1) = 0;
for i=2:length(r2_nt_theory.H2f)
    if savedata==1
        h_res_nt((Lt==nt).*(stats_res_nt==r2_nt_theory.S2(i))==1) = r2_nt_theory.H2f(i);
    end
    r2_nt_theory.V2f(i) = sum((max(r2_nt_theory.H2f(1:i))-r2_nt_theory.H2f(1:i)).* ...
        (r2_nt_theory.A2(1:i)-[0;r2_nt_theory.A2(1:i-1)]))/1e6;% [hm3]
end
r2_nt_theory.V2f = r2_nt_theory.V2f(1:i);

i_end_nt = i;
clear i

end

% ========================================================================