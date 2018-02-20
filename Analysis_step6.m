%=========================================================================
% MONITORING RESERVOIRS STORAGE IN A BASIN
% FINAL ANALYSIS FOR EACH RESERVOIR (STEP 6)
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

clear variables
close all

savedata = 1; % Save data & results?
time = '';


%% Parameters for computation
pSa = 0.99; % Minimum percentage of main area indivisble to be saved


%% Source files
dir_region = '.\region_mask\';
    
fprintf('-> Loading data\n')
load('statistics','Res_tot')
load(['statistics_2',time],'h_res','stats_res','img2','r2')
% dams
[~,~,dams] = xlsread([dir_region,'dams_info.xlsx'],'A2:C3');

Iclass = geotiffinfo([dir_region,'mask.tif']);

[~,Lt,Nt] = bwboundaries(Res_tot,'noholes');

[ndd,~] = size(dams);% nb of dams according to studies
% dam number in the xls file for each reservoir in the image
IMtoDAM = zeros(1,Nt);
% reservoir number on the image for each dam
DAMtoIM = zeros(1,ndd);
for nd=1:ndd
    if ~isnan(dams{nd,14})% if the location of the dam is known
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
dir_landsat = '.\Landsat_results\';
listing = dir([dir_landsat,'\L*']);
L = length(listing);


%% Analysis for each reservoir
% data to consider in reservoirs' statistics
r3 = struct('stats',cell(1,Nt),'DS',cell(1,Nt));
for nt=1:Nt
    if ~isempty(r2(nt).theory)
        r2(nt).stats = struct('i',NaN(1,L),'j',NaN(1,L),'x',NaN(1,L),'y',NaN(1,L), ...
            'z',NaN(1,L),'doy',NaN(1,L),'yr',NaN(1,L),'H',NaN(1,L), 'H_min',[],'H_max',[], ...
            'A',NaN(1,L),'V',NaN(1,L),'V_min',[],'V_max',[],'trust',NaN(1,L),'T',[]);
        r2(nt).pct = struct('Lmask_i',NaN(1,L),'Wmask_i',NaN(1,L),'CSmask_i',NaN(1,L), ...
            'Smask_i',NaN(1,L),'Cmask_i',NaN(1,L),'Omask_i',NaN(1,L), ...
            'Res_i',NaN(1,L),'Res_f',NaN(1,L));
        r3(nt).stats = struct('doy',[],'yr',[],'T',[],'H',[],'A',[],'V',[],'trust',[]);
    end
end
clear nt

DateString = cell(L,1);
fprintf('\n-> Calculating reservoirs statistics\n0%%')
progress1 = L/10;
progress2 = L/40;
for l=1:L% for each acquisition date
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
    
    l_img = length(img2(l).stats);
    list_nt = zeros(1,l_img);
    i = zeros(1,l_img);
    j = zeros(1,l_img);
    x = zeros(1,l_img);
    y = zeros(1,l_img);
    z = zeros(1,l_img);
    doy = zeros(1,l_img);
    yr = zeros(1,l_img);
    A = zeros(1,l_img);
    data_trust = zeros(1,l_img);
    pct_Lmask_i = zeros(1,l_img);
    pct_Wmask_i = zeros(1,l_img);
    pct_CSmask_i = zeros(1,l_img);
    pct_Smask_i = zeros(1,l_img);
    pct_Cmask_i = zeros(1,l_img);
    pct_Omask_i = zeros(1,l_img);
    pct_Res_i = zeros(1,l_img);
    pct_Res_f = zeros(1,l_img);
    
    for n=1:l_img% for each reservoir in current image
        list_nt(n) = Lt(img2(l).stats(n).iBottom,img2(l).stats(n).jBottom);
        
        i(n) = img2(l).stats(n).iBottom;
        j(n) = img2(l).stats(n).jBottom;
        x(n) = img2(l).stats(n).xBottom;
        y(n) = img2(l).stats(n).yBottom;
        z(n) = img2(l).stats(n).zBottom;
        doy(n) = img2(l).yd(2);
        yr(n) = img2(l).yd(1);
        A(n) = img2(l).stats(n).Area;
        data_trust(n) = img2(l).stats(n).trust;
        
        pct_Lmask_i(n) = img2(l).pct(n).Lmask_i;
        pct_Wmask_i(n) = img2(l).pct(n).Wmask_i;
        pct_CSmask_i(n) = img2(l).pct(n).CSmask_i;
        pct_Smask_i(n) = img2(l).pct(n).Smask_i;
        pct_Cmask_i(n) = img2(l).pct(n).Cmask_i;
        pct_Omask_i(n) = img2(l).pct(n).Omask_i;
        pct_Res_i(n) = img2(l).pct(n).Res_i;
        pct_Res_f(n) = img2(l).pct(n).Res_f;
    end
    clear n l_img
    while nnz(list_nt)
        nt = nonzeros(list_nt);
        nt = nt(1);
        id_1 = find(list_nt==nt);
        list_nt(id_1(A(id_1)/sum(A(id_1))<pSa)) = 0;
        id_1 = find(list_nt==nt);
        if length(i(id_1))==1 && ~isempty(r2(nt).theory)
            r2(nt).stats.i(l) = i(id_1);
            r2(nt).stats.j(l) = j(id_1);
            r2(nt).stats.x(l) = x(id_1);
            r2(nt).stats.y(l) = y(id_1);
            r2(nt).stats.z(l) = z(id_1);
            r2(nt).stats.doy(l) = doy(id_1);
            r2(nt).stats.yr(l) = yr(id_1);
            r2(nt).stats.A(l) = A(id_1);
            r2(nt).stats.trust(l) = data_trust(id_1);
            
            r2(nt).pct.Lmask_i(l) = pct_Lmask_i(id_1);
            r2(nt).pct.Wmask_i(l) = pct_Wmask_i(id_1);
            r2(nt).pct.CSmask_i(l) = pct_CSmask_i(id_1);
            r2(nt).pct.Smask_i(l) = pct_Smask_i(id_1);
            r2(nt).pct.Cmask_i(l) = pct_Cmask_i(id_1);
            r2(nt).pct.Omask_i(l) = pct_Omask_i(id_1);
            r2(nt).pct.Res_i(l) = pct_Res_i(id_1);
            r2(nt).pct.Res_f(l) = pct_Res_f(id_1);
            % index of closest values to the area detected
            idA_lower = find(r2(nt).theory.A2<=A(id_1),1,'last');
            idA_higher = find(r2(nt).theory.A2>A(id_1),1);
            if isempty(idA_higher)
                r2(nt).stats.H(l) = r2(nt).theory.H2f(idA_lower);
                r2(nt).stats.V(l) = r2(nt).theory.V2f(idA_lower);
            else
                r2(nt).stats.H(l) = (A(id_1)-r2(nt).theory.A2(idA_lower))* ...% (A-Amin)*
                    (r2(nt).theory.H2f(idA_higher)-r2(nt).theory.H2f(idA_lower))/ ...% (Hmax-Hmin)/
                    (r2(nt).theory.A2(idA_higher)-r2(nt).theory.A2(idA_lower))+ ...% (Amax-Amin)+
                    r2(nt).theory.H2f(idA_lower);% Hmin
                r2(nt).stats.V(l) = (A(id_1)-r2(nt).theory.A2(idA_lower))* ...% (A-Amin)*
                    (r2(nt).theory.V2f(idA_higher)-r2(nt).theory.V2f(idA_lower))/ ...% (Vmax-Vmin)/
                    (r2(nt).theory.A2(idA_higher)-r2(nt).theory.A2(idA_lower))+ ...% (Amax-Amin)+
                    r2(nt).theory.V2f(idA_lower);% Vmin
            end
            clear idA*
        end
        list_nt(id_1) = 0;
        clear id_1 nt
    end
    clear list_nt i j x y z doy yr A data_trust
    E = eomday(img2(l).yd(1),1:12);
    E_sum = zeros(1,12);
    for k=1:12
        E_sum(k) = sum(E(1:k));
    end
    clear k
    year_yd = img2(l).yd(1);
    month_yd = find(E_sum>=img2(l).yd(2),1);
    day_yd = img2(l).yd(2)-E_sum(month_yd)+E(month_yd);
    if month_yd<10
        m_yd = ['0',num2str(month_yd)];
    else
        m_yd = num2str(month_yd);
    end
    if day_yd<10
        d_yd = ['0',num2str(day_yd)];
    else
        d_yd = num2str(day_yd);
    end
    DateString{l} = [num2str(year_yd),'-',m_yd,'-',d_yd];
    clear E E_sum *_yd
end
clear l BasicArea pSa progress*

formatIn = 'yyyy-mm-dd';
T = datenum(DateString,formatIn);
clear formatIn
[T,I] = sort(T);
DateString = DateString(I);

fprintf('\n-> Saving reservoirs statistics\n0%%')
progress1 = Nt/10;
progress2 = Nt/40;
for nt=1:Nt
    % progression bar
    if nt==Nt
        fprintf('100%%\n')
    elseif nt>progress1
        fprintf('%.0f%%',progress1/Nt*100)
        progress1 = progress1+Nt/10;
        progress2 = Nt/40;
    elseif nt>progress1-Nt/10+progress2
        fprintf('.')
        progress2 = progress2+Nt/40;
    end
    
    if ~isempty(r2(nt).theory)
        r2(nt).stats.i = r2(nt).stats.i(I);
        r2(nt).stats.j = r2(nt).stats.j(I);
        r2(nt).stats.x = r2(nt).stats.x(I);
        r2(nt).stats.y = r2(nt).stats.y(I);
        r2(nt).stats.z = r2(nt).stats.z(I);
        r2(nt).stats.doy = r2(nt).stats.doy(I);
        r2(nt).stats.yr = r2(nt).stats.yr(I);
        r2(nt).stats.A = r2(nt).stats.A(I);
        r2(nt).stats.trust = r2(nt).stats.trust(I);
        r2(nt).stats.H = r2(nt).stats.H(I);
        r2(nt).stats.V = r2(nt).stats.V(I);
        r2(nt).pct.Lmask_i = r2(nt).pct.Lmask_i(I);
        r2(nt).pct.Wmask_i = r2(nt).pct.Wmask_i(I);
        r2(nt).pct.CSmask_i = r2(nt).pct.CSmask_i(I);
        r2(nt).pct.Smask_i = r2(nt).pct.Smask_i(I);
        r2(nt).pct.Cmask_i = r2(nt).pct.Cmask_i(I);
        r2(nt).pct.Omask_i = r2(nt).pct.Omask_i(I);
        r2(nt).pct.Res_i = r2(nt).pct.Res_i(I);
        r2(nt).pct.Res_f = r2(nt).pct.Res_f(I);

        idw = find(~isnan(r2(nt).stats.H));
        r2(nt).stats.i = r2(nt).stats.i(idw);
        r2(nt).stats.j = r2(nt).stats.j(idw);
        r2(nt).stats.x = r2(nt).stats.x(idw);
        r2(nt).stats.y = r2(nt).stats.y(idw);
        r2(nt).stats.z = r2(nt).stats.z(idw);
        r2(nt).stats.doy = r2(nt).stats.doy(idw);
        r2(nt).stats.yr = r2(nt).stats.yr(idw);
        r2(nt).stats.A = r2(nt).stats.A(idw);
        r2(nt).stats.trust = r2(nt).stats.trust(idw);
        r2(nt).stats.H = r2(nt).stats.H(idw);
        r2(nt).stats.V = r2(nt).stats.V(idw);
        r2(nt).pct.Lmask_i = r2(nt).pct.Lmask_i(idw);
        r2(nt).pct.Wmask_i = r2(nt).pct.Wmask_i(idw);
        r2(nt).pct.CSmask_i = r2(nt).pct.CSmask_i(idw);
        r2(nt).pct.Smask_i = r2(nt).pct.Smask_i(idw);
        r2(nt).pct.Cmask_i = r2(nt).pct.Cmask_i(idw);
        r2(nt).pct.Omask_i = r2(nt).pct.Omask_i(idw);
        r2(nt).pct.Res_i = r2(nt).pct.Res_i(idw);
        r2(nt).pct.Res_f = r2(nt).pct.Res_f(idw);
        r2(nt).stats.T = T(idw)';
        r2(nt).DS = DateString(idw);
        clear idw
        
        l_T = length(r2(nt).stats.T);
        t_tmp = 1;
        id_f = 1;
        while t_tmp<=l_T % When several images at the same date
            id_t = find(r2(nt).stats.T==r2(nt).stats.T(t_tmp));
            l_id_t = length(id_t);
            r2_stats_i(id_f) = r2(nt).stats.i(id_t(1));
            r2_stats_j(id_f) = r2(nt).stats.j(id_t(1));
            r2_stats_x(id_f) = r2(nt).stats.i(id_t(1));
            r2_stats_y(id_f) = r2(nt).stats.y(id_t(1));
            r2_stats_z(id_f) = r2(nt).stats.z(id_t(1));
            r2_stats_doy(id_f) = r2(nt).stats.doy(id_t(1));
            r2_stats_yr(id_f) = r2(nt).stats.yr(id_t(1));
            r2_stats_A(id_f) = mean(r2(nt).stats.A(id_t));
            r2_stats_trust(id_f) = min(r2(nt).stats.trust(id_t));
            r2_stats_H(id_f) = mean(r2(nt).stats.H(id_t));
            r2_stats_V(id_f) = mean(r2(nt).stats.V(id_t));
            r2_pct_Lmask_i(id_f) = mean(r2(nt).pct.Lmask_i(id_t));
            r2_pct_Wmask_i(id_f) = mean(r2(nt).pct.Wmask_i(id_t));
            r2_pct_CSmask_i(id_f) = mean(r2(nt).pct.CSmask_i(id_t));
            r2_pct_Smask_i(id_f) = mean(r2(nt).pct.Smask_i(id_t));
            r2_pct_Cmask_i(id_f) = mean(r2(nt).pct.Cmask_i(id_t));
            r2_pct_Omask_i(id_f) = mean(r2(nt).pct.Omask_i(id_t));
            r2_pct_Res_i(id_f) = mean(r2(nt).pct.Res_i(id_t));
            r2_pct_Res_f(id_f) = mean(r2(nt).pct.Res_f(id_t));
            r2_stats_T(id_f) = r2(nt).stats.T(id_t(1));
            r2_DS{id_f} = r2(nt).DS{id_t(1)};
            id_f = id_f+1;
            t_tmp = t_tmp+l_id_t;
        end
        clear id_f t_tmp id_t l_id_t l_T
        
        r2(nt).stats.i = r2_stats_i;
        r2(nt).stats.j = r2_stats_j;
        r2(nt).stats.x = r2_stats_x;
        r2(nt).stats.y = r2_stats_y;
        r2(nt).stats.z = r2_stats_z;
        r2(nt).stats.doy = r2_stats_doy;
        r2(nt).stats.yr = r2_stats_yr;
        r2(nt).stats.A = r2_stats_A;
        r2(nt).stats.trust = r2_stats_trust;
        r2(nt).stats.H = r2_stats_H;
        r2(nt).stats.V = r2_stats_V;
        pxls_Amax = max(r2(nt).stats.A/900);
        r2(nt).pct.Lmask_i = r2_pct_Lmask_i/pxls_Amax;
        r2(nt).pct.Wmask_i = r2_pct_Wmask_i/pxls_Amax;
        r2(nt).pct.CSmask_i = r2_pct_CSmask_i/pxls_Amax;
        r2(nt).pct.Smask_i = r2_pct_Smask_i/pxls_Amax;
        r2(nt).pct.Cmask_i = r2_pct_Cmask_i/pxls_Amax;
        r2(nt).pct.Omask_i = r2_pct_Omask_i/pxls_Amax;
        r2(nt).pct.Res_i = r2_pct_Res_i/pxls_Amax;
        r2(nt).pct.Res_f = r2_pct_Res_f/pxls_Amax;
        r2(nt).stats.T = r2_stats_T;
        r2(nt).DS = r2_DS';
        clear r2_* pxls_Amax

        if savedata==1
            T_T_1 = r2(nt).stats.T-[r2(nt).stats.T(1),r2(nt).stats.T(1:end-1)];
            IDW0 = find(T_T_1>180);
            r2A = r2(nt).stats.A;
            r2trust =r2(nt).stats.trust;
            r2H = r2(nt).stats.H;
            r2V = r2(nt).stats.V;
            r2yr = r2(nt).stats.yr;
            r2doy = r2(nt).stats.doy;
            r2T = r2(nt).stats.T;
            r2DS = r2(nt).DS;
            for i = 2:length(r2(nt).stats.V)-1
                if abs(r2(nt).stats.V(i)-mean([r2(nt).stats.V(i-1) r2(nt).stats.V(i+1)])) > ...
                        min([r2(nt).stats.V(i) mean([r2(nt).stats.V(i-1) r2(nt).stats.V(i+1)])]) && ...
                        T_T_1(i+1)+T_T_1(i)<=40% anomaly: Diff in V too big while Diff in T <=40
                    r2V(i) = NaN;
                end
            end
            clear i
            
            r3(nt).stats.doy = r2doy;
            r3(nt).stats.yr = r2yr;
            r3(nt).stats.H = r2H;
            r3(nt).stats.A = r2A;
            r3(nt).stats.V = r2V;
            r3(nt).stats.trust = r2trust;
            r3(nt).stats.T = r2T;
            r3(nt).DS = r2DS;
            
            i = 0;
            for idw0 = IDW0% when more than 6 months separate measurements
                r2A = [r2A(1:idw0+i-1),NaN,r2A(idw0+i:end)];
                r2trust = [r2trust(1:idw0+i-1),NaN,r2trust(idw0+i:end)];
                r2H = [r2H(1:idw0+i-1),NaN,r2H(idw0+i:end)];
                r2V = [r2V(1:idw0+i-1),NaN,r2V(idw0+i:end)];
                r2yr = [r2yr(1:idw0+i-1),NaN,r2yr(idw0+i:end)];
                r2doy = [r2doy(1:idw0+i-1),NaN,r2doy(idw0+i:end)];
                r2T = [r2T(1:idw0+i-1),NaN,r2T(idw0+i:end)];
                r2DS = [r2DS(1:idw0+i-1);NaN;r2DS(idw0+i:end)];
                i = i+1;
            end
            clear i idw0 IDW0 T_T_1
            
            lyr = max(r2yr)-min(r2yr);
            yearly_data_yr = cell(lyr,1);
            yearly_data = NaN(lyr,2);
            for y = min(r2yr):max(r2yr)-1
                yearly_data_yr{y-min(r2yr)+1} = [num2str(y),'-',num2str(y+1)];
                % water year
                if nnz((r2yr==y).*(r2doy>=278)+(r2yr==y+1).*(r2doy<278))
                    yearly_data(y-min(r2yr)+1,:) = ...
                        [min(r2V((r2yr==y).*(r2doy>=278)+(r2yr==y+1).*(r2doy<278)==1)) ...
                        max(r2V((r2yr==y).*(r2doy>=278)+(r2yr==y+1).*(r2doy<278)==1))];
                end
            end
            clear y

            lw = length(r2H);
            % saving storage variations into an Excel file
            xlswrite([dir_region,'dams_data',time,'.xlsx'],r2DS,dams{IMtoDAM(nt),1},['A2:A',num2str(lw+1)])
            xlswrite([dir_region,'dams_data',time,'.xlsx'],[r2trust' r2H' r2A'/1e6 r2V'], ...
                dams{IMtoDAM(nt),1},['B2:E',num2str(lw+1)])
            xlswrite([dir_region,'dams_data',time,'.xlsx'],yearly_data_yr,dams{IMtoDAM(nt),1},['G2:G',num2str(lyr+1)])
            xlswrite([dir_region,'dams_data',time,'.xlsx'],yearly_data,dams{IMtoDAM(nt),1},['H2:I',num2str(lyr+1)])
            clear r2A r2trust r2H r2V r2yr r2doy r2DS r2T lw lyr yearly_data*
        end
    end
end
clear nt I dir_dams progress*

if savedata==1
    % data saving:
    % Res_tot ~ all detected reservoirs merged in a single raster
    % img ~ statistics of all acquisition dates
    % r ~ statistics of all reservoirs
    fprintf('\n-> Saving global statistics\n')
    save(['statistics_2',time],'h_res','stats_res','img2','r2','DateString','r3')
end

%=========================================================================