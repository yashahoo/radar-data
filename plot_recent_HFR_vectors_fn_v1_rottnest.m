function plot_recent_HFR_vectors_fn_v1_rottnest(input_directory,node,date_wind, whereput)

% plot_recent_HFR_vectors_fn_v1_rottnest(input_directory,node,date_wind, whereput)

% % clear all;clc;close all;
% addpath ( '/Users/00068592/Documents/MATLAB/m-files/m_map');
% % addpath ( '/Users/00068592/Documents/MATLAB/m-files/curvvec');
% addpath('altmany-export_fig-76c775b') % this is contained in codar_wera_merge folder on github yasha-development branch
% addpath('/opt/homebrew/Cellar/ghostscript/10.0.0/bin')


% for rottnest
addpath(genpath('/home2/yasha/matlab/m-files/yashatools/'))
addpath(genpath('/home2/yasha/matlab/m-files/yasha/m_map/'))
addpath(genpath('/home2/yasha/matlab/m-files/altmany-export_fig-76c775b'));


% inputs to function

% date_wind= 1 %days 1 hour = 1/24
% input_directory='/Users/00068592/Dropbox/01.2017-2022/Shared_HFR_documents/HFR_data_processing/AODN_statistics'
% ,node,date_wind, whereput

% ----------------------------------------------------------------------- %
% -----------------General settings for ALL------------------------------ %
% ----------------------------------------------------------------------- %
save_plots='Y'
% whereput=['/home2/yasha/radar/data/latest24'];
% whereput_data='Output_averaged24'
get_date_range = 'N' % leave thsi as N for single plot (preferred)
% ----------------------------------------------------------------------- %
% -----------------Processing settings ALL------------------------------- %
% ----------------------------------------------------------------------- %
% velocity threshold
remove_fast='Y' % 'Y' to remove speeds > cutoff;  'N' do nothing
cutoff=2; % remove data with velocity > cutoff

apply_QC='Y' % 'Y' to use qc flag cutoff below 3) 'N' do not use QC flags at all
qlevel=3  %3 % remove data with Q flag > qlevel

flag_outliers='Y'

%  GDOP bounds is 'GDOP' above       orig  30-150 to 20-160 or 15-165
apply_GDOP='Y' % 'Y' or 'N'
GDOPmin=30 %15 %15 %30%15%15 %30 (outer areas)
GDOPmax=150 %135%140%65%35%165 %150 (coastal area)

plotback='none'
plot_labels='Y' % place names on map?
plot_drifters = 'N'
% -----------------------vector plot settings---------------------------- %
% settings for vectors with no background ( very sensitive)
% color the vectors according to speed? [only for plotback case 'none' above]
color_vecs='N' %'Y' % 'Y' or 'N'
shaftwid=.07;%.04;%.02
headlen=.25; %.15; %.1

min_num_vecs=5; % minimum number of vectorsto plot, will skip timestep if not more datapoints than this
clims=[ 0 .8]; % color caxis limits

% -----------------------scale vector settings---------------------------- %
% add scale vector
add_scale_vector='Y' % 'Y' or 'N'


%% ----------------------------------------------------------------------- %
% -----------------------define for each site---------------------------- %
% ----------------------------------------------------------------------- %

switch node
    
    case 'NWA'
        site_num='01'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 30;
        west=112.3; east=114.5; south=-23;north= -21; %big%
        % west=112.5; east=114.5; south=-22.6;north= -21.3; %zoom 2
        % west=112; east=114.5; south=-23;north= -20.5; %big wide
        vec_scale=100 %35
        % scale_pos=[114.2 -23.4]; %big NWA
        scale_pos=[112.5 -22.7]; %big NWA
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %2 NWA removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        % Mask data inside Exmouth Gulf
        mask_gulf='Y' % 'Y' or 'N'
        
        %       Place labels
        place(1).name='Jurabi (JTC) ';    place(1).location=[114.0999 -21.807 ];
        place(2).name='   Point Billie (PTB)';  place(2).location=[113.6903 -22.5428];
        place(3).name=' Exmouth'; place(3).location=[114.13 -21.954848];
        
        
        
        
    case 'CORL'
        site_num='02'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 0;
        west=113; east=115.5; south=-31;north= -28.5; %CORL
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        
        place(1).name='Green Head';  place(1).location=[114.96606514837579, -30.071444636915817];
        place(2).name='Dongara';    place(2).location=[114.92353125420638, -29.281149902121523];
        mask_gulf='N'; % only used for NWA
        
        
        scale_pos=[113.25 -30.7]; %
        
        
        
    case 'TURQ'
        site_num='03'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 30;
        west=113; east=116; south=-32.25;north= -29.5;% TURQ
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        place(1).name='Lancelin'; place(1).location=[115.32884071237285, -31.026203600884845];
        place(2).name='Green Head';  place(2).location=[114.96606514837579, -30.071444636915817];
        mask_gulf='N'; % only used for NWA
        
        scale_pos=[113.35 -32.0]; %
        
    case 'ROT'
        site_num='04'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 0;
        west=114.2; east=116; south=-32.5;north= -31; %ROT
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        place(1).name='Fremantle'; place(1).location=[115.74709564688033, -32.03129069711441];
        place(2).name='Guilderton';  place(2).location=[115.4905732617174, -31.342357642778055];
        mask_gulf='N'; % only used for NWA
        
        scale_pos=[115.33 -32.4]; 
        
    case 'SAG'
        site_num='05'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 0;
        %         west=134.0; east=135.5; south=-34;north= -36; %SAG
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        place(1).name='Cape Spencer'; place(1).location=[136.881316, -35.297302]
        place(2).name='Cape Wiles'; place(2).location=[135.684682,-34.944694];
        mask_gulf='N'; % only used for NWA
        
        scale_pos=[112.5 -22.7]; %big NWA
        
    case 'NEWC'
        site_num='06'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 0;
        west=151.0; east=154.0; south=-34.25;north= -32; %NEWC
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        place(1).name='Seal Rocks'; place(1).location=[152.537872, -32.440928];
        place(2).name='Red Head'; place(2).location=[151.7268, -33.0109167];
        place(3).name='Newcastle'; place(3).location=[151.770209, -32.930713];
        mask_gulf='N'; % only used for NWA
        
        scale_pos=[153.0 -34.2]; %
        
        
    case 'COF'
        site_num='07'
        ncref=[input_directory filesep node '_recent.nc'];
        site_mins= 30;
        west=152.9; east=154.25; south=-30.9;north= -29.7; %COF
        
        vec_scale=100 %50
        scale_speed=0.5;
        
        % Filter outliers before plotting?
        flag_outliers='Y'
        std_factor=3; %removes speed outliers > [std_factor]*standard deviations from median. default is 3, set very high to remove no data (now shoudl not remove if flag_outliers set to 'N'
        
        place(1).name='Red Rock'; place(1).location=[153.2311, -29.98388];
        place(2).name='Nambucca Heads'; place(2).location=[153.0111, -30.624166];
        place(3).name='Coffs Harbour';place(3).location=[153.130527,-30.312904];
        mask_gulf='N'; % only used for NWA
        
        scale_pos=[112.5 -22.7]; %big NWA
        
 
end
%% read the nc

try
    tim=ncread(ncref,'TIME');
    time=tim+datenum(1950,1,1); %TIME:units = "days since 1950-01-01 00:00:00 UTC" ;
    gU=permute(ncread(ncref,'UCUR'),[3,2,1]);
    gV=permute(ncread(ncref,'VCUR'),[3,2,1]);
    Uqual=permute(ncread(ncref,'UCUR_quality_control'),[3,2,1]);
    DOP=permute(ncread(ncref,'GDOP'),[2,1]);
    DOP3(1,:,:) = DOP;
    GDOP=repmat(DOP3,length(time),1,1);
    
    lat=ncread(ncref,'LATITUDE');
    lon=ncread(ncref,'LONGITUDE');
    
    
    if isvector(lon)
        glat=repmat(lat,1,length(lon));
        glon=repmat(lon',length(lat),1);
        
    else
        lat=permute(ncread(ncref,'LATITUDE'),[2,1]); % to be consistent  with NWA for otehr sites i need to add this permute even though logically it doesnt make sense
        lon=permute(ncread(ncref,'LONGITUDE'),[2,1]);
        glat=lat;
        glon=lon;
    end
    
    % get file version for naming and titles
    file_version=ncreadatt(ncref,'/','file_version');
    switch file_version
        case 'Level 1 - Quality Controlled data'
            FV='FV01'
        case 'Level 0 - Raw data with RT Quality Control'
            FV='FV00'
        case 'Level 0 - Real Time Quality Controlled data'
            FV='FV00'
        case 'Level 0 - RT Quality Controlled data'
            FV='FV00'
        case 'Level 0 - Raw data'
            FV='FV00'
    end
    
    
    % get map bounds if not given above
    if ~exist('west','var') & ~exist('east','var') & ~exist('north','var') & ~exist('south','var')
        
        west=min(lon(:))-.25;
        east=max(lon(:))+.25;
        south=min(lat(:))-.25;
        north=max(lat(:))+.25;
    end
    
    
    %%
    gspd=sqrt(gU.^2+gV.^2);
    
    % apply quality control flags and/or GDOP
    QCi=[];GDOPi=[]; badi=[];
    switch apply_QC
        case 'Y'
            switch apply_GDOP
                case 'N'
                    QCi=find(Uqual>qlevel);
                    gU(QCi)=NaN;
                    gV(QCi)=NaN;
                    gspd(QCi)=NaN;
                    disp('Applying QC but not GDOP')
                case 'Y'
                    QCi=find(Uqual>qlevel); GDOPi=find(GDOP>GDOPmax | GDOP<GDOPmin);
                    badi=intersect(QCi,GDOPi);
                    gU(badi)=NaN;
                    gV(badi)=NaN;
                    gspd(badi)=NaN;
                    disp('Applying QC and GDOP')
                    %                 deal with case where GDOP is narrower than 30-150
                    if (GDOPmax<150 | GDOPmin<30)
                        GDOPi=find(GDOP>GDOPmax | GDOP<GDOPmin);
                        gU(GDOPi)=NaN;
                        gV(GDOPi)=NaN;
                        gspd(GDOPi)=NaN;
                        disp('fixing case where GDOP is narrower than 30-150')
                    end
            end
            
        case 'N'
            qlevel=10; % if not using Quality flags, assign unrealistic QF to keep all data
            switch apply_GDOP
                case 'Y'
                    badi=find(GDOP>GDOPmax | GDOP<GDOPmin);
                    gU(badi)=NaN;
                    gV(badi)=NaN;
                    gspd(badi)=NaN;
                    disp('Applying  GDOP only')
                    
            end
            
    end
    
    % remove speeds above threshold
    switch remove_fast
        case 'Y'
            bad=find(gspd>cutoff);
            gU(bad)=nan;
            gV(bad)=nan;
            gspd(bad)=nan;
    end
    
    
    
    
    mkdir(whereput)
%     mkdir(whereput_data)
    
    % get string for proper labels and titles
    if date_wind<1 && date_wind>=0.0417
        date_wind_str=[sprintf('%02.0f',date_wind*24) '_hour_mean'];
        title_date_wind_str=[sprintf('%02.0f',date_wind*24) ' hour mean'];
    elseif date_wind>=1
        date_wind_str=[sprintf('%02.0f',date_wind) '_day_mean'];
        title_date_wind_str=[sprintf('%02.0f',date_wind) ' day mean'];
        
    elseif date_wind<0.0417
        date_wind_str=[sprintf('%02.0f',date_wind*24*60) '_minute_mean'];
        title_date_wind_str=[sprintf('%02.0f',date_wind*24*60) ' minute mean'];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% make  HFR averages
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timestart=time(end)-date_wind;
    timestop=time(end);
    tstepi=find(time>timestart);
    num_tsteps_available=length(tstepi)
    datestr(time(tstepi))
    
    
    
    switch get_date_range
        case 'Y'
            tsteps=length(date_range);
            umean=nan(length(tsteps),size(glon,1),size(glon,2));
            vmean=nan(length(tsteps),size(glon,1),size(glon,2));
            avgspd=nan(length(tsteps),size(glon,1),size(glon,2));
            hftime=nan(length(tsteps),1);
        case 'N'
            date_range=(timestart+timestop)/2;
            tsteps=1;
            umean=nan(1,size(glon,1),size(glon,2));
            vmean=nan(1,size(glon,1),size(glon,2));
            avgspd=nan(1,size(glon,1),size(glon,2));
            hftime=nan(length(tsteps),1);
            
    end
    
    % get ready
    t=time;
    t(isnan(t))=[];
    
    for i=1:length(date_range)
        
        switch get_date_range
            case 'Y'
                timestart=date_range(i)-date_wind./2;
                timestop=date_range(i)+date_wind./2;
                disp([num2str(i) '  of  ' num2str(length(date_range)) '  ' date_wind_str])
            case 'N'
                disp([num2str(i) '.  Getting single average over ' num2str(timestop-timestart)  '  days using timestart and timestop'])
        end
        
        ti=find(t>timestart); % & t<timestop
        umean(i,:,:)=squeeze(mean(gU(ti,:,:),1,'omitnan'));
        vmean(i,:,:)=squeeze(mean(gV(ti,:,:),1,'omitnan'));
        avgspd(i,:,:)=sqrt(umean(i,:,:).^2+vmean(i,:,:).^2);
%         hftime(i)=mean(t(ti)); % for mid time
        hftime(i)=t(max(ti));
        
    end
    
    % get rid of empty timesteps
    if sum(isnan(hftime))>0
        umean(isnan(hftime),:,:)=[];
        vmean(isnan(hftime),:,:)=[];
        avgspd(isnan(hftime),:,:)=[];
        hftime(isnan(hftime))=[];
        disp(['Removed ' num2str(sum(isnan(hftime))) ' empty timesteps'])
    end
    
    
    % rename variables to fit code
    LO=glon;
    LA=glat;
    %%
    
    %% get sst data to plot if available
    switch plotback
        
        case 'sst'
            disp('DOWNLOADING IMOS L3S sst multisensor data....')
            addpath('/Users/00068592/Documents/MATLAB/m-files/yashatools/SST_L3S_SSTAARS_toolbox_yh')
            %      [outfname]=get_IMOS_L3S_openddap_AODN_v8loop_fn(years,months,days,range,mode,plotit,south,north,east,west,sitename,use_multisensor)
            
            str2num(datestr(hftime(1),'yyyymmdd'))
            %                             [outfname]=get_IMOS_L3S_openddap_AODN_v8loop_fn(str2num(datestr(hftime(1),'yyyy')),str2num(datestr(hftime(1),'mm')),[str2num(datestr(hftime(1),'dd')): str2num(datestr(hftime(end),'dd'))],sst_range,'night','N',south,north,east,west,node,'Y');
            
            % %                             if ~exist('sst_outfname','var')
            [sst_outfname]=get_IMOS_L3S_openddap_AODN_v9loop_fn(str2num(datestr(hftime(1),'yyyy')),str2num(datestr(hftime(1),'mm')),[str2num(datestr(hftime(1),'dd')): str2num(datestr(hftime(end),'dd'))],sst_range,'night','N',south,north,east,west,node,'Y',4);
            % %                             else
            % %                             disp('SST data file available locally... not downloading!')
            % %                             end
            
            sst_data_available='Y';
            try
                sst=load(sst_outfname);
            catch
                disp('No sst data available')
                sst_data_available='N';
            end
            
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot with m_map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %set the projection
    m_proj('mercator','lat',[south north],'lon',[west east]);
    
    f1=figure;
    f1.Color='w';
    f1.Position=[680 80 900 900]; % change figure size here
    f1.Units='inches'; % yh test
    f1.Position=[14.8750 3.4722 8 8] %yh test
    
    
    la=LA;lo=LO;
    
    % gridLON,gridLAT,uu,vv
    
    % make equivalent 2d version of lons and lats corresponding to
    gridLON=reshape(glon,size(glon,1)*size(glon,2),1);
    gridLAT=reshape(glat,size(glat,1)*size(glat,2),1);
    mask=load('ExGulf_mask'); % this is polygon for exmouth gulf
    
    for ind=1:length(hftime) %%60 %1000:1002%:length(dailytime) %50:60
        
        %  get subset for each timestep
        speed_dat=sqrt(squeeze(umean(ind,:,:)).^2+squeeze(vmean(ind,:,:)).^2);
        datu=squeeze(umean(ind,:,:));
        datv=squeeze(vmean(ind,:,:));
        sp=reshape(speed_dat,size(datu,1)*size(datu,2),1);
        uu=reshape(datu,size(datu,1)*size(datu,2),1);
        vv=reshape(datv,size(datu,1)*size(datu,2),1);
        maxmag=nanmax(nanmax(speed_dat))
        
        
        
        
        
        %---------------------mask exmouth gulf points----------------------------%
        switch mask_gulf
            case 'Y'
                IN = inpolygon(gridLON,gridLAT,mask.X,mask.Y) ;
                sp(IN)=NaN;
                uu(IN)=NaN;
                vv(IN)=NaN;
        end
        %--------------------- add scale vector-------------------------------%
        switch add_scale_vector
            case 'Y'
                uu=[uu;scale_speed];
                vv=[vv;0];
                sp=[sp;NaN];
                
                if length(gridLON)< prod(size(glon))+1
                    gridLON=[gridLON;scale_pos(1)];
                    gridLAT=[gridLAT;scale_pos(2)];
                end
                
                sh=m_text(scale_pos(1),scale_pos(2)+.01,[sprintf('%0.2f',scale_speed) ' m s^-^1']);
                sh.FontSize=10;sh.VerticalAlignment='bottom';
                sh.HorizontalAlignment='center';
                
        end
        
        %---------------------flag outliers as nans-------------------------------%
        switch flag_outliers
            case 'Y'
                %   lat/lon grid
                gbadi=find(speed_dat>(median(speed_dat(:),'omitnan')  +std_factor*std(speed_dat(:),'omitnan')));
                speed_dat(gbadi)=NaN; datu(gbadi)=NaN; datv(gbadi)=NaN;
                %   position format
                badi=find(sp>(median(sp,'omitnan')+std_factor*std(sp,'omitnan')));
                sp(badi)=nan;uu(badi)=nan;vv(badi)=nan;
        end
        %-------------------------------------------------------------------------%
        
        if sum(sum(~isnan(uu)))>min_num_vecs; % skip the timesteps with no (or almost no) data %sum(sum(~isnan(datu)))>min_num_vecs;
            hold on
            %---------------------plot background pcolor------------------------------%
            switch plotback
                case 'speed'
                    hp=m_pcolor(lo,la,speed_dat); % plot speed background
                    shading flat
                    [hq, HT]=m_vec(2,LO+.03,LA-.03,datu,datv,'k','shaftwidth', 0.5, 'headlength', 4,'edgeclip','on','centered','yes'); %[HP, HT]
% [hq, HT]=m_vec(2,LO+.03,LA-.03,datu,datv,'k','shaftwidth', 0.5, 'headlength', 4,'centered','yes'); %[HP, HT]
                    %                                                      hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,'k','shaftwidth', shaftwid, 'headlength', headlen,'edgeclip','on','centered','yes');
                case 'sst'
                    
                    switch sst_data_available
                        case 'Y'
                            %                 get closest sst data if available
                            sstime=sst.L3S_time
                            datestr(hftime(ind))
                            dst=sstime-hftime(ind);
                            [mdst mdsi]=min(abs(dst)); % clostest satellite time to radar time. min dt and corresponding index
                            if mdst<3
                                sdat=squeeze(sst.L3S_sst(mdsi,:,:));
                                slon=sst.L3S_lon;
                                slat=sst.L3S_lat;
                                
                                
                                hp=m_pcolor(slon,slat,sdat); % plot sst background
                                shading flat
                                [hq, HT]=m_vec(1,LO+.03,LA-.03,datu,datv,'k','shaftwidth', 0.5, 'headlength', 4,'edgeclip','on','centered','yes'); %[HP, HT]
% [hq, HT]=m_vec(1,LO+.03,LA-.03,datu,datv,'k','shaftwidth', 0.5, 'headlength', 4,'centered','yes'); %[HP, HT]
                                %                                                                 hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,'k','shaftwidth', shaftwid, 'headlength', headlen,'edgeclip','on','centered','yes');
                            end
                        case 'N'
                            sst_data_available
                            
                            hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,'k','shaftwidth', shaftwid, 'headlength', headlen,'edgeclip','on','centered','yes');
% hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,'k','shaftwidth', shaftwid, 'headlength', headlen,'centered','yes');
                            %                       [hq, HT]=m_vec(2,LO+.03,LA-.03,datu,datv,'k','shaftwidth', 0.5, 'headlength', 4,'edgeclip','on','centered','yes');
                    end
                    
                case 'none'
                    disp('Not plotting pcolor background')
                    
                    %---------------------------plot vectors----------------------------------%
                    switch color_vecs
                        case 'Y'
                            hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,sp,'shaftwidth', shaftwid, 'headlength', headlen,'edgeclip','on','centered','yes');
                            
% hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,sp,'shaftwidth', shaftwid, 'headlength', headlen,'centered','yes');
                            
                            % % % % % % %                  hq = m_curvvec(LO,LA,datu,datv)
                            % % % % % % %                 hq = m_curvvec(LO,LA,datu,datv,'color',jet,'minmag',5,...
                            % % % % % % %             'maxmag',50,'thin',3)
                            
                            cbh=colorbar;
                            caxis(clims)
                            ylabel(cbh,'Velocity (m s^-^1)')
                        case 'N'
                            
% % % %                             %  hq = m_vec(vec_scale, gridLON,gridLAT,uu,vv,'k','shaftwidth',shaftwid, 'headlength', headlen,'edgeclip','on','centered','yes');
%                             hq = m_vec(vec_scale,gridLON,gridLAT,uu,vv,'k','shaftwidth',shaftwid,'headlength',headlen,'edgeclip','on','centered','yes') %normal good plotting way
%                             keyboard
                             hq = m_quiver(gridLON,gridLAT,uu,vv,'k') ;
% % % % %                             hq = m_vec(vec_scale,gridLON,gridLAT,uu,vv,'k','shaftwidth',shaftwid,'headlength',headlen,'centered','yes') %rottnest not have this version to edgeclip
% % % %                             %                                              hq = m_vec(5, gridLON(1:2:end,1:2:end),gridLAT(1:2:end,1:2:end),uu(1:2:end,1:2:end),vv(1:2:end,1:2:end),'k','edgeclip','on','centered','yes'); %COF
% % % %                             %                                               hq = m_vec(5, gridLON(1:3:end,1:3:end),gridLAT(1:3:end,1:3:end),uu(1:3:end,1:3:end),vv(1:3:end,1:3:end),'k','edgeclip','on','centered','yes'); % COF 2
                    end
            end
            
            
            %--------------plot depth contours (not implimented)----------------------%
            % % plot ROMS bathy
            % hc=m_contour(lon_rho,lat_rho,h,[0 100 ],'k','linewidth',2); %[0 10 20 50 100 200][50 100 500 1000 2000 3000]
            
            %------------------------plot map stuff----------------------------------%
            %plot the data (projected)
            m_usercoast('gshhs_h_AUSNZ','patch',[.6 .6 .6],'Linewidth', .25);
            %     % m_grid('box','fancy','tickdir','in');
            m_grid('box','on','tickdir','in','linestyle','none','fontsize',10);
            %     % caxis([0 1.5])
            colormap jet
            box on
            %         m_ruler([.7 .9],.95,2,'fontsize',10,'tickdir','in') % make a scale bar :-)
            %         m_ruler([.7 .95],.06,2,'fontsize',10,'tickdir','in') % make a scale bar :-) % for big region
            m_ruler([.1 .4],.06,2,'fontsize',10,'tickdir','in') % make a scale bar :-) % for small region
            % title(datestr(dailytime(ind),'dd-mmm-yyyy'))
            %     title(['HFR average (' datestr(t(ti(1)),'dd-mmm-yyyy') ' - ' datestr(t(ti(end)),'dd-mmm-yyyy') ')'])
%             keyboard
            switch get_date_range
                case 'Y'
                    switch apply_QC
                        case 'Y'
                            switch apply_GDOP
                                case 'Y'
%                                     title([FV ' QF<' num2str(qlevel+1) ' ' num2str(GDOPmin) ' <GDOP< ' num2str(GDOPmax) ' HFR ' title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                      title(['QF<' num2str(qlevel+1)  ' HFR ' title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])

                                case 'N'
                                    title([FV ' QF<' num2str(qlevel+1)  ' HFR ' title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
%                                     title([FV ' QF<' num2str(qlevel+1)  ' HFR ' title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                    
                            end
                        case 'N'
                            switch apply_GDOP
                                case 'Y'
                                    
                                    title([FV ' No QC flags '   num2str(GDOPmin) ' <GDOP< ' num2str(GDOPmax)  ' HFR '  title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                
                                case 'N'
                                    title([FV ' No QC flags HFR '  title_date_wind_str ' ending on  (' datestr(hftime(ind),'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                            end
                    end
                    
                case 'N'
                    switch apply_QC
                        case 'Y'
                            switch apply_GDOP
                                case 'Y'
%                                     title([FV ' QF<' num2str(qlevel+1) ' ' num2str(GDOPmin) ' <GDOP< ' num2str(GDOPmax) ' HFR average (' datestr(timestart,'dd-mmm-yyyy HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                    title(['QF<' num2str(qlevel+1) ' HFR average (' datestr(timestart,'dd-mmm-yyyy HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])

                                case 'N'
%                                     title([FV ' QF<' num2str(qlevel+1) ' HFR average (' datestr(timestart,'dd-mmm-yyyy  HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy  HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                      title(['QF<' num2str(qlevel+1) ' HFR average (' datestr(timestart,'dd-mmm-yyyy  HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy  HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                            end
                        case 'N'
                            switch apply_GDOP
                                case 'Y'
                                    title([FV ' No QC flags ' num2str(GDOPmin) ' <GDOP< ' num2str(GDOPmax) ' HFR average (' datestr(timestart,'dd-mmm-yyyy  HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy  HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                                case 'N'
                                    title([FV ' No QC flags HFR average (' datestr(timestart,'dd-mmm-yyyy  HH:MMz') ' - ' datestr(timestop,'dd-mmm-yyyy  HH:MMz') ') # timesteps= ' num2str(num_tsteps_available)])
                            end
                    end
            end
            
            %     end
            
% ----------------------------------------------------------------------- %    
%             add text -  time since latest data
% ----------------------------------------------------------------------- %
% IF DATE IS MORE THAN 1 day ago make it RED
days_since_last_data=((now-8/24)-hftime(ind));
% hax=m_text(east-0.4,north-0.1,[sprintf('%.1f',days_since_last_data*24) ' hours since last data']);
hax=m_text(west+0.8,north-0.1,[sprintf('%.1f',days_since_last_data*24) ' hours since last data']);
hax.HorizontalAlignment='center';
hax.FontSize=10;




if days_since_last_data>28/24
    tax=gca;
    tax.Title.Color='r';
    hax.Color='r';
end
% ----------------------------------------------------------------------- %    
            
            
            %     cbh=colorbar;
            switch plotback
                case 'speed'
                    cbh=colorbar;
                    caxis(clims)
                    ylabel(cbh,'Velocity (m s^-^1)')
                case 'sst'
                    switch sst_data_available
                        case 'Y'
                            cbh=colorbar;
                            caxis([csst(1) csst(2)])
                            %                                             caxis auto % usually looks good
                            ylabel(cbh,'SST (^oC)')
                    end
            end
            % pause(.5)
            
            %----------------option to plot text labels-------------------------------%
            switch plot_labels
                case 'Y'
                    for i=1:length(place)
                        ph(i)=m_scatter(place(i).location(1),place(i).location(2));
                        th(i)=m_text(place(i).location(1)+2/60,place(i).location(2),place(i).name)
                        th(i).FontSize=10;th(i).VerticalAlignment='Top'; th(i).HorizontalAlignment='left';
                        ph(i).Marker='o';ph(i).MarkerFaceColor='k';ph(i).MarkerEdgeColor='k';
                    end
%                     th(1).HorizontalAlignment='right';
% move text to left side if east coast
if strcmpi(node,'COF') | strcmpi(node,'NEWC')
        disp(node)
        delete(th)
        for i=1:length(place)
            th(i)=m_text(place(i).location(1)-2/60,place(i).location(2),place(i).name)
            th(i).FontSize=10;th(i).VerticalAlignment='Top'; th(i).HorizontalAlignment='right';
            th(i).HorizontalAlignment='right';
        end     
end


            end
            
            set(gca,'fontsize',10)
            
            
 %-----------------------plot drifter tracks-------------------------------%

 %  get centered on same range
 switch plot_drifters
     case 'Y'
         switch get_date_range
             case 'Y'
                 try
                     dhfi=find(dtim<hftime(ind) &  dtim>(hftime(ind)-drift_window_days));
                     
                     m_plot(dlon(dhfi),dlat(dhfi),'b')
                     m_plot(dlon(dhfi(end)),dlat(dhfi(end)),'r*')
                     m_text(dlon(dhfi(end)),dlat(dhfi(end)),datestr(dtim(dhfi(end))),'dd-mmm')
                 catch
                     disp('No drifter data in given window...')
                 end
                 
             case 'N'
                 
                 try
                     dhfi=find(  dtim>timestart  & dtim<timestop );
                     
                     m_plot(dlon(dhfi),dlat(dhfi),'b')
                     m_plot(dlon(dhfi(end)),dlat(dhfi(end)),'r*')
                     m_text(dlon(dhfi(end)),dlat(dhfi(end)),datestr(dtim(dhfi(end))),'dd-mmm')
                     
                 catch
                     disp('No drifter data in given window...')
                 end
         end
 end
 %---------------------------save figures----------------------------------%

% make simple names  
fname=[whereput filesep site_num '_' node '_' FV '_recent']

            
%             switch get_date_range
%                 case 'Y'
%                     switch apply_QC
%                         case 'Y'
%                             switch apply_GDOP
%                                 case 'Y'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(hftime(ind),'yyyymmdd_HHz_') date_wind_str '_qlevel_' num2str(qlevel) '_' num2str(GDOPmin) 'GDOP' num2str(GDOPmax)]
%                                 case 'N'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(hftime(ind),'yyyymmdd_HHz_') date_wind_str '_qlevel_' num2str(qlevel)]
%                             end
%                             
%                             
%                         case 'N'
%                             switch apply_GDOP
%                                 case 'Y'
%                                     fname=[whereput '/' FV '_'  node '_'  datestr(hftime(ind),'yyyymmdd_HHz_') date_wind_str '_NoQC_' num2str(GDOPmin) 'GDOP' num2str(GDOPmax)]
%                                 case 'N'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(hftime(ind),'yyyymmdd_HHz_') date_wind_str '_NoQC']
%                             end
%                     end
%                     
%                 case 'N'
%                     switch apply_QC
%                         case 'Y'
%                             switch apply_GDOP
%                                 case 'Y'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(timestart,'yyyymmdd') '-' datestr(timestop,'yyyymmdd') '_mean_qlevel_' num2str(qlevel) '_' num2str(GDOPmin) 'GDOP' num2str(GDOPmax)]
%                                 case 'N'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(timestart,'yyyymmdd') '-' datestr(timestop,'yyyymmdd') '_mean_qlevel_' num2str(qlevel)]
%                             end
%                             
%                             
%                         case 'N'
%                             switch apply_GDOP
%                                 case 'Y'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(timestart,'yyyymmdd') '-' datestr(timestop,'yyyymmdd') '_mean_NoQC' num2str(qlevel) '_' num2str(GDOPmin) 'GDOP' num2str(GDOPmax)]
%                                 case 'N'
%                                     fname=[whereput '/' FV '_' node '_'  datestr(timestart,'yyyymmdd') '-' datestr(timestop,'yyyymmdd') '_mean_NoQC' num2str(qlevel)]
%                             end
%                     end
%             end
            
            % %                             pause
            switch save_plots
                case 'Y'
set(gca,'SortMethod','ChildOrder')
% %                     export_fig([fname '_test1'],'-pdf','-painters'); %'-c[70 80 100 80]' %,'-transparent'
%                     export_fig([fname '_v1'],'-png','-r300');
                    
%                     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
%                     print([fname '_9x6in'],
                    
%                     print([fname '_test2'],'-depsc','-painters');
                    print([fname '_v2'],'-dpng','-r300');

                    % pause(.2)
                    % % delete(hp); delete(hq);delete(hc); delete(hr1);  delete(hr2);
                    % % close;
                    
%                     screen2jpeg([fname '_screen'])
% % % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
% % % export_fig(fname,'-png','-r300');
                    
                case 'N'
                    disp('NOT SAVING PLOTS')
            end
            %-------------------------------------------------------------------------%
                                        pause(5)
            %                                                                     pause
            clf
            
        else
            disp(['No data to plot on: ' datestr(hftime(ind))]);
            disp(['Done with ' num2str(ind) ' of ' num2str(length(hftime))]);
            clf % yh
        end
    end
    
    close all % yh
    
catch
    
%     keyboard
%     lasterr
    disp([ncref ' not available...trying next...'])
    pause(2)
    
end
%             end
%         end
%     end
% end


end