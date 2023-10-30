% driver_plot_recent_HFR_vectors_fn_v1_rottnest.m
% driver to plot latest 24 hours of HFR vector map data
% yasha hetzel
% 2023-06-08

% for rottnest
addpath(genpath('/home2/yasha/matlab/m-files/yashatools/'))
addpath(genpath('/home2/yasha/matlab/m-files/yasha/m_map'))
addpath(genpath('/home2/yasha/matlab/m-files/altmany-export_fig-76c775b'));

%%
% inputs 
nodes={'NWA','CORL','TURQ','ROT','SAG','NEWC','COF'};
date_wind= 1 %days 1 hour = 1/24
input_directory='/home2/yasha/radar/data/latest24'
whereput='/home2/yasha/radar/data/latest24'; 

[SUCCESS,MESSAGE,MESSAGEID] =mkdir(whereput);

% remove old
delete([whereput filesep '*recent*.pdf']);
delete([whereput filesep '*recent*.png']);


for i=1:length(nodes)
    node=nodes{i};
    disp(['Plotting.......' node]);
plot_recent_HFR_vectors_fn_v1_rottnest(input_directory,node,date_wind, whereput);
end



%% combine the pdfs (using matlab - retains vector graphics)
% 
% delete([whereput filesep 'combined_recent_HFR.pdf']);
% clear flist
% d=dir([whereput filesep '0*recent*pdf']);for i=1:length(d);flist{i}=[whereput filesep d(i).name];end
% % flist
% append_pdfs([whereput filesep 'combined_recent_HFR.pdf'],flist{:})

%% or do it with bash imagemagick

% 
% % use pdfs (slower) - need to use this on rottenst
% convert -density 150 0*recent.pdf combined_recent_HFR.pdf

% % use pngs (best for easy quick)
% magick 0*recent.png combined_recent_HFR.pdf (fast)
% 
% or 

% keyboard

% !convert -density 150 /home2/yasha/radar/data/latest24/0*recent*v2.png /home2/yasha/radar/data/latest24/combined_recent_HFR_v2.pdf
disp('DONE creating figures... now combine in terminal')


%% DO HOURLY NOW

disp('Do latest hourly map now')

% inputs 
nodes={'NWA','CORL','TURQ','ROT','SAG','NEWC','COF'};
date_wind= 1/24 %days 1 hour = 1/24
input_directory='/home2/yasha/radar/data/latest24'
whereput='/home2/yasha/radar/data/latest1hr'; 

[SUCCESS,MESSAGE,MESSAGEID] =mkdir(whereput);

% remove old
delete([whereput filesep '*recent*.pdf']);
delete([whereput filesep '*recent*.png']);


for i=1:length(nodes)
    node=nodes{i};
    disp(['Plotting.......' node]);
plot_recent_HFR_vectors_fn_v1_rottnest(input_directory,node,date_wind, whereput);
end


%% DO 6 HOUR NOW

disp('Do latest 6 hour map now')

% inputs 
nodes={'NWA','CORL','TURQ','ROT','SAG','NEWC','COF'};
date_wind= 6/24 %days 1 hour = 1/24
input_directory='/home2/yasha/radar/data/latest24'
whereput='/home2/yasha/radar/data/latest6hr'; 

[SUCCESS,MESSAGE,MESSAGEID] =mkdir(whereput);

% remove old
delete([whereput filesep '*recent*.pdf']);
delete([whereput filesep '*recent*.png']);


for i=1:length(nodes)
    node=nodes{i};
    disp(['Plotting.......' node]);
plot_recent_HFR_vectors_fn_v1_rottnest(input_directory,node,date_wind, whereput);
end

disp('Done with 6 hourly...')