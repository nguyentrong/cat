function hypercolumn_v2%(filename)
% Version infomation:
% Weighting the matrix in the connectivity function

    %% close and clear all
    close all; 
    clc;
    warning off;  
    addpath('input');
    
    disp('Starting...');
      
    %% Define global variables
    global  x y nmass te ti he hi smax e0 r0 v0 n sz neqn tstart tend ...
        d p path signal_cut T
     
    %% Load initial data file  
    % te=p(1); ti=p(2); he=p(3); hi=p(4); 
    % smax=p(5); d=p(6); x=p(7); y= p(8); sz = p(9);
    % c1s=p(10); c1n=p(11); c1l=p(12);
    % c2s=p(13); c1n=p(14); c1l=p(15);
    disp('Loading initial parameter...');
    %load(['param/' filename '.mat']);
    %load([filename '.mat']);
    load p1x1.mat
    disp(['# te=' num2str(p(1)) 'ms; ti=' num2str(p(2)) 'ms; he=' num2str(p(3)) 'mV; hi=' num2str(p(4)) 'mV;']);
    disp(['# smax=' num2str(p(5)) '; d=' num2str(p(6)) '; x=' num2str(p(7)) '; y=' num2str(p(8)) '; sz=' num2str(p(9)) ]);
    disp(['# c1s=' num2str(p(10)) '; c1n=' num2str(p(11),'%.3f') '; c1l=' num2str(p(12),'%.3f')]);
    disp(['# c2s=' num2str(p(13)) '; c2n=' num2str(p(14),'%.3f') '; c2l=' num2str(p(15),'%.3f')]);
    disp(['# tstart=' num2str(p(16)) '; tend=' num2str(p(17))]);
       
    %% Parameters
    % synaptic and dendrite parameter
    te       = p(1);% 10  % ms
    ti       = p(2);% 20; % ms
    he       = p(3);% 3.25;  % mV
    hi       = p(4);% -22;% mV
    smax     = p(5);% 10;
    d        = p(6);
    % parameter of sigmoid function
    e0       = 2.5e-3;% 1/ms
    r0       = 0.56;  % 1/mV
    v0       = 6;     % mV
    % others
    x        = p(7);% 31 field size
    y        = p(8);% 21 field size
    n        = x*y;
    nmass    = n*3*8;
    neqn     = 2*nmass;
    sz       = p(9); % 0.1; % ms
    % integration time
    tstart   = p(16);   % ms
    tend     = p(17); % ms
    t        = tstart:sz:tend;
    signal_cut = 0;
    T        = [];
       
    % Result path
    path     = 'output';
    % Load map    
    %load(['input/map' num2str(x) 'x' num2str(y) '.mat'],'map');
    %map = [map{1}(:) map{2}(:) map{3}(:) map{4}(:) map{5}(:) map{6}(:) map{7}(:) map{8}(:) ];
    load input/cmap8OMs.mat   % cmap
    
    %% Solving the ode system
    [OM] = VC_ode();

    %% Postprocessing and plot and save
    %postproc(OM);
    save tmp_om.mat OM
    %OM = OM(:,end-2000:end);
    %OM = OM - repmat(mean(OM,2),1,size(OM,2));
    for i = 1:8
       subplot 421; hold on; plot(OM(i,:),'color',cmap(i,:)/255); axis tight; title('Pyramidal');
       subplot 423; hold on; plot(OM(i+8,:),'color',cmap(i,:)/255);axis tight; title('Excitatory');
       subplot 425; hold on; plot(OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('Inhibitory');
       subplot 427; hold on; plot(OM(i,:)+OM(i+8,:)+OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('EEG');
    end
    OM = OM(:,end-4000:end);
    for i = 1:8
       subplot 422; hold on; plot(OM(i,:),'color',cmap(i,:)/255); axis tight; title('Pyramidal');
       subplot 424; hold on; plot(OM(i+8,:),'color',cmap(i,:)/255);axis tight; title('Excitatory');
       subplot 426; hold on; plot(OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('Inhibitory');
       subplot 428; hold on; plot(OM(i,:)+OM(i+8,:)+OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('EEG');
    end    
    
    disp('Finished!');
end

function [OM] = VC_ode()
%% ODE solver using Euler method with constant stepsize 

global te ti he hi e0 r0 v0 nmass tstart sz tend d p smax signal_cut

% The connectivity and the weighting matrix
disp('Loading the connectivity matrix...');
[W] = localConnectivity();
W = W/smax^2;
%W = wtM*W;
%clear wtM
    
% modify the connectivity
ke = he/te; ki = hi/ti;
K  = [repmat(ke,2*nmass/3,1); repmat(ki,nmass/3,1)];
T1 = [repmat(-2/te,2*nmass/3,1); repmat(-2/ti,nmass/3,1)];
T2 = [repmat(-1/te^2,2*nmass/3,1); repmat(-1/ti^2,nmass/3,1)];

trunk = 1000;
n_trunk = ceil(tend/trunk/sz);

%%
%% Solver main loop
%%
% wtM = ones(nmass,1);
OM = [];  % Orientation map
Y0 = [zeros(nmass,1); zeros(nmass,1)];% Initial condition
Y0(1) = 100; Y0(nmass+1) = 100;
dv = zeros(nmass,1);
du = dv;

fprintf('Solver main loop...\n');

% additional info
trunkidx = 1; 
tic;
fprintf('[%.2f %.2f %.2f %.2f](%d):\n',te,ti,he,hi,n_trunk);

for idx = 1:length(tstart:sz:tend)      
    % Update half of Y0   
    u = Y0(1:nmass); 
    v = Y0(nmass+1:2*nmass);
   
    % sigmoid function of membrane potentials    
    S = 2*e0 ./ (1+exp(r0*(v0-u))); 
    
    % Weight S by the quotient h/tau  
    S = W*S;
       
    % ODEs 
    du = v;
    dv = T1.*v + T2.*u + K.*S;          
    
    % no noise
    % dY = [du; dv]
    Y1 = Y0 + sz * [du; dv];   
 
    % Update Y
    Y0 = Y1;
       
    % Save result in each step   
    if ~mod(idx,d) && idx > signal_cut 
        %OM = [OM single(det_om(u))];
        %OM = [OM u(1:end/3)+u(end/3+1:2*end/3)+u(2*end/3+1:end)];
        OM = [OM u];
    end
      
    if ~mod(idx,trunk)
       % display the informations
       estTime = toc;
       fprintf('%.3d(%.1fs) ',trunkidx,estTime);
       if ~mod(trunkidx,5); fprintf('\n'); end;
       trunkidx = trunkidx + 1;
       tic;
    end

end

%% End of solver loops
    
end

function postproc(OM)
    global x y te ti he hi smax path p 
   
    disp('Postprocessing...');
  
    currentTime = clock;
    currentTime = [num2str(currentTime(1)) '-'  num2str(currentTime(2)) '-'  num2str(currentTime(3)) '-' ...
                   num2str(currentTime(4)) '-'  num2str(currentTime(5)) '-'  num2str(currentTime(6),'%.0f')];      
    fname = ['Spontan_' num2str(x) 'x' num2str(y) '_te=' num2str(te) ...
        '_ti=' num2str(ti) '_he=' num2str(he) '_hi=' num2str(hi) '_smax=' num2str(smax) ...
        '_c1s=' num2str(p(10)) '_c1n=' num2str(p(11)) '_c1l=' num2str(p(12)) ...
        '_c2s=' num2str(p(13)) '_c2n=' num2str(p(14)) '_c2l=' num2str(p(15)) ];%...       
        %'_' currentTime ];
        
    load input/map31x21.mat map
    map = [map{1}(:) map{2}(:) map{3}(:) map{4}(:) map{5}(:) map{6}(:) map{7}(:) map{8}(:)];
    T = [];
    for ii = 1:size(OM,2) 
        om = OM(:,ii);
        tmp = [om map];
        [C] = corrcoef(tmp);
        [C Cidx] = sort(abs(C(1,2:end)),2,'descend'); 
        T = [T; C Cidx];
    end
    
    save([path '/mat/' fname '.mat'],'p','OM','T');   
      
    load input/cmap.mat    
    scrsz = get(0,'ScreenSize');
    h = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],'visible','off');    
    for ii = 1:1:size(T,1); 
        hold on; subplot 211; plot(ii,T(ii,1),'s','Color',cmap(T(ii,9),:)); 
    end
    axis tight;
    
    % Plot the distribution
    distr = []; 
    for ii = 1:8; 
        tmp = T(T(:,9) == ii);
        distr = [distr size(tmp,1)];  
    end    
    load input/Fig3KenetData.mat relOcc
    
    % square absolute errors
    err= distr/max(distr) - relOcc(:,2)';
    err = sqrt(sum(err.*err));    
    
    data = [ceil(relOcc(:,1)) relOcc(:,2) distr(:)/max(distr(:))];
    subplot 212; plot(data(:,1)-1,data(:,2),'s-','MarkerSize',12,'MarkerFaceColor','r','LineWidth',3)
    hold on;  subplot 212; plot(data(:,1)-1,data(:,3),'o-','MarkerSize',12,'MarkerFaceColor','g','LineWidth',3), hold off;   
    legend({'Data','Model'});  set(gca,'XTick',0:7,'XTickLabel',{'0','22.5','45','67.5','90','112.5','135','157.5'});
    ylabel('Rel. Occ.');   xlabel('Angle (degree)');   grid on;   ylim([0 1.1]);xlim([0 7]);    
    saveas(gcf,[path '/png/' fname '.png'], 'png');
    
    save([path '/mat/' fname '.mat'],'err','-append'); 
    
    close all;
    clearvars -global
end

function [out] = det_om(Y)
    global x y
    tmp = Y(1:end/3)+Y(end/3+1:2*end/3)+Y(2*end/3+1:end);
    out = reshape(tmp,x*y,8);
    out = sum(out,2);
end

function [W] = localConnectivity()
global p
    
%% The weighting matrix    
    c1s = p(10); c1n = p(11); c1l = p(12);
    c2s = p(13); c2n = p(14); c2l = p(15);
    
    wtMp = [c1s c1l c1l c1l c1l c1l c1l c1n;
            c1n c1s c1l c1l c1l c1l c1l c1l;
            c1l c1n c1s c1l c1l c1l c1l c1l;
            c1l c1l c1n c1s c1l c1l c1l c1l;
            c1l c1l c1l c1n c1s c1l c1l c1l;
            c1l c1l c1l c1l c1n c1s c1l c1l;
            c1l c1l c1l c1l c1l c1n c1s c1l;
            c1l c1l c1l c1l c1l c1l c1n c1s ];
        
     wtMi = [c2s c2l c2l c2l c2l c2l c2l c2n;
            c2n c2s c2l c2l c2l c2l c2l c2l;
            c2l c2n c2s c2l c2l c2l c2l c2l;
            c2l c2l c2n c2s c2l c2l c2l c2l;
            c2l c2l c2l c2n c2s c2l c2l c2l;
            c2l c2l c2l c2l c2n c2s c2l c2l;
            c2l c2l c2l c2l c2l c2n c2s c2l;
            c2l c2l c2l c2l c2l c2l c2n c2s ];  
        
%% The connectivity
    angle_min = 0;
    angle_max = 180;
    d_angle   = 22.5;
    angle     = angle_min:d_angle:angle_max-d_angle;
    
    sel_mat = 1500;%1650.36;
    
    % Orientation differences
    orientation_difference = [];
    for i = 1:8
        tmp = [];
       for j = 1:8
           tmp = [tmp abs(angle(i)-angle(j))];
       end
       orientation_difference = [orientation_difference; tmp];
    end
    orientation_difference(orientation_difference==157.5) = 22.5;
    orientation_difference(orientation_difference==135) = 45;
    orientation_difference(orientation_difference==112.5) = 67.5;    
    
    % There is only local connection between the masses
    diffAngleLocal = orientation_difference;
    diffAngleLocal(diffAngleLocal>40) = NaN;
    %Wpp_local = sel_mat*0.01.*exp(-0.0072*diffAngleLocal); Wpp_local(isnan(Wpp_local)) = 0; 
    Wpp_local = sel_mat*0.002.*exp(-0.0072*diffAngleLocal); Wpp_local(isnan(Wpp_local)) = 0; 

    % P->eN = local + distal
    Wpe_local = Wpp_local;
    % P->iN = local + distal
    Wpi_local = Wpp_local;
    % Wp = [Wpp Wpe Wpi];
    Wp = [wtMp.*Wpp_local wtMp.*Wpe_local wtMi.*Wpi_local];

    % Excitatory interneuron 
    % eN->P = local
    Wep = Wpp_local;
    % eN->eN = local
    Wee = Wpp_local;
    % eN->iN = local
    Wei = Wpp_local;
    % We = [Wep Wee Wei];
    We = [wtMp.*Wep wtMp.*Wee wtMi.*Wei];

    % Inhibitory interneuron 
    % iN->P = local
    diffAngleLocal = orientation_difference;
    diffAngleLocal(diffAngleLocal>60) = NaN;

    % Wip_local = sel_mat*0.08.*exp(-0.038*diffAngleLocal); 
    Wip_local = sel_mat*5.*exp(-0.038*diffAngleLocal);
    Wip_local(isnan(Wip_local)) = 0; 
    Wip = Wip_local; 

    % iN-eN = local
    Wie_local = Wip_local;

    % iN-iN = local
    Wii_local = sel_mat*0.04.*exp(-0.038*diffAngleLocal); 
    Wii_local(isnan(Wii_local)) = 0; 

    Wi = [wtMp.*Wip wtMp.*Wie_local wtMi.*Wii_local];         
   
    W = [Wp; We; Wi];
    W = W';
end

