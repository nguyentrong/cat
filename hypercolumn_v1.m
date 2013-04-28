function hypercolumn%(filename)
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
        map d p path signal_cut T
     
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
    signal_cut = 100/sz;
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
    OM = OM(:,2000:end);
    %OM = OM - repmat(mean(OM,2),1,size(OM,2));
    for i = 1:8
       subplot 411; hold on; plot(OM(i,:),'color',cmap(i,:)/255); axis tight; title('Pyramidal');
       subplot 412; hold on; plot(OM(i+8,:),'color',cmap(i,:)/255);axis tight; title('Excitatory');
       subplot 413; hold on; plot(OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('Inhibitory');
       subplot 414; hold on; plot(OM(i,:)+OM(i+8,:)+OM(i+16,:),'color',cmap(i,:)/255); axis tight; title('EEG');
    end
    
    
    disp('Finished!');
end

function [OM] = VC_ode()
%% ODE solver using Euler method with constant stepsize 

global te ti he hi e0 r0 v0 nmass tstart sz tend d p smax signal_cut

% The connectivity and the weighting matrix
disp('Loading the connectivity matrix...');
[W wtM] = localConnectivity();
W = W/smax^2;
%W = wtM*W;
%clear wtM
    
% modify the connectivity
ke = he/te; ki = hi/ti;
K  = [repmat(ke,2*nmass/3,1); repmat(ki,nmass/3,1)];
T1 = [repmat(-2/te,2*nmass/3,1); repmat(-2/ti,nmass/3,1)];
T2 = [repmat(-1/te^2,2*nmass/3,1); repmat(-1/ti^2,nmass/3,1)];

trunk = 50;
n_trunk = ceil(tend/trunk/sz);

%%
%% Solver main loop
%%
% wtM = ones(nmass,1);
OM = [];  % Orientation map
Y0 = [-70*ones(nmass,1); zeros(nmass,1)];% Initial condition
dv = zeros(nmass,1);
du = dv;

fprintf('Solver main loop...\n');

% additional info
trunkidx = 1; 
tic;
c1s = p(10); c1n = p(11); c1l = p(12);
c2s = p(13); c2n = p(14); c2l = p(15);
fprintf('[%.2f %.2f %.2f %.2f](%d):\n',te,ti,he,hi,n_trunk);

for idx = 1:length(tstart:sz:tend)      
    % Update half of Y0   
    u = Y0(1:nmass); 
    v = Y0(nmass+1:2*nmass);
   
    % sigmoid function of membrane potentials    
    S = 2*e0 ./ (1+exp(r0*(v0-u))); 
    
    % Weight S by the quotient h/tau  
    S = W*S;
    if idx > signal_cut
        tmp = u(1:nmass/3)+u(nmass/3+1:2*nmass/3)+u(2*nmass/3+1:nmass);
        [maxC_val maxC_idx] = max(tmp);
        wtM_ex = repmat(c1l,nmass/3,1); % all states = last state = min
        wtM_in = repmat(c2l,nmass/3,1); % all states = last state = max            
        wtM_ex(maxC_idx) = c1s; % current state = 1
        wtM_in(maxC_idx) = c2s; % current state = 1
        if maxC_idx < 8            
            wtM_ex(maxC_idx+1) = c1n; % next state = max
            wtM_in(maxC_idx+1) = c2n; % next state = min      
        else            
            wtM_ex(1) = c1n; % next state = max
            wtM_in(1) = c2n; % next state = min 
        end 
        wtM = [wtM_ex; wtM_ex; wtM_in];        
    else
        wtM = ones(nmass,1);
    end
    S = wtM.*S;
    S = K.*S;  
       
    % ODEs 
    du = v;
    %dv(1:2*end/3)     = -2/te * v(1:2*end/3)     - 1/te^2 * u(1:2*end/3)     + S(1:2*end/3);
    %dv(2*end/3+1:end) = -2/ti * v(2*end/3+1:end) - 1/ti^2 * u(2*end/3+1:end) + S(2*end/3+1:end);          
    dv = T1.*v + T2.*u + S;          
    
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

function [W wtM] = localConnectivity()
global p
%% The connectivity
    angle_min = 0;
    angle_max = 180;
    d_angle   = 22.5;
    angle     = angle_min:d_angle:angle_max-d_angle;
    
    sel_mat = 255*255;
    
    % Orientation differences
    orientation_difference = [];
    for i = 1:8
        tmp = [];
       for j = 1:8
           tmp = [tmp abs(angle(i)-angle(j))];
       end
       orientation_difference = [orientation_difference; tmp];
    end
    
    % There is only local connection between the masses
    diffAngleLocal = orientation_difference;
    diffAngleLocal(diffAngleLocal>40) = NaN;
    Wpp_local = sel_mat*0.01.*exp(-0.0072*diffAngleLocal); Wpp_local(isnan(Wpp_local)) = 0; 

    % P->eN = local + distal
    Wpe_local = Wpp_local;
    % P->iN = local + distal
    Wpi_local = Wpp_local;
    % Wp = [Wpp Wpe Wpi];
    Wp = [Wpp_local Wpe_local Wpi_local];

    % Excitatory interneuron 
    % eN->P = local
    Wep = Wpp_local;
    % eN->eN = local
    Wee = Wpp_local;
    % eN->iN = local
    Wei = Wpp_local;
    % We = [Wep Wee Wei];
    We = [Wep Wee Wei];

    % Inhibitory interneuron 
    % iN->P = local
    diffAngleLocal = orientation_difference;
    diffAngleLocal(diffAngleLocal>60) = NaN;

    Wip_local = sel_mat*0.08.*exp(-0.038*diffAngleLocal); 
    Wip_local(isnan(Wip_local)) = 0; 
    Wip = Wip_local; 

    % iN-eN = local
    Wie_local = Wip_local;

    % iN-iN = local
    Wii_local = sel_mat*0.04.*exp(-0.038*diffAngleLocal); 
    Wii_local(isnan(Wii_local)) = 0; 

    Wi = [Wip Wie_local Wii_local];   
    
    W = [Wp; We; Wi];
    W = W';
    
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
        
     wtMi = [c1s c1l c1n c1n c1n c1n c1n c1n;
            c1n c1s c1l c1n c1n c1n c1n c1n;
            c1n c1n c1s c1l c1n c1n c1n c1n;
            c1n c1n c1n c1s c1l c1n c1n c1n;
            c1n c1n c1n c1n c1s c1l c1n c1n;
            c1n c1n c1n c1n c1n c1s c1l c1n;
            c1n c1n c1n c1n c1n c1n c1s c1l;
            c1n c1n c1n c1n c1n c1n c1n c1s ];    
        
    wtM = blkdiag(wtMp, wtMp, wtMi);    
    

end

function [W wtM] = connectivity()
    global x y p

    %disp('CONNFUN::Start determining the connectivity...')
    disp('# the connectivity');
    % Sensitive angle
    angle_min = 0;
    angle_max = 180;
    d_angle   = 22.5;
    angle     = angle_min:d_angle:angle_max-d_angle;

    % Sample factor
    sample_factor = 1/3;

    %% Load images and estimate the position of pixel
    %disp('CONNFUN::Load images and estimate the position of pixel');
    map = cell(length(angle),1);
    for ii = 1:length(angle); 
        im = imread(['grafic/opm/' num2str(angle(ii)) '.png']); 
        im = im(:,:,1);
        % cut off the images border
        im1 = double(im(11:103,10:72,1));
        %im1 = im1(1:93,1:63);
        % downsample
        im2 = zeros(size(im1)*sample_factor);
        z = 0;
        for ridx = 1:size(im2,1)      % row indices 
            for cidx = 1:size(im2,2)  % column indices               
                im2(ridx,cidx) = im1(1/sample_factor*(ridx-1)+1,1/sample_factor*(cidx-1)+1); 
            end        
        end
        map{ii} = im2;      
    end


    %% Formula of the connectivity 
    % (see our paper for references)
    % P: Pyramidal
    % eN: excitatory interneuron
    % iN: inhibitory interneuron

    %% No spatial dependence
    % wle = s_x*s_y*0.01*exp(-0.23*abs(dangle));   Long-range: P -> P; P -> eN
    % wli = s_x*s_y*0.03*exp(-0.023*abs(dangle));  Long-range: P -> iN
    % wie = s_x*s_y*0.08*exp(-0.038*abs(dangle));  Short-range: iN -> P; iN -> eN
    % wii = s_x*s_y*0.04*exp(-0.038*abs(dangle));  Short-range: iN -> iN
    % wex = s_x*s_y*0.01*exp(-0.0072*abs(dangle)); Short-range: P -> P;
    % s_x and s_y are sensitivity at (x,y) 
    % dangle = the difference between 2 preffered orientation
    % r = distance between 2 neurons

    %% Blumenfeld said: The result will not depend on if we take the spatial dependences 
    %% into account or not (I also believe it) that's why I only computed the 2nd case
    %% I will test this hypothese for the first case but not now :)

    % Number of column and row
    [nrow ncol] = size(im2);
    n = nrow * ncol; 
    N = 3 * length(angle) * nrow * ncol; % Number of total masses = Number of pyramidal mass + No. exc. interneurons and inh. interneurons

    %% Calculate the selectivity matrix nmass x nmass and angle diffenrence
    %disp(['CONNFUN::Calculate the selectivity matrix ' num2str(N) ' x ' num2str(N) ' and angle diffenrence']);
    sel = [map{1}(:); map{2}(:); map{3}(:); map{4}(:); map{5}(:); map{6}(:); map{7}(:); map{8}(:)];
    len = length(sel);
    sel_mat = sel * sel';

    clear sel

    %% Estimate the matrix of orientation differences
    %disp('CONNFUN::Estimate the matrix of orientation differences');    
    tmp1 = []; for ii = 1:8; tmp1 = [tmp1; repmat(angle(ii),x*y,len)] ; end 
    orientation_difference = abs(tmp1-tmp1');

    % orientation_difference_local = [];  % Orientation difference for local connections
    % orientation_difference_distal = []; % Orientation difference for distal connections
    % 
    % for ii = 1:length(angle)
    %     tmp1 = []; tmp2 = [];
    %     for jj = 1:length(angle)
    %         % local connections (within a hypercolumn)
    %         tmpMat = NaN(n);
    %         tmpMat(logical(eye(n))) = abs(angle(ii)-angle(jj));
    %         tmp1 = [tmp1 tmpMat];
    %         % interhypercolumn connections
    %         tmpMat = orientation_difference((ii-1)*n+1:ii*n,(jj-1)*n+1:jj*n);
    %         tmpMat(logical(eye(n))) = NaN;
    %         tmp2 = [tmp2 tmpMat];
    %     end
    %     orientation_difference_local = [orientation_difference_local; tmp1];
    %     orientation_difference_distal = [orientation_difference_distal; tmp2];
    % end

    % TRK: ALternative
    orientation_difference_local = orientation_difference;
    orientation_difference_distal = orientation_difference;

    for ia = 1:length(angle) 
        for ic = 1:n
            for ja = 1:length(angle) 
                for jc = 1:n
                    if ic == jc
                        orientation_difference_distal((ia-1)*n+ic,(ja-1)*n+jc) = NaN;
                    else
                        orientation_difference_local((ia-1)*n+ic,(ja-1)*n+jc) = NaN;
                    end
                end
            end
        end
    end


    % Constrain the connectivity 
    %orientation_difference_local(orientation_difference_local>40) = NaN;
    %orientation_difference_distal(orientation_difference_distal>60) = NaN;
    clear im* tmp* orientation_difference

    %% Calculate the matrices
    %disp('CONNFUN::Calculate the matrices'); 

    %% Richtig
    % wle = s_x*s_y*0.01*exp(-0.23*abs(dangle));   Long-range: P -> P; P -> eN (>60? -> 0)
    % wli = s_x*s_y*0.03*exp(-0.023*abs(dangle));  Long-range: P -> iN (>60? -> 0)
    % wie = s_x*s_y*0.08*exp(-0.038*abs(dangle));  Short-range: iN -> P; iN -> eN (>60? -> 0)
    % wii = s_x*s_y*0.04*exp(-0.038*abs(dangle));  Short-range: iN -> iN (>60? -> 0)
    % wex = s_x*s_y*0.01*exp(-0.0072*abs(dangle)); Short-range: P -> P; (>40? -> 0)

    %% Pyramidal
    % P->P = local + distal: 
    diffAngleLocal = orientation_difference_local;
    diffAngleDistal = orientation_difference_distal;
    diffAngleLocal(diffAngleLocal>40) = NaN;
    diffAngleDistal(diffAngleDistal>60) = NaN;

    Wpp_local = sel_mat*0.01.*exp(-0.0072*diffAngleLocal); Wpp_local(isnan(Wpp_local)) = 0; 
    Wpp_distal = sel_mat*0.01.*exp(-0.23*diffAngleDistal); Wpp_distal(isnan(Wpp_distal)) = 0; 
    Wpp = Wpp_distal + Wpp_local;

    % P->eN = local + distal
    Wpe = Wpp;

    % P->iN = local + distal
    Wpi_local = Wpp_local;
    Wpi_distal = sel_mat*0.03.*exp(-0.023*diffAngleDistal); Wpi_distal(isnan(Wpi_distal)) = 0; 
    Wpi = Wpi_local + Wpi_distal;

    % Wp = [Wpp Wpe Wpi];
    Wp = [Wpp Wpe Wpi];

    clear Wpp_distal Wpp Wpe Wpi_local Wpi_local

    %% Excitatory interneuron 
    % eN->P = local
    Wep = Wpp_local;
    % eN->eN = local
    Wee = Wpp_local;
    % eN->iN = local
    Wei = Wpp_local;
    % We = [Wep Wee Wei];
    We = [Wep Wee Wei];
    clear Wep Wee Wei

    %% Inhibitory interneuron 
    % iN->P = local
    diffAngleLocal = orientation_difference_local;
    diffAngleLocal(diffAngleLocal>60) = NaN;

    Wip_local = sel_mat*0.08.*exp(-0.038*diffAngleLocal); 
    Wip_local(isnan(Wip_local)) = 0; 
    Wip = Wip_local; 

    % iN-eN = local
    Wie_local = Wip_local;

    % iN-iN = local
    Wii_local = sel_mat*0.04.*exp(-0.038*diffAngleLocal); 
    Wii_local(isnan(Wii_local)) = 0; 

    Wi = [Wip Wie_local Wii_local];

    clear Wie_local Wii_local Wip Wpi Wpp Wpe orientation_difference* sel_mat

    W = [Wp; We; Wi]; clear Wp We Wi

    W = W';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%% Weighting Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('# the weighting matrix');  
%     c1s = p(10); c1n = p(11); c1l = p(12);
%     c2s = p(13); c2n = p(14); c2l = p(15);   
%     % Pyramidal  
%     tmp1 = [repmat(c1s,x*y,x*y) repmat(c1l,x*y,x*y) repmat(c1l,x*y,5*x*y) repmat(c1n,x*y,x*y)];
%     wtMp = tmp1;
%     tmp2 = [repmat(c2s,x*y,x*y) repmat(c2l,x*y,x*y) repmat(c2l,x*y,5*x*y) repmat(c2n,x*y,x*y)];
%     wtMi = tmp2;
%     for i = 1:7
%         wtMp = [wtMp; circshift(tmp1,[0 i*x*y])];
%         wtMi = [wtMi; circshift(tmp2,[0 i*x*y])];
%     end
%    
%     % Put all together
%     %tmp = zeros(size(wtMp));
%     %wtM = [ wtMp              tmp               tmp; 
%     %        tmp               wtMp              tmp; 
%     %        tmp               tmp               wtMi  ];
%     wtM = blkdiag(wtMp,wtMp,wtMi);
%     clear wtMp wtMi tmp*
%         
    %wtM = wtM';
end




