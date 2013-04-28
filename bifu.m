function bifu(p_filename,res_filename)
% Version infomation:
% 1 - Weighting the matrix in the connectivity function
% 2 - Wichten the selectivity smax unterschiedlich fuer excitatory und inhibitory
% 3 - 09.05.13 Change structure of p 

    %% close and clear all
    close all; 
    warning off;  
    addpath('input');
         
    %% Define global variables
    global  x y nmass te ti he hi e0 r0 v0 n sz neqn tend ...
        c1s c1n c1l c2s c2n c2l f1 f2 f3 f4 f5 p sample_factor path 
     
    %% Load initial data file  
    % te=p(1); ti=p(2); he=p(3); hi=p(4); 
    % sz = p(5); sample_factor = p(6)
    % c1s=p(7); c1n=p(8); c1l=p(9);
    % c2s=p(10); c2n=p(11); c2l=p(12);
    % f1 = p(13); f2 = p(14); f3 = p(15);
    % f4 = p(16); f5 = p(17); tend = p(18);
    % disp('Loading initial parameter...');
    %load(['param/' filename '.mat']);
    load([p_filename '.mat']);
    fprintf('[%.3f %.3f %.3f %.3f %.3f %.3f] %s %s %.4f %.4f %.4f %.4f %.4f\n',...
        p(7),p(8),p(9),p(10),p(11),p(12),p_filename,res_filename,p(13),p(14),p(15),p(16),p(17));     
    
    %% Parameters
    % synaptic and dendrite parameter
    te      = p(1);% 10  % ms
    ti      = p(2);% 20; % ms
    he      = p(3);% 3.25;  % mV
    hi      = p(4);% -22;% mV
    sz      = p(5);
    c1s     = p(7);
    c1n     = p(8);
    c1l     = p(9);
    c2s     = p(10);
    c2n     = p(11);
    c2l     = p(12);
    f1      = p(13);
    f2      = p(14);
    f3      = p(15);
    f4      = p(16);
    f5      = p(17);
    tend    = p(18);
    sample_factor = p(6);
    % parameter of sigmoid function
    e0       = 2.5e-3;% 1/ms
    r0       = 0.56;  % 1/mV
    v0       = 6;     % mV
    % others
    x        = sample_factor*90;% 31 field size
    y        = sample_factor*60;% 21 field size
    n        = x*y;
    nmass    = n*3*8;
    neqn     = 2*nmass;
       
    % Result path
    path     = 'output';

    %% Solving the ode system
    [OM] = VC_ode();
    fname = [path '/c1s=' num2str(p(10)) '_c1n=' num2str(p(11)) '_c1l=' num2str(p(12))...
        '_c2s=' num2str(p(13)) '_c2n=' num2str(p(14)) '_c2l=' num2str(p(15)) '_f1='...
        num2str(f1) '_f2=' num2str(f2) '_f3=' num2str(f3) '_f4=' num2str(f4) '_f5=' num2str(f5)];
    %save([fname '.mat'],'OM');
    plotLenin(fname,OM);  
    
    OM = OM(:,end-200:end);
    OM = [ sum(OM(1:n,:));       sum(OM(n+1:2*n,:)); sum(OM(2*n+1:3*n,:)); ...
           sum(OM(3*n+1:4*n,:)); sum(OM(4*n+1:5*n,:)); sum(OM(5*n+1:6*n,:)); ...
           sum(OM(6*n+1:7*n,:)); sum(OM(7*n+1:8*n,:)) ];

    minOM = min(OM,[],2);
    maxOM = max(OM,[],2);
    res = [f1 f2 f3 f4 f5 minOM(:)' maxOM(:)'];  
    
    fid = fopen([res_filename '.txt'],'a+');    
    for i = 1:length(res)
        fprintf(fid,'%.6f\t',res(i)); 
    end
    fprintf(fid,'\n');
    fclose(fid);
    
    for i = 6:8; fprintf('%.9f\t',res(i)); end;  fprintf('\n');
    for i = 14:16; fprintf('%.9f\t',res(i)); end; fprintf('\n');
  
    clear, clearvars -global
    %disp('Finished!');
end

function [OM] = VC_ode()
%% ODE solver using Euler method with constant stepsize 
global te ti he hi e0 r0 v0 nmass sz tend

% The connectivity and the weighting matrix
%disp('Loading the connectivity matrix...');
[W] = connectivity(); 
    
% modify the connectivity
ke = he/te; ki = hi/ti;
K  = [repmat(ke,2*nmass/3,1); repmat(ki,nmass/3,1)];
T1 = [repmat(-2/te,2*nmass/3,1); repmat(-2/ti,nmass/3,1)];
T2 = [repmat(-1/te^2,2*nmass/3,1); repmat(-1/ti^2,nmass/3,1)];

trunk = 200;
n_trunk = ceil(tend/trunk/sz);

%%
%% Solver main loop
%%
% wtM = ones(nmass,1);
OM = [];  % Orientation map
%Y0 = [zeros(nmass,1); zeros(nmass,1)];% Initial condition
%Y0(1) = -70; Y0(nmass+1) = -70;
Y0 = rand(2*nmass,1);
dv = zeros(nmass,1);
du = dv;

%fprintf('Solver main loop...\n');

% additional info
trunkidx = 1; 

fprintf('Processing (total %d):    \t',n_trunk);

for idx = 1:length(0:sz:tend)      
    % Update half of Y0   
    u = Y0(1:nmass); 
    v = Y0(nmass+1:2*nmass);
   
    % sigmoid function of membrane potentials    
    S = 2*e0 ./ (1+exp(r0*(v0-u))); 
    
    % Weight S by the quotient h/tau  
    S = W*S;  
       
    % ODEs 
    du = v;
    %dv(1:2*end/3)     = -2/te * v(1:2*end/3)     - 1/te^2 * u(1:2*end/3)     + S(1:2*end/3);
    %dv(2*end/3+1:end) = -2/ti * v(2*end/3+1:end) - 1/ti^2 * u(2*end/3+1:end) + S(2*end/3+1:end);          
    dv = T1.*v + T2.*u + K.*S;          
    
    % no noise
    % dY = [du; dv]
    Y1 = Y0 + sz * [du; dv];   
 
    % Update Y
    Y0 = Y1;
       
    % Save result in each step   
    OM = [OM u(1:end/3)+u(end/3+1:2*end/3)+u(2*end/3+1:end)];
      
    if ~mod(idx,trunk)
       % display the informations
       fprintf('\b\b\b\b%04d',trunkidx);
       trunkidx = trunkidx + 1;
    end

end

%% End of solver loops
fprintf('\n');    
end


function [W] = connectivity()
% Changed 23.04.13: Die Winkeldifferenz zwischen 157.5 und 0 is 22.5 
% (nicht 157.5 wie vorher)

    global x y f1 f2 f3 f4 f5 c1s c1n c1l c2s c2n c2l sample_factor 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%% Weighting Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp('# the weighting matrix');   
    % Pyramidal and excitatory  
    tmp1 = [repmat(c1s,x*y,x*y) repmat(c1l,x*y,6*x*y) repmat(c1n,x*y,x*y)];
    wtMp = tmp1;
    % Inhibitory  
    tmp2 = [repmat(c2s,x*y,x*y) repmat(c2l,x*y,6*x*y) repmat(c2n,x*y,x*y)];
    wtMi = tmp2;
    for i = 1:7
        wtMp = [wtMp; circshift(tmp1,[0 i*x*y])];
        wtMi = [wtMi; circshift(tmp2,[0 i*x*y])];
    end
    clear tmp1 tmp2
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%% Connectivity Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp('# the connectivity');
    % Sensitive angle
    angle_min = 0;
    angle_max = 180;
    d_angle   = 22.5;
    angle     = angle_min:d_angle:angle_max-d_angle;

    %% Load images and estimate the position of pixel
    %disp('CONNFUN::Load images and estimate the position of pixel');
    map = cell(length(angle),1);
    for ii = 1:length(angle); 
        im = imread(['grafic/opm/' num2str(angle(ii)) '.png']); 
        im = im(:,:,1);
        % cut off the images border
        im1 = double(im(11:103,10:72,1));
        im1 = im1(2:91,2:61);
        % downsample
        im2 = zeros(size(im1)*sample_factor);
        for ridx = 1:size(im2,1)      % row indices 
            for cidx = 1:size(im2,2)  % column indices               
                %im2(ridx,cidx) = im1(1/sample_factor*(ridx-1)+1,1/sample_factor*(cidx-1)+1); 
                im2(ridx,cidx) = sum(sum(im1(1/sample_factor*(ridx-1)+1:1/sample_factor*ridx,...
                    1/sample_factor*(cidx-1)+1:1/sample_factor*cidx)))/(1/sample_factor*1/sample_factor); 
            end        
        end
        map{ii} = im2;      
    end

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
    orientation_difference(orientation_difference==157.5) = 22.5;
    orientation_difference(orientation_difference==135)   = 45;
    orientation_difference(orientation_difference==112.5) = 67.5;

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

    %% Formula
    % wle = s_x*s_y*0.01*exp(-0.23*abs(dangle));   Long-range: P -> P; P -> eN (>60? -> 0)
    % wli = s_x*s_y*0.03*exp(-0.023*abs(dangle));  Long-range: P -> iN (>60? -> 0)
    % wie = s_x*s_y*0.08*exp(-0.038*abs(dangle));  Short-range: iN -> P; iN -> eN (>60? -> 0)
    % wii = s_x*s_y*0.04*exp(-0.038*abs(dangle));  Short-range: iN -> iN (>60? -> 0)
    % wex = s_x*s_y*0.01*exp(-0.0072*abs(dangle)); Short-range: P -> P; (>40? -> 0)
    
    %% 07.05.2013
    % use 5 weighting factor for 5 types of connectivity
    % wle = f1*s_x*s_y*0.01*exp(-0.23*abs(dangle));   Long-range: P -> P; P -> eN (>60? -> 0)
    % wli = f2*s_x*s_y*0.03*exp(-0.023*abs(dangle));  Long-range: P -> iN (>60? -> 0)
    % wie = f3*s_x*s_y*0.08*exp(-0.038*abs(dangle));  Short-range: iN -> P; iN -> eN (>60? -> 0)
    % wii = f4*s_x*s_y*0.04*exp(-0.038*abs(dangle));  Short-range: iN -> iN (>60? -> 0)
    % wex = f5*s_x*s_y*0.01*exp(-0.0072*abs(dangle)); Short-range: P -> P; (>40? -> 0)
    
    %% Pyramidal
    % P->P = local + distal: 
    diffAngleLocal = orientation_difference_local;
    diffAngleDistal = orientation_difference_distal;
    diffAngleLocal(diffAngleLocal>40) = NaN;
    diffAngleDistal(diffAngleDistal>60) = NaN;

    Wpp_local = 1/f5*sel_mat*0.01.*exp(-0.0072*diffAngleLocal); Wpp_local(isnan(Wpp_local)) = 0; 
    Wpp_distal = 1/f1*sel_mat*0.01.*exp(-0.23*diffAngleDistal); Wpp_distal(isnan(Wpp_distal)) = 0; 
    Wpp = Wpp_distal + Wpp_local;

    % P->eN = local + distal
    Wpe = Wpp;

    % P->iN = local + distal
    Wpi_local = Wpp_local;
    Wpi_distal = 1/f2*sel_mat*0.03.*exp(-0.023*diffAngleDistal); Wpi_distal(isnan(Wpi_distal)) = 0; 
    Wpi = Wpi_local + Wpi_distal;

    % Wp = [Wpp Wpe Wpi];
    Wp = [wtMp.*Wpp wtMp.*Wpe wtMi.*Wpi];
    clear Wpp_distal Wpp Wpe Wpi_local Wpi_local

    %% Excitatory interneuron 
    % eN->P = local
    Wep = Wpp_local;
    % eN->eN = local
    Wee = Wpp_local;
    % eN->iN = local
    Wei = Wpp_local;
    % We = [Wep Wee Wei];
    We = [wtMp.*Wep wtMp.*Wee wtMi.*Wei];
    clear Wep Wee Wei

    %% Inhibitory interneuron 
    % iN->P = local
    diffAngleLocal = orientation_difference_local;
    diffAngleLocal(diffAngleLocal>60) = NaN;

    Wip_local = 1/f3*sel_mat*0.08.*exp(-0.038*diffAngleLocal); 
    Wip_local(isnan(Wip_local)) = 0; 
    Wip = Wip_local; 

    % iN-eN = local
    Wie_local = Wip_local;

    % iN-iN = local
    Wii_local = 1/f4*sel_mat*0.04.*exp(-0.038*diffAngleLocal); 
    Wii_local(isnan(Wii_local)) = 0; 

    % Weight Wi by the asymmetric coefficients
    Wi = [wtMi.*Wip wtMi.*Wie_local wtMi.*Wii_local];

    clear Wie_local Wii_local Wip Wpi Wpp Wpe orientation_difference* sel_mat

    W = [Wp; We; Wi]; clear Wp We Wi

    W = W';

end



