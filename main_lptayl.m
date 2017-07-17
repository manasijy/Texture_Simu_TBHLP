function [gmatrix0] = main_lptayl(filename,incremental_strain,n_iteration,varargin)

%% 
% Syntax: main_lptayl(filename,incremental_strain,n_iteration,varargin)
% main_lptayl_function('myfile.txt',0.05,20,'criteria',0,...
% 'ConstraintsRelaxed',1,'slip_type','fp')
% inputs:
% filename: it must be a text file e.g. 'myfile.txt'
% incremental_strain: the strain applied in each iteration- 0.01,0.5 etc
% n_iteration: no of iteration to be run
% varagin:
%     slip_type:  f (full), p (only partial), fp (full+partial), ft(full+twin)
%     criteria: criteria to solve taylor ambiguity
%     0 (minVar), 1(maxVar), 2(minPlasticSpin), 3(wintenberger):2,3 to do
%     cryst_str:  f(fcc) and b(bcc) :To do
%     ConstraintsRelaxed: 0 (Full)\t 1 (Lath)\t 2 (PanCake): To do

%%

cs = crystalSymmetry('-43m'); 
ss = specimenSymmetry('mmm');
g_file = filename; 
namesplit = strsplit(filename,'.');
name = namesplit(1); % namesplit(2) is extension
gmatrix0 = txtfile2ori(g_file,cs,ss);
plotTexture(gmatrix0,[1]);%(gmatrix0,cs,ss,[1]);
if isempty(gmatrix0), return, end

%% Get other inputs or set defaults

if check_option(varargin,'slip_type')
  slip_type = get_option(varargin,'slip_type');
else
  slip_type = 'f';
end

if check_option(varargin,'criteria')
  criteria = get_option(varargin,'criteria');
else
  criteria = 0;
end

if check_option(varargin,'ConstraintsRelaxed')
  ConstraintsRelaxed = get_option(varargin,'ConstraintsRelaxed');
else
  ConstraintsRelaxed = 0;
end

if check_option(varargin,'cryst_str')
  cryst_str = get_option(varargin,'cryst_str');
else
  cryst_str = 0;
end

tic

%% Initializing the cost matrix according to type of slip mode
c = ones(1,48);
% full slip: n = 1-12,25-36 % partial slip: n=13-24,37-48 % twin:n=13-24

switch slip_type
    case {'F', 'f'} % only perfect slip
        n = [13:24,37:48];
        c(n)=1e4; % All partial modes are made costlier        
    case {'P','p'} % only partial slip
        n = [1:12,25:36];
        c(n)=1e4; % All perfect modes are made costlier
    case {'FP','fp','Fp','fP'} %both perfect and partial slips
        
    case {'FT','ft','Ft','fT'}
        n = 37:48;
        c(n)=1e4; % All - partial modes are made costlier        
    otherwise
        warning('Unexpected slip type.\n {Default perfect slip is applied');        
 end
        


%% Initialization

% Criteria will be given as input in later development. But for present
% minvar is taken as the single unique solution selection criteria. Others
% are 'maxvar' and 'minplasticspin'
% criteria = 'minvar';
% criteria = 'maxvar';
% criteria = 'minplasticspin';

SlipSystem = SS_setFCC_fp_function(); % need to cover for other crystal systems
A = getA(SlipSystem);
gmatrix = gmatrix0;

lg =  length(gmatrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
new_gmatrix(:,1) = gmatrix;

%% Calculation Block

for j=1:1:n_iteration

for i=1:1:lg 
    DCMtrx = matrix(gmatrix(i));  
    e_grain=DCMtrx'*e_ext*DCMtrx; 
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    AllSolutions = calcLPSolution(b,c,A); 
    [spin, UniqueSolution] = calculate_spin(AllSolutions,criteria,SlipSystem);
    spin = incremental_strain*spin;
    gmatrix(i) = applySpin(-spin,DCMtrx); 
    new_gmatrix(i,j+1) = gmatrix(i);
    ActiveSS(i,j).ss= UniqueSolution.B;
    ActiveSS(i,j).gamma = UniqueSolution.xb;
end
end
toc
%% Output Block
prompt = 'Do you want to save results; yes or no \n';
reply = input(prompt,'s');
if strcmp(reply,'yes')
    disp('Active slip systems with respective shears \nwill be stored in ActiveSSData matrix in the work space\n')
    disp('Texture data at each iteration will be stored\n in the folder filename in the current path\n')
    folderpath = resmet_input(name,new_gmatrix);
    ActiveSSData = [folderpath '\' 'ActiveSSData.mat'];
    save(ActiveSSData,'ActiveSS');    
    %%%Following lines to save o/p euler angles in file
    %     eu_gmat = Euler(new_gmatrix)/degree;
    %     euler_angles(name,new_gmatrix);
% else
%     return
end
%%
prompt = 'Do you want to plot results; yes or no \n';
reply = input(prompt,'s');
if strcmp(reply,'yes')
   disp('Which iterations you want to plot');
   whichOnes = input('For only last iteration: type l, for all enter 0\n else enter iteration numbers e.g. 1,2,3\n');
   switch whichOnes
       case 'l'
           iteration_list = n_iteration+1;
       case '0'
           iteration_list = 1:(n_iteration+1);
       otherwise
           iteration_list = sscanf(whichOnes,strcat("%d",','))';
           iteration_list = [1,iteration_list+1];
   end  
   plotTexture(new_gmatrix,iteration_list) 
end
