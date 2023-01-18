function [qA, thresh, oFile] = getHydroVolcSweepSummary(sweepFile,MCflag,saveOut)
% [qA, thresh, oFile] = getHydroVolcSweepSummary(sweepFile,MC,saveOut)
% Retrieves summary data for a set of hydroVolc runs. Note: this script
% assumes some parameters do not change between runs inside the same file
%   sweepFile = either path to sweep output file, or the struct itself
%   MC        : true = Monte Carlo run, false = sweep. Only needed if
%               sweepFile input is a struct (not a path)
%     saveOut = [true/false] - will save output file to the same directory
%                           as input
%             = 'char'  - saveOut is used as the output save directory


if nargin<3
    saveOut = false;
end
if nargin<2
    MCflag = [];
end
if isempty(MCflag)
    MCflag = false;
end

if ischar(saveOut)
    oDir    = saveOut;
    saveOut = true;
end

oFileName = 'outputSummary_';

%%
if ischar(sweepFile)
    T = load(sweepFile);
    
    % Fiddling with model versioning/var names a bit... var name could also
    % be an input?
    if isfield(T,'qSet')
        dat = T.qSet;
    elseif isfield(T,'dat')
        dat = T.dat;
    end
    MCflag = isfield(T,'randPars'); % Overwrites user input

    clear T
    
elseif isstruct(sweepFile)
    dat = sweepFile;
end

% if isfield(dat(1).cI,'proxy') % Currently assuming all sweep runs are proxied or not. Eventually need to generalize the array function to accommodate either
%     if dat(1).cI.proxy
%         useProxy = true;
%     else
%         useProxy = false;
%     end
% else
%     useProxy = false;
% end
    
if isrow(dat)
    dat = dat';
end
%%
 % Currently assuming all sweep runs are proxied or not. Eventually need to generalize the array function to accommodate either
[qA, thresh] = getCoupledSweepArrays(dat,MCflag);

if saveOut
    [oDir,oName,oext] = fileparts(sweepFile);
    oFile = fullfile(oDir,[oFileName oName oext]);
    fprintf('Saving output summary:\n\t%s\n',oFile)
%     save(oFile,'qA','M_thresh','K_thresh','LdZ_thresh','Ld0_thresh','fPlume_thresh','LdLj_thresh','Clps_lo','Clps_hi')
    save(oFile,'qA','thresh')
else
    oFile = '';
end
    

% if useProxy
%     qA = getProxySweepArrays(dat);
%     thresh = [];
%     
%     % Proxy conduit
%     if saveOut
%         [oDir,oName,oext] = fileparts(sweepFile);
%         oFile = fullfile(oDir,[oFileName oName oext]);
%         fprintf('Saving output summary:\n\t%s\n',oFile)
%     %     save(oFile,'qA','M_thresh','K_thresh','LdZ_thresh','Ld0_thresh','fPlume_thresh','LdLj_thresh','Clps_lo','Clps_hi')
%         save(oFile,'qA')
%     else
%         oFile = '';
%     end
%     
%     
% else
%     [qA, thresh] = getCoupledSweepArrays(dat);
%     
%     % Full
%     if saveOut
%         if ~exist('oDir','var')
%             [oDir,oName,oext] = fileparts(sweepFile);
%         else
%             [~,oName,oext] = fileparts(sweepFile);
%         end
%         oFile = fullfile(oDir,[oFileName oName oext]);
%         fprintf('Saving output summary:\n\t%s\n',oFile)
%         save(oFile,'qA','thresh')
%     else
%         oFile = '';
%     end
% 
% end

%%



end