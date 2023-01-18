% MER sweep sandboxing, conduit V6
%  ---> Shooting runs to find the range of viable conduit solutions for
%  "dry" runs, aka Hw=0
%       --> Conduit stop conditions: Pm~Pf, U~0, M>~.999
%       --> Require Z(end)~0 AND either (a) choked or (b) pm(end) = pf
%           --> So watch out for runs where stop condition is M>=.999
%  ---> Test a range of vent pressures (~water depths) on the base case to
%   assess response, adjusting MER or radius (phi0 or Z0 as options?)
%
% PROCEDURE (should write as accessible functions to refine full look-up table):
%     (1) "runCoarse": Coarse brute force sweep with a coarse step size to find
%           approximate bounding locations for conduit radii.
%           --> Generate a coarse lookup table in R,MER,Pf
%           --> CURRENTLY SAVES ALL CONDUIT OUTPUT, FILE CAN BE VERY LARGE
%
%     (2) "getParams": Get summary values from coarse search, sufficient to
%           run fine search.
%           --> Could combine with (1) if you don't care to save coarse
%           output
%
%     (3) "runRefine": Use initial coarse table to:
%           (a) refine precision bounds on R for each MER,Pf 
%               --> (or go the other way and refine precision MER for given R,Pf)?
%
%   NOT YET FULLY IMPLEMENTED
%           (b) interpolate guesses for arbitrary R, MER, Pf from fine
%           search lookup table.
%           (c) run a fine search down to a tolerance from interpolations.
%               '-> This should be a fast-ish function to do a final,
%                   repeatable shoot (e.g. get MER for given input R, Pf)
%
%
% C Rowell, May 2021
% Uses Hajimirza conduit model, V6


clear all; close all
addpath(genpath('C:\Users\crowell\Documents\GitHub\glaciovolc\glaciovolc-dev\'))


%% ============================== INPUT ===================================


% ---- COARSE SEARCH -----
% Test params for guessing R
% T0 = 850+273.15;     %Gas Temp
% rho_melt = 2400;
% n0       = .035;
rho      = 110;
u        = 140;

% Also used for fine search
ZfailScale = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .04;  % Over-/underpressure threshold


% Test search
% logMER = [7.5 8]; %6:0.1:9;
% rscale = [0.8:0.1:1.5]; %(0.6:.05:3.5);
% pf     = [1e5 1e6];
% Full search
logMER = (5:.1:9);  %6:0.1:9;
rscale = (0.6:0.1:3.6); %(0.6:.05:3.5);
pf     = (1e5:1e5:5.1e6);

% ---- FINE SEARCH -----
% rscaleF = [0.9:0.01:1.1];
% Rmin = 0.05;
dRminScale = 0.2e-3;
maxIter = 10;

% ---- INTERPOLATION -----
% Spacing in R


% Data directory
oDir        = '/Users/crrowell/Kahuna/data/glaciovolc/conduitSweeps/';
% oDir        = 'C:\Users\crowell\Kahuna\data\glaciovolc\conduitSweeps\';


% ------ DIRECTORIES AND SWITCHES  -------
% Generally best to run these one at a time.
runCoarse  = false;
    coarseName = 'conduitV6_coarseSweep_n%i_%s'; % Name template for sweep file
    % Coarse sweep output file name:
    coarseFile = fullfile(oDir,'BigFatCoarseSweep_2021-05\conduitV6_coarseSweep_n64821_21-05-16.mat');

% (setting any to false attempts to load existing files instead)
getParams  = false;
    % Reduced file:
    compFile = fullfile(oDir,'conduitV6_coarseSweep_n64821_21-05-16_compressedV4.mat');
    compFile = fullfile(oDir,'conduitV6_coarseSweep_n64821_21-05-16_compressedV3.mat');
    % compFile = fullfile(oDir,'conduitV6_coarseSweep_n32_21-05-10_compressed.mat');


runRefine  = false;
    verbose = 1; % 0 = none (progress bar), 1 = single line each run, 2 = full search output
    % Fine search output file:
    fineFile = fullfile(oDir,'conduitV6_fineSweep_n64821_21-05-16_compressedV4_fineV2.mat');

runInterp  = true;

runPlots   = false;


% =========================================================================
%% (1) RUN COARSE FULL SWEEP
MER     = 10.^logMER;
rvGuess = (MER./(u*rho*pi)).^(1/2);

if runCoarse
    nr = length(MER)*length(rscale)*length(pf);
    fprintf('Coarse search, %i total runs:\n',nr)

    % Default start conditions
    % cI0.pf = pf(1);
    cI0 = getConduitSource;

    % Coarse run setup
    cIc(length(rscale),length(MER),length(pf)) = table2struct(cI0);
    Zpass = false(size(cIc)); % Top height falls below ~0 (at MER>=.999, or U=0)
    Upass = false(size(cIc)); % Not choked at top height
%     Ppass = false(size(cIc)); % Not pressure balanced at vent (only needed if Upass fails)

    ct = 0;
    textprogressbar('')
    for kk=1:length(pf)
        for jj=1:length(MER)
            for ii=1:length(rscale)
                ct = ct+1;
                textprogressbar(ct/nr*100)

                C                = cI0;
                C.Q              = MER(jj);
                C.conduit_radius = rvGuess(jj)*rscale(ii);
                C.pf             = pf(kk);

                cOc(ii,jj,kk) = Conduit_flow_with_nucleation_V6(C);
                cIc(ii,jj,kk) = table2struct(C);

                Zpass(ii,jj,kk) = cOc(ii,jj,kk).Z(end)<zFailThresh;
                Upass(ii,jj,kk) = cOc(ii,jj,kk).M(end)>Mfailthresh;
            end
        end
    end
    textprogressbar(' --> Done')
    oFile = fullfile(oDir,sprintf(coarseName,nr,datestr(now,'YY-mm-dd')));
    save(oFile,'cIc','cOc','MER','rvGuess','rscale','pf','Zpass','Upass')

end

%% Slimmed down output files for easier porting/loading/processing
if getParams
    load(coarseFile) % use with caution
    
    nr = numel(cIc);
    Mc  = zeros(size(cIc)); % Mach #
    Ztc = Mc;               % top depth
    pmc = Mc;               % Final pressure
    ric = Mc;               % Initial radius
    pfc = Mc;               % Surface pressure
    avc = Mc;               % Final vent radius
    uc  = Mc;               % Velocity
    rhoc= Mc;               % Bulk density

    textprogressbar('')
    for ii=1:nr
        textprogressbar(ii/nr*100)
        Mc(ii) = cOc(ii).M(end);
        Ztc(ii) = cOc(ii).Z(end);
        pmc(ii) = cOc(ii).pm(end);
        pfc(ii) = cIc(ii).pf;
        ric(ii) = cIc(ii).conduit_radius;
        avc(ii) = cOc(ii).a(end);
        uc(ii)  = cOc(ii).U(end);
        rhoc(ii)= cOc(ii).rho_magma(end);
    end
    Ppass = pmc<(pfc.*Pfailthresh);
    textprogressbar(' --> Done')
    
    disp('Saving...')
    [~,oname,~] = fileparts(coarseFile);
    oname = fullfile(oDir,[oname '_compressedV4']);
    fprintf('Saving:\n\t%s\n',oname)
    save(oname,'cIc','MER','rvGuess','rscale','pf','Zpass','Upass','Mc','Ztc','uc','rhoc','pmc','pfc','ric','avc','Ppass')
    clear cOc
end

%% Fine search
load(compFile)
if runRefine

    % Temp selection
    qi = find(log10(MER)==8);
    
%     nr = numel(cIc(:,qi,:));
    nSearches = prod(size(cIc,[2 3]));
    fprintf('Fine search, %i total runs:\n',nSearches )

    % Selection params for initial bounds
        % -> Z<Zthresh % Find highest R fail
        % -> Mach #: for Z>Zthresh, find first value below Uthresh
    
    dd = zeros(size(Mc,2),size(Mc,3));
    cS = struct('Rmax',dd,'Rmin',dd,'Ztop',dd,'M',dd,'pm',dd,'rho_b',dd,...
        'U',dd,'amax',dd,'dPdt_fr',dd,'phi_g',dd,'Z_fr',dd,'U_fr',dd,...
        'pfr',dd,'success',false(size(dd)),'rSrchBounds',repmat(dd,[1 1 2]),...
        'zFailThresh',dd);
    cI(size(cIc,2),size(cIc,3)) = cIc(1);
    
    % ----- Temp fix from updated name hw-->Zw
    cI(size(cIc,2),size(cIc,3)).Zw = cIc(1).hw;
    cI = rmfield(cI,'hw');
    % -------------
    
    cIf = getConduitSource;
    ct = 0;
    if verbose==0
        textprogressbar('')
        textprogressbar(0)
    end
    if verbose ==2
        full_verbose = true;
    else
        full_verbose = false;
    end

    for kk=1:length(pf)
        for jj=1:length(MER)
%             for ii=1:length(rscale)
            ct = ct+1;
            tic
            cIf.pf  = pf(kk);
            cIf.Q   = MER(jj);
            Rconflict = false;
            
            % Tests to bracket the success range
            % === NEW TESTS ===
            zPass = Ztc(:,jj,kk)<(ric(:,jj,kk)*ZfailScale);
            upPass= ( (pmc(:,jj,kk) + uc(:,jj,kk).^2.*rhoc(:,jj,kk)/2) ./ pfc(:,jj,kk) ) > (1-Pfailthresh);
            pPass = and((pmc(:,jj,kk) ./ pfc(:,jj,kk)) < (1+Pfailthresh) , upPass);
            uPass = and( Mc(:,jj,kk)>Mfailthresh, Mc(:,jj,kk) < 1.1);
            fChk  = abs((avc(:,jj,kk)-ric(:,jj,kk))./ric(:,jj,kk))>.005; % .5% difference threshold for now
            pass  = and(zPass,upPass) & or(pPass,uPass);
            
            % =================

            % Check that top is reached, not underpressured
            PUnder = ~(zPass & upPass);

            % Get search bounds and any valid results
            if any(pass)
                iValid = [find(pass,1,'first') find(pass,1,'last')];
                Rvalid = unique(ric([find(pass,1,'first') find(pass,1,'last')],jj,kk));
%                 target2 = find(pass,1,'first');
                rIdx   = [min(iValid)-1 max(iValid)+1];
                
                % Check bounds
                if rIdx(1)<1
                    Rbounds(1) = ric(1,jj,kk)-diff(ric(1:2,jj,kk));
                else
                    Rbounds(1) = ric(rIdx(1),jj,kk);
                end                
                if rIdx(2)>size(ric,1)
                    Rbounds(2) = ric(end,jj,kk)+diff(ric(end-1:end,jj,kk));
                else
                    Rbounds(2) = ric(rIdx(2),jj,kk);
                end
%                 Rbounds = ric(([min(iValid)-1 max(iValid)+1]),jj,kk);
            else
                Rvalid = [];
                zf1   = find(zPass,1,'first'); % First Z success
%                 iMin  = max([zf1-1 find(PUnder,1,'last')]);
                iMin  = nanmax([zf1-1 find(~upPass,1,'last')]);
                iMax  = find((zPass & ~uPass & ~fChk),1,'first');
                if isempty(iMax)
                    iMax  = find(((1:length(zPass))'>=zf1 & ~uPass & ~fChk),1,'first');
                end
                rIdx    = [iMin iMax];
                Rbounds = ric(rIdx,jj,kk);
                
                % Check iMin<iMax
                while (Rbounds(2)-Rbounds(1))<(mean(diff(ric(:,jj,kk)))*.95)
                    if (Rbounds(2)<Rbounds(1)); Rconflict = true; end
                    Rbounds(1) = Rbounds(1)-mean(diff(ric(:,jj,kk)))/2;
                    Rbounds(2) = Rbounds(2)+mean(diff(ric(:,jj,kk)))/2;
                end
                
            end
            
            % Dynamic Zfailthresh to within ~1 conduit radii
            cS.zFailThresh(jj,kk) = round(mean(Rbounds))*ZfailScale;
            [Rrange,cIi,cO,success] = conduitRadiusFromQ(cIf,Rbounds,'Rvalid',Rvalid,'maxIter',maxIter,...
                'Zfailthresh',cS.zFailThresh(jj,kk),'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
                'dRminScale',dRminScale,'verbose',full_verbose,'output',true);
            cI(jj,kk)       = table2struct(cIi);
            cS.success(jj,kk) = success;
            cS.pf(jj,kk)    = pf(kk);
            cS.Q(jj,kk)     = MER(jj);
            cS.Rmax(jj,kk)  = max(Rrange);
            cS.Rmin(jj,kk)  = min(Rrange);
            cS.Ztop(jj,kk)  = cO.Z(end);
            cS.M(jj,kk)     = cO.M(end);
            cS.pm(jj,kk)    = cO.pm(end);
            cS.rho_b(jj,kk) = cO.rho_magma(end);
            cS.U(jj,kk)     = cO.U(end);
            cS.amax(jj,kk)  = cO.a(end);
            cS.phi_g(jj,kk) = cO.porosity(end);
            cS.Z_fr(jj,kk)  = cO.Par.Zf;          
            if cO.Par.frag
                fi              = find(cO.Z==cO.Par.Zf,1,'first');
                cS.dPdt_fr(jj,kk)= cO.dPdt(fi); % Peak dPdt at fragmentation
                cS.pfr(jj,kk)   = cO.pm(fi);
                cS.U_fr(jj,kk)   = cO.U(fi);     % Velocity at fragmentation
            else
                cS.dPdt_fr(jj,kk)= max(cO.dPdt); % Peak dPdt at fragmentation
                cS.pfr(jj,kk)   = NaN; %cO.pm(fi);
                cS.U_fr(jj,kk)   = NaN; %cO.U(fi);     % Velocity at fragmentation
            end
            
            if verbose==1
                fprintf(' %i/%i ... Q: %.2e, Pf: %.2f, Rb_lo: %.3f\tRb_hi: %.3f\tR_confl: %i\tRf_lo: %.4f\tRf_hi: %.4f\tV?: %i\n',...
                    ct,nSearches,MER(jj),pf(kk)/1e6,Rbounds(1),Rbounds(2),Rconflict,min(Rrange),max(Rrange),success)
            elseif ~verbose
                textprogressbar(ct/nSearches*100)
            end
        end
    end
    
    if ~verbose
        textprogressbar('-> done')
    else
        disp('-> done')
    end
    
    MER = MER';
    [~,oname,~] = fileparts(compFile);
    oname = strrep(oname,'coarse','fine');
    oname = fullfile(oDir,[oname '_fineV2']);
    fprintf('Saving:\n\t%s\n',oname)
    save(oname,'cI','cS','MER','pf','rscale')
else
    load(fineFile)
end

%% (3) Interpolate + fine search for arbitrary R, MER, Pf
    % --> Function for array input of R, Pf, or MER, Pf?

if runInterp
%     load(compFile)
    load(fineFile)
    
    Rmax = cS.Rmax;
    Rmin = cS.Rmin;
    Rmax(~cS.success) = NaN;
    [P,Q] = meshgrid(pf,MER);
    
    RmaxN = griddata(P(cS.success),Q(cS.success),Rmax(cS.success),P(~cS.success),Q(~cS.success));
    Rmax(~cS.success) = RmaxN;
    
    RminN = griddata(P(cS.success),Q(cS.success),Rmin(cS.success),P(~cS.success),Q(~cS.success));
    Rmin(~cS.success) = RminN;
    
%     for kk=1:size(Rmax,2)
%         rn = isnan(Rmax(:,kk));
%         if any(rn)
%             Rmax(rn,kk) = interp1(MER(~rn),Rmax(~rn,kk),Rmax(rn,kk),'pchip','extrap');
%         end
%     end
    
    [~,oname,~] = fileparts(fineFile);
    oname = [oname '_rInterpV1'];
    oname = fullfile(oDir,oname);
    fprintf('Saving:\n\t%s\n',oname)
    save(oname,'cI','cS','MER','pf','rscale','Rmax','Rmin')

    % Test vals
%     jj = 3;
%     kk = 1;
%     
%     cIt = cI(jj,kk);
%     cO = Conduit_flow_with_nucleation_V6(cIt);
%     [Zpass,UPpass,Mpass,Ppass,Fcheck,valid] = checkConduitResult(cO,ZfailScale*cIt.conduit_radius,Mfailthresh,Pfailthresh);
%     
%     if valid
%         Rvalid  = cIt.conduit_radius;
%     else
%         Rvalid = [];
%     end
%     [~,Rbi] = mink(abs(ric(:,jj,kk) - cIt.conduit_radius), 2);
%     Rbounds = sort(ric(Rbi,jj,kk));
%     
%     [Rrange,cIi,cO,success] = conduitRadiusFromQ(cIt,Rbounds,'Rvalid',Rvalid,'maxIter',maxIter,...
%         'Zfailthresh',ZfailScale*cIt.conduit_radius,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
%         'dRminScale',dRminScale,'verbose',true,'output',true);

%     G = griddedInterpolant(find(cS.success),Rmax(cS.success))
%     Rmax = interp2(P,Q,Rmax); % cS.pf, cS.
    
end
%% Plot (sample) sweep results
% 3d array: 1st dim: rscale, 2nd dim: MER, 3rd dim: pf

if runPlots
    % load(fineFile)
    pk = length(pf); % Which pressure to plot

    % Show min depth and mach number for all MER and rscale for 2 pf vals
    figure
    subplot(2,2,1)
    imagesc(log10(MER),rscale,Ztc(:,:,1)/1e3)
    set(gca,'YDir','normal')
    xlabel('log_{10}(MER)')
    ylabel('rscale')
    title(sprintf('P = %.2f MPa',pf(1)/1e6))
    cb=colorbar;
    cb.Label.String = 'Z_{top} (km)';

    subplot(2,2,2)
    imagesc(log10(MER),rscale,Mc(:,:,1))
    set(gca,'YDir','normal')
    xlabel('log_{10}(MER)')
    ylabel('rscale')
    title(sprintf('P = %.2f MPa',pf(1)/1e6))
    cb=colorbar;
    cb.Label.String = 'M';

    subplot(2,2,3)
    imagesc(log10(MER),rscale,Ztc(:,:,pk)/1e3)
    set(gca,'YDir','normal')
    xlabel('log_{10}(MER)')
    ylabel('rscale')
    title(sprintf('P = %.2f MPa',pf(pk)/1e6))
    cb=colorbar;
    cb.Label.String = 'Z_{top} (km)';

    subplot(2,2,4)
    imagesc(log10(MER),rscale,Mc(:,:,pk))
    set(gca,'YDir','normal')
    xlabel('log_{10}(MER)')
    ylabel('rscale')
    title(sprintf('P = %.2f MPa',pf(pk)/1e6))
    cb=colorbar;
    cb.Label.String = 'M';
end

%% Show initial course search over R with refined bounds for a sample MER,pf
if runPlots
    pk = 21;
    iMER = 1; %qi;

    % Plot Mach#, ZMin depth, and Overpressure vs radius for up to 3 MER and 1 pf
    figure
    sp(1)=subplot(4,1,1);
    plot(ric(:,iMER,pk),Ztc(:,iMER,pk)/1e3,'.-')
    hold on
    set(gca,'ColorOrderIndex',1)
    % plot([1;1]*cS.Rmax(iMER,pk)',[0;max(Ztc(:))]*[1 1],'--')
    set(gca,'ColorOrderIndex',1)
    % plot([1;1]*cS.Rmin(iMER,pk)',[0;max(Ztc(:))]*[1 1],':')
    scatter(cS.Rmax(iMER,pk),cS.Ztop(iMER,pk)/1e3,50,'r')
    ylabel('Z_{top} (km)')
    axis tight
    ylim([0 .1])

    sp(2) = subplot(4,1,2);
    plot(ric(:,iMER,pk),Mc(:,iMER,pk),'.-')
    hold on
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmax(iMER,pk)',[0;1]*[1 1],'--')
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmin(iMER,pk)',[0;1]*[1 1],':')
    scatter(cS.Rmax(iMER,pk),cS.M(iMER,pk),50,'r')
    ylabel('M')
    legend({'7.5','8'})
    axis tight
    title(sprintf('MER = %.3e kg/s, P_f = %.3f MPa',MER(iMER),pf(pk)/1e6))

    sp(3)=subplot(4,1,3);
    plot(ric(:,iMER,pk),pmc(:,iMER,pk)./pfc(:,iMER,pk),'.-');
    hold on
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmax(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],'--')
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmin(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],':')
    scatter(cS.Rmax(iMER,pk),cS.pm(iMER,pk)./pf(pk),50,'r')
    ylabel('K')
    xlabel('a (m)')
    axis tight

    sp(4) = subplot(4,1,4);
    plot(ric(:,iMER,pk),avc(:,iMER,pk)./ric(:,iMER,pk),'.-');
    hold on
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmax(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],'--')
    set(gca,'ColorOrderIndex',1)
    plot([1;1]*cS.Rmin(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],':')
    scatter(cS.Rmax(iMER,pk),cS.pm(iMER,pk)./pf(pk),50,'r')
    ylabel('a_v/a_0')
    xlabel('a (m)')
    axis tight

    linkaxes(sp,'x')
end

%% Plot up the parameter space
if runPlots
    fs = 12;
    logQ = log10(MER);
    dx = 0.04;
    dy = 0.08;
    ppads = [.04 .03 0.8 0.05];

    figure('position',[50 50 1800 1200])
    % success, Rmax, Rmin, Ztop, M, K, rho
    subplot(2,4,1)
    imagesc(MER,pf/1e6,double(cS.success'))
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'FontSize',fs)
    set(gca,'Xscale','log')
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Success';

    subplot(2,4,2)
    imagesc(MER,pf/1e6,cS.Rmax')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'R_{max} (m)';

    subplot(2,4,3)
    imagesc(MER,pf/1e6,cS.Rmin')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'R_{min} (m)';

    subplot(2,4,4)
    imagesc(MER,pf/1e6,cS.Ztop')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Z_{min} (m)';
    caxis([0 max(cS.Rmax(:))*2])

    subplot(2,4,5)
    imagesc(MER,pf/1e6,cS.U')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Velocity (m/s)';

    subplot(2,4,6)
    imagesc(MER,pf/1e6,cS.M')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Mach #';

    subplot(2,4,7)
    imagesc(MER,pf/1e6,log10(cS.pm./cS.pf)')
    % view([0 0 1])
    % axis tight
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Overpressure';

    subplot(2,4,8)
    imagesc(MER,pf/1e6,cS.rho_b')
    shading flat
    colormap jet
    set(gca,'YDir','normal')
    set(gca,'Xscale','log')
    set(gca,'FontSize',fs)
    xlabel('MER (kg/s)')
    ylabel('P_a (MPa)')
    cb = colorbar;
    cb.Label.String = 'Density (kg/m^3)';
end
%% 
% Rmax (and Rmin) vs pf for 1 or a few MERS
% MER vs pf for 1 or a few Rmax?
% if size(O.Rmax,1)>3
%     iMER = round(linspace(1,size(ric,2),3));
% else
%     iMER = 1:size(O.Rmax,1);
% end

% figure
% subplot
% plot(
%%
% Rmax (and Rmin) vs MER for a few pf's

if runPlots
    if size(cS.Q,2)>3
        ip = round(linspace(1,size(cS.Q,2),3));
    else
        ip = 1:size(cS.Q,2);
    end
    pi = pf(ip);
    for ii=1:length(ip)
        pl{ii} = sprintf('P_a = %.2f MPa',pi(ii)/1e6);
    end
    figure
    rp = cS.Rmax(:,ip);
    semilogx(MER,rp,'--');
    rp(~cS.success(:,ip))=NaN;
    hold on
    set(gca,'ColorOrderIndex',1)
    pp=semilogx(MER,rp,'.-','LineWidth',1.5);
    legend(pp,pl)
end
%% Show Rmax surface as f'n of MER, PEF
if runPlots
    figure
    surf(MER,pf/1e6,cS.Rmax','FaceAlpha',0)
    cSgood = cS.Rmax;
    cSgood(~cS.success) = NaN;
    hold on
    surf(MER,pf/1e6,cSgood','EdgeAlpha',0)
    set(gca,'XScale','log')
    xlabel('MER')
    ylabel('P_a')
    cb=colorbar;
end