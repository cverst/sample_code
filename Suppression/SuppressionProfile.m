function Q = SuppressionProfile(D, GainBounds, PhaseBounds) %#ok<INUSD>
% SuppressionProfile - suppression analysis of ZW data with varied spectrum
%   SuppressionProfile(S, GainBounds, PhaseBounds)
%   Inputs are 
%             S: struct array from SuppressionData
%    GainBounds: array holding boundaries for gain contour plot.
%                Default: -47.5:5:7.5
%                
%   PhaseBounds: array holding boundaries for phase contour plot.
%                Default: -.45:.1:.55
%
%   See also SuppressionDisplacement, SuppressionSensitivity, SuppressionSimilarity,
%   SuppressionSPL, SuppressionStruct, SuppressionTuning.
%

% Only plotting?
if isequal(fieldnames(D).',{'Gain' 'Freq' 'Phase' 'Lsupp' 'TitleStr' 'Gbnds' 'Pbnds' 'CF' 'Fsupp'}), local_plot(D), return, end

% Defaults
if nargin < 2, GainBounds = -47.5:5:7.5; end
if nargin < 2, PhaseBounds = -.45:.1:.55; end

% Parameters
Nfreq = 100; % # freq components for uniform freq spacing

% Select data
S = D([D.is_supp_data]); % suppression data
allFsupp = [S.Fsupp];
uFsupp = unique(allFsupp);
nFsupp = numel(uFsupp);

% Reference
Sref = D([D.is_ref_data]); % reference data
refStr = cellfun(@char,[Sref.irec],'UniformOutput',false);

% Boundary freqs
allFprim = cat(1,S.Fprim);
Fmin = max(allFprim(:,1));
Fmax = min(allFprim(:,end));
Freq = linspace(Fmin,Fmax,Nfreq);

% Process data
warning('off','MATLAB:interp1:NaNinY')
for ii = 1:nFsupp
    cS = S(allFsupp == uFsupp(ii)); % current Fsupp
    frq = cat(1,cS.Fprim);
    Lsupp = samesize(cat(1,cS.Lsupp),frq);
%     msk = pmask(cat(1,cS.Rsignif));
    cSstr = cellfun(@char,[cS.irec_refX],'UniformOutput',false);
    indx = cellfun(@(chs) strmatch(chs,refStr,'exact'),cSstr,'UniformOutput',false);
    indx = [indx{:}];
    msk = pmask(cat(1,cS.Rsignif)) + pmask(~cat(1,cS.qPlateau)) + pmask(~cat(1,Sref(indx).qPlateau));
    gn = cat(1,cS.Gain_refX);
    phs = cat(1,cS.Phase_refX);
    %
    Q(ii).ExpName = S(1).ExpName;
    [XI YI] = meshgrid(Freq,Lsupp(:,1));
    Q(ii).Gain = interp2(frq,Lsupp,gn+msk,XI,YI); % interpolation is not really necessary, but it definitely makes the plots look better
    phasor = exp(2*pi*i*(phs)+msk);
    phasor = interp2(frq,Lsupp,phasor,XI,YI);
    Q(ii).Freq = XI;
    Q(ii).Phase = cunwrap(cangle(phasor),1);
    Q(ii).Lsupp = YI;
    Q(ii).TitleStr = ['F_{supp} = ' num2str(uFsupp(ii)/1e3,3) ' kHz'];
    Q(ii).Gbnds = GainBounds;
    Q(ii).Pbnds = PhaseBounds;
    Q(ii).CF = S(1).CF;
    Q(ii).Fsupp = uFsupp(ii);
end
warning('on','MATLAB:interp1:NaNinY')

local_plot(Q)


function local_plot(Q)
% Plot routine

% Number of plots
nQ = numel(Q);

% Plot gain
for ii = 1:nQ
    cQ = Q(ii);
    if mod(ii-1,9) == 0
        figure('Position',[1 31 1680 953])
    end
    subplot(3,3,mod(ii-1,9)+1)
    contourf(cQ.Freq/1e3,cQ.Lsupp,cQ.Gain,cQ.Gbnds)
    title(cQ.TitleStr)
    caxis(gca,[min(cQ.Gbnds) max(cQ.Gbnds)])
    YL = ylim;
    xplot(cQ.CF,YL(1)+diff(YL)/40,'w^','MarkerFaceColor','k','MarkerSize',7)
    xplot(cQ.Fsupp/1e3,YL(1)+diff(YL)/40,'wo','MarkerFaceColor','k','MarkerSize',7)
    if mod(ii,9) == 0 || ii == nQ
        hcb = colorbar('Location','East','Position',[.936 .11 .016 .815]);
        title(hcb,'Suppression (dB)')
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'Suppressor level (dB SPL)','Visible','on');
        xlabel(hax,'Stimulus frequency (kHz)','Visible','on');
        title([Q(1).ExpName ' -- GAIN'],'Visible','on')
        colormap(flipud(jet(numel(cQ.Gbnds)-1)))
        uistack(hax,'bottom')
    end
end

% Plot phase
for ii = 1:nQ
    cQ = Q(ii);
    if mod(ii-1,9) == 0
        figure('Position',[1 31 1680 953])
    end
    subplot(3,3,mod(ii-1,9)+1)
    contourf(cQ.Freq/1e3,cQ.Lsupp,cQ.Phase,cQ.Pbnds)
    title(cQ.TitleStr)
    caxis(gca,[min(cQ.Pbnds) max(cQ.Pbnds)])
    YL = ylim;
    xplot(cQ.CF,YL(1)+diff(YL)/40,'w^','MarkerFaceColor','k','MarkerSize',7)
    xplot(cQ.Fsupp/1e3,YL(1)+diff(YL)/40,'wo','MarkerFaceColor','k','MarkerSize',7)
    if mod(ii,9) == 0 || ii == nQ
        hcb = colorbar('Location','East','Position',[.936 .11 .016 .815]);
        title(hcb,'Phase change (cycle)','Visible','on');
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'Suppressor level (dB SPL)','Visible','on');
        xlabel(hax,'Stimulus frequency (kHz)','Visible','on');
        title([Q(1).ExpName ' -- PHASE'],'Visible','on')
        colormap(jet(numel(cQ.Pbnds)-1))
        uistack(hax,'bottom')
    end
end
