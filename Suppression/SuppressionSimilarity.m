function SuppressionSimilarity(raw,iFsuppFigs,iFsuppPlots)
% SuppressionSimilarity - shows how similar different suppression data are
%   SuppressionSimilarity(raw,iFsuppFigs,iFsuppPanels) takes the second
%   (pooled) output from SuppressionProfile, raw, and plots response
%   curves. Each figure shows in all panels the same trace for a single
%   Fsupp at 40 dB SPL, with Fsupp determined by input iFsuppFigs. Each
%   panel shows for a different Fsupp all responses of different SPLs. This
%   way, different suppression-curve shapes can be compared.
%
%   The number of figures generated is equal to numel(iFsuppFigs). The
%   number of panels per figure is equal to numel(iFsuppPanels) and should
%   be equal to 4 or 9.
%   If iFsuppFigs and/or iFsuppPanels are left empty, 5 figures with 4
%   panels are generated. The Fsupp used are evenly distributed through the
%   response.
%
%   See also SuppressionDisplacement, SuppressionProfile, SuppressionSensitivity,
%   SuppressionSPL, SuppressionStruct, SuppressionTuning.
%

% Number of Fsupp
allFsupp = [raw.Fprofile];
uFsupp = unique(allFsupp);
nFsupp = numel(uFsupp);

% Defaults: select random Fprofiles for figures and subplots
% Figures
if nargin < 2 || isempty(iFsuppFigs)
    nFigs = 5;
    iFsuppFigs = unique(round(1:nFsupp/nFigs:nFsupp));
end
if nargin < 3 || isempty(iFsuppPlots)
    % Subplots
    nSubplots = 4;
    iFsuppPlots = unique(round(1:nFsupp/nSubplots:nFsupp));
else
    nSubplots = numel(iFsuppPlots);
    if ~isequal(nSubplots,4) && ~isequal(nSubplots,9)
        error('Number of queried subplots should be 4 or 9.')
    end
end

% Create matrices with gain for each Fsupp (Fprofile)
for ii = 1:nFsupp
    cRaw = raw(allFsupp == uFsupp(ii));
    S(ii).Fsupp = uFsupp(ii)/1e3;
    S(ii).Fprim = cat(1,cRaw.Fprim).'/1e3;
    S(ii).Gain = (cat(1,cRaw.Gain)+pmask(cat(1,cRaw.qAlpha))).';
    S(ii).SPL = [cRaw.SPLjump];
end

% Plotting
for ii = 1:numel(iFsuppFigs)
    Sfig = S(iFsuppFigs(ii)); % S for this figure
    figure('Position',[1 31 1680 953])
    for jj = 1:numel(iFsuppPlots)
        Splot = S(iFsuppPlots(jj)); % S for this subplot
        subplot(sqrt(nSubplots),sqrt(nSubplots),jj)
        hold on
        iSPL = round((1+numel(Sfig.SPL))/2);
        for kk = iSPL
            plot(Sfig.Fprim(:,kk),Sfig.Gain(:,kk)+94,'ko-','LineWidth',2)
        end
        hpl = plot(Splot.Fprim,Splot.Gain+94,'*-');
        title(['Compare against Fsupp = ' num2str(Splot.Fsupp,3) ' kHz'])
    end
    for jj = 1:numel(Splot.SPL), legstr{jj} = [int2str(Splot.SPL(jj)+raw(1).baseSPL) ' dB SPL']; end
    legend(hpl,legstr)
    % Set labels using secondary axes
    hax = axes('Visible','off','FontSize',14);
    pos = get(hax,'Position');
    set(hax,'Position',pos.*[1 1 1 1.04])
    ylabel(hax,'Magnitude re 1 mm/s/Pa (dB)','Visible','on');
    xlabel(hax,'F_{probe} (kHz)','Visible','on');
    title(hax,[raw(1).ExpName ' -- CF ' num2str(raw(1).CF,3) ' kHz -- Fsupp ' num2str(Sfig.Fsupp,3) ' kHz @' int2str(Sfig.SPL(iSPL)) ' dB SPL'],'Visible','on');
end
