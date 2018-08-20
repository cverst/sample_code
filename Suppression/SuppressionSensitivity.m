function SuppressionSensitivity(raw)
% SuppressionSensitivity - plots sens. change vs. SPL for suppression data
%   SuppressionSensitivity(raw) takes the second (pooled) output from
%   SuppressionProfile, raw, and plots sensitivity change vs. SPL, where
%   each curve describes a different suppressor frequency and each panel
%   describes a different probe frequency.
%
%   See also SuppressionDisplacement, SuppressionProfile, SuppressionSimilarity,
%   SuppressionSPL, SuppressionStruct, SuppressionTuning.
%

% Issue error/warning
if ~isequal(raw.Fprim,raw.refFprim)
    warning('Fprim of suppression data and reference data not (always) equal; difference ignored.')
end

% Number of plots
allFsupp = [raw.Fprofile];
uFsupp = unique(allFsupp);
nFsupp = numel(uFsupp);

% Number of curves
allFprobe = cat(1,raw(1).refFprim); % assuming that recordings were done properly and reference Fprobe (and Fprobe) are always identical
nFprobe = numel(allFprobe);

% Colormap
clrs = jet(nFprobe); % colors for plotting of different frequencies
for ii = 1:nFsupp
    % Data preparation
    cFsupp = uFsupp(ii); % current Fsupp
    qFsupp = allFsupp == cFsupp; % logical index to current Fsupp's
    cRaw = raw(qFsupp); % raw data for current Fsupp
    GG = cat(1,cRaw.Gain); % gain
    rGG = cat(1,cRaw.refGain); % reference gain
    qsign = cat(1,cRaw.qsignif); % rayleigh-significant
    supp = GG-rGG+pmask(qsign); % suppression
    SPL = [cRaw.SPLjump].'+[cRaw.baseSPL].';
    % Figure generation
    if mod(ii-1,9) == 0 % create figure when necessary
        figure('Position',[1 31 1680 953])
    end
    % Actual plotting
    hpl(ii) = subplot(3,3,mod(ii-1,9)+1);
    hold on
    title(['F_{supp} = ' num2str(cFsupp/1e3,3)])    
    for jj = 1:nFprobe
%         cFprobe = allFprobe(jj); % current probe frequency
        plot(SPL,supp(:,jj),'color',clrs(jj,:))
    end
    % Axis labelling
    if mod(ii,9) == 0 || ii == nFsupp
        hcb = colorbar('Location','East','Position',[.936 .11 .016 .815]);
        colormap(clrs)
        set(hcb,'YTick',(1:5:numel(clrs(:,1)))+.5,'YTickLabel',round(allFprobe(1:5:end)/1e2)/10)
        title(hcb,'F_{probe} (kHz)')
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'Sensitivity change (dB)','Visible','on');
        xlabel(hax,'L_{supp} (dB SPL)','Visible','on');
        title(hax,{[raw(1).ExpName ' -- CF ' num2str(raw(1).CF,3) ' kHz -- ' int2str(raw(1).baseSPL) ' dB SPL -- irec ' strrep(int2str(unique([raw.irec])),'  ',' ')]},'Visible','on');
    end
end
YL = cell2mat(get(hpl,'YLim'));
YL = [min(YL(:,1)) max(YL(:,2))];
set(hpl,'YLim',YL)
