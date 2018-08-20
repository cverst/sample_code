function SuppressionTuning(raw,SupCrit)
% SuppressionTuning - plots suppression tuning curves
%   SuppressionTuning(raw,crit) takes the second (pooled) output from
%   SuppressionProfile, raw, and plots suppression tuning curves. Each
%   curve describes a different probe frequency. crit is the criterion
%   below which suppression must be (default = -3 dB).
%
%   See also SuppressionDisplacement, SuppressionProfile, SuppressionSensitivity,
%   SuppressionSimilarity, SuppressionSPL, SuppressionStruct.
%

% Criterion for minimum of suppression in dB
if nargin <2, SupCrit = -3; end

% Issue warning
if ~isequal(raw.Fprim,raw.refFprim)
    warning('Fprim of suppression data and reference data not (always) equal; difference ignored.')
end

% Data preparation
allFsupp = [raw.Fprofile];
uFsupp = unique(allFsupp);
nFsupp = numel(uFsupp);
Fprobe = raw(1).Fprim;
nFprobe = numel(Fprobe); % also number of plots
thrSup = NaN(nFsupp,nFprobe);
% Data interpolation and find last diff gain above supp threshold
warning('off','MATLAB:interp1:NaNinY')
for ii = 1:nFsupp
    cRaw = raw(allFsupp==uFsupp(ii)); % raw data for current Fsupp
    qsign = cat(1,cRaw.qsignif); % rayleigh-significant
    SS = [zeros(1,nFprobe); cat(1,cRaw.Gain)-cat(1,cRaw.refGain)+pmask(qsign)]; % suppression data
    SPL = [cRaw.SPLjump]+[cRaw.baseSPL];
    fSPL = linspace(0,max(SPL),1e3); % fine spaced SPL
    fSS = INTERP1([0 SPL].',SS,fSPL.'); % interpolated suppression
    for jj = 1:nFprobe
        indx = find(fSS(:,jj)>=SupCrit,1,'last');
        if ~isempty(indx)
            thrSup(ii,jj) = fSPL(indx);
        end
    end
end
warning('on','MATLAB:interp1:NaNinY')

% Plotting
for ii = 1:nFprobe
    % Figure generation
    if mod(ii-1,9) == 0 % create figure when necessary
        figure('Position',[1 31 1680 953])
    end

    % Actual plotting
    hpl(ii) = subplot(3,3,mod(ii-1,9)+1);
    hold on
    plot(uFsupp/1e3,thrSup(:,ii),'*-')
    
    % Axis labelling
    if mod(ii,9) == 0 || ii == nFprobe
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,['Suppression threshold @' num2str(SupCrit) ' dB (dB SPL)'],'Visible','on');
        xlabel(hax,'F_{supp} (kHz)','Visible','on');
        title(hax,{[raw(1).ExpName ' -- CF ' num2str(raw(1).CF,3) ' kHz -- ' int2str(raw(1).baseSPL) ' dB SPL -- irec ' strrep(int2str(unique([raw.irec])),'  ',' ')]},'Visible','on');
    end
end

% Axes operations
% YL = cell2mat(get(hpl,'YLim'));
% YL = [min(YL(:,1)) max(YL(:,2))];
YL = [20 80];
set(hpl,'YLim',YL)
for ii = 1:numel(hpl)
    plot(hpl(ii),raw(1).CF,YL(1)+diff(YL)/40,'w^','MarkerFaceColor','k','MarkerSize',10)
    plot(hpl(ii),Fprobe(ii)/1e3,YL(1)+diff(YL)/40,'wo','MarkerFaceColor','k','MarkerSize',10)
    title(hpl(ii),['F_{probe} = ' num2str(Fprobe(ii)/1e3,3) ' kHz'])
end
