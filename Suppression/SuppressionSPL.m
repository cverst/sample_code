function SuppressionSPL(SP)
% SupSPL - plots (pooled) data of suppression per SPL
%   SupSPL(SP) takes the first (pooled) output from SuppressionProfile, SP,
%   and plots suppression for Fprofile vs. Fprobe per SPL.
%
%   See also SuppressionDisplacement, SuppressionProfile, SuppressionSensitivity,
%   SuppressionSimilarity, SuppressionStruct, SuppressionTuning.
%

% Number of plots
allSPLjump = {SP.SPLjump};
SPLjump = unique(vertcat(allSPLjump{:}));
nSPL = numel(SPLjump);

for ii = 1:nSPL
    % Data preparation
    ind2SPL = cellfun(@(spl) find(spl == SPLjump(ii)),allSPLjump,'UniformOutput',false); % cell array containing indices to the current SPLjump
    ind2Fprofile = find(~cell2mat(cellfun(@isempty,ind2SPL,'UniformOutput',false))); % cell array containing indices to Fprofiles that actually have the current SPLjump
    count = 1;
    for jj = ind2Fprofile
        SPLindx = ind2SPL{jj}; % index to current SPLjump for current Fprofile
        Sup(count,:) = mean(SP(jj).Gain(SPLindx,:),1); % suppression for current SPL and Fprofile || Taking mean here; sometimes an SPL was repeated
        Fsup(count) = unique(SP(jj).Fprofile(SPLindx)); % current Fprofile, now named Fsup || if this goes wrong, lose the recording..
        count = count+1;
    end
    Fprobe = SP(jj).Freq/1e3; % probe frequencies
    
    % Actual plotting
    if mod(ii-1,9) == 0 % create figure when necessary
        figure('Position',[1 31 1680 953])
    end
    subplot(3,3,mod(ii-1,9)+1)
    hold on
    contourf(Fprobe,Fsup/1e3,Sup)
    caxis(gca,[min(SP(jj).GainBounds) max(SP(jj).GainBounds)]) % just take last jj; that should do it
    YL = ylim;
    XL = xlim;
    plot(SP(jj).CF,YL(1)+diff(YL)/40,'w^','MarkerFaceColor','k','MarkerSize',10)
    plot(XL(1)+diff(XL)/40,SP(jj).CF,'w>','MarkerFaceColor','k','MarkerSize',10)
    plot(YL,YL,'w-','LineWidth',2)
    plot(YL,YL,'k-','LineWidth',1)
    title(['L_{sup} ' int2str(SPLjump(ii)+SP(jj).baseSPL(1)) ' dB SPL'])
    if mod(ii,9) == 0 || ii == nSPL
        hcol = colorbar('Location','East','Position',[.936 .11 .016 .815]);
        title(hcol,'Suppression (dB)')
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'F_{sup} (kHz)','Visible','on');
        xlabel(hax,'F_{probe} (kHz)','Visible','on');
        title(hax,[SP(jj).TitleStr{1}(1:7) ' -- ' SP(jj).TitleStr{1}(end-27:end-16) ' -- L_{probe} ' int2str(SP(jj).baseSPL(1)) ' dB SPL'],'Visible','on');
        colormap(flipud(colormap))
    end
end
