function SuppressionDisplacement(D,iFprobe)
% SuppressionDisplacement - plots probe displacement vs. supp. displacement
%   SuppressionDisplacement(D) takes the output from SuppressionStruct and
%   plots probe displacement vs. suppressor displacement. Each curve
%   describes a different probe frequency and each panel shows a different
%   suppressor frequency.
%   SuppressionSensitivity(D,iFprobe) limits the number of plots to only
%   the Fprobes specified by iFprobe.
%
%   See also SuppressionProfile, SuppressionSensitivity, SuppressionSimilarity,
%   SuppressionSPL, SuppressionStruct, SuppressionTuning.
%

% Default
if nargin < 2, iFprobe = []; end

% Reduce data
Dflat = D([D.is_flat_data]);
Dref = D([D.is_ref_data]);
D = D([D.is_supp_data]);

% Number of plots
allFprobe = cat(1,D(1).Fprim); % assuming that recordings were done properly and reference Fprobe (and Fprobe) are always identical
nFprobe = numel(allFprobe);

% Number of curves
allFsupp = [D.Fsupp];
uFsupp = unique(allFsupp);
nFsupp = numel(uFsupp);

% Data preparation
FsuppCF = cat(1,D.Fsupp_eq_CF);
DD = cat(1,D.Disp); % displacement
PH = cat(1,D.Phase); % phase
allisupp = cat(1,D.isupp);
allPlateauStart = cat(1,D.PlateauStart);
qsign = cat(1,D.Rsignif); % rayleigh-significant
msk = pmask(qsign);% + pmask(~cat(1,D.qPlateau)); % and not in pateau...
DDprobe = DD+msk; % significant displacement
PHprobe = PH+msk; % significant phase
for ii = 1:numel(D)
    DDsupp(ii) = D(ii).Disp(D(ii).isupp); %#ok<AGROW> % suppressor displacements
    DDprobe(ii,D(ii).isupp) = NaN; % set suppressor displacements to NaN
    PHprobe(ii,D(ii).isupp) = NaN; % set suppressor phase to NaN
end
% Probe alone
PA = Dflat([Dflat.baseSPL] == D(1).baseSPL);
PAdisp = PA.Disp+pmask(PA.Rsignif);

% Plotting
if isempty(iFprobe), iter = 1:nFprobe; else iter = iFprobe; end
count = 1;
for ii = iter
    % Figure generation
    if mod(count-1,9) == 0 % create figure when necessary
        figure('Position',[1 31 1680 953])
    end
    cFprobe = allFprobe(ii); % current Fprobe
    cDDprobe = DDprobe(:,ii); % displacements for current Fprobe
    hpl(count) = subplot(3,3,mod(count-1,9)+1); %#ok<AGROW>
    hold on
    title(['F_{probe} = ' num2str(cFprobe/1e3,3)])
    for jj = 1:nFsupp
        qFsupp = allFsupp == uFsupp(jj);
        cisupp = unique(allisupp(qFsupp)); % index to current suppressor
        cPlateauStart = allPlateauStart(qFsupp); % current start index of plateau
        notPlateau = cisupp < cPlateauStart; % logical index to non-plateau values
        if ~any(FsuppCF(qFsupp)), lwdth = 1; else lwdth = 3; CFsupp = uFsupp(jj); end
        plot(DDsupp(qFsupp)+120,cDDprobe(qFsupp)+120+pmask(notPlateau),'color',local_color(uFsupp(jj)),'LineWidth',lwdth)
    end
    set(hpl(count),'UserData',uFsupp)
    % Probe alone
    pPA = db2p(PAdisp(ii)+120);
    paX = p2db(linspace(db2p(-100),pPA,1e3));
    paY = p2db(pPA-db2p(paX));
    paY(end) = -400; % set -inf to -400
    plot(paX,paY,'k:')
    % Axis labelling
    if mod(count,9) == 0 || ii == iter(end)
        local_colorbar(uFsupp,CFsupp)
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'Probe displacement (dB re 1 nm)','Visible','on');
        xlabel(hax,'Suppressor displacement (dB re 1 nm)','Visible','on');
        title([D(1).ExpName ' -- DISPLACEMENT'],'Visible','on')
        uistack(hax,'bottom')
    end
    count = count+1;
end
set(hpl,'YLim',[-70 -5],'XLim',[-60 50])

% Plotting PHASE
count = 1;
for ii = iter
    % Figure generation
    if mod(count-1,9) == 0 % create figure when necessary
        figure('Position',[1 31 1680 953])
    end
    cFprobe = allFprobe(ii); % current Fprobe
    cPHprobe = PHprobe(:,ii); % phase for current Fprobe
    hpl(count) = subplot(3,3,mod(count-1,9)+1); %#ok<AGROW>
    hold on
    title(['F_{probe} = ' num2str(cFprobe/1e3,3)])
    for jj = 1:nFsupp
        qFsupp = allFsupp == uFsupp(jj);
        cisupp = unique(allisupp(qFsupp)); % index to current suppressor
        cPlateauStart = allPlateauStart(qFsupp); % current start index of plateau
        notPlateau = cisupp < cPlateauStart; % logical index to non-plateau values
        if ~any(FsuppCF(qFsupp)), lwdth = 1; else lwdth = 3; CFsupp = uFsupp(jj); end
        ph = cunwrap(cPHprobe(qFsupp));
        plot(DDsupp(qFsupp)+120,ph-round(nanmean(ph))+pmask(notPlateau),'color',local_color(uFsupp(jj)),'LineWidth',lwdth)
    end
    % Axis labelling
    if mod(count,9) == 0 || ii == iter(end)
        local_colorbar(uFsupp,CFsupp)
        % Set labels using secondary axes
        hax = axes('Visible','off','FontSize',14);
        pos = get(hax,'Position');
        set(hax,'Position',pos.*[1 1 1 1.04])
        ylabel(hax,'Probe phase (cycle)','Visible','on');
        xlabel(hax,'Suppressor displacement (dB re 1 nm)','Visible','on');
        title([D(1).ExpName ' -- PHASE'],'Visible','on')
        uistack(hax,'bottom')
    end
    count = count+1;
end
set(hpl,'YLim',[-.5 .5],'XLim',[-60 50])

% Second figure: RMS displacement vs. Lsupp
figure
for ii = 1:nFsupp
    qFsupp = allFsupp == uFsupp(ii);
    cD = D(qFsupp);
    [Lsupp dspl] = sortaccord([cD.Lsupp],[cD.SignDispPower],[cD.Lsupp]);
    if ~cD(1).Fsupp_eq_CF, lwdth = 1; else lwdth = 3; end
    xplot(Lsupp,dspl+120,'.-','Color',local_color(uFsupp(ii)),'LineWidth',lwdth)
end
xlabel('L_{supp} (dB SPL)')
ylabel('RMS displacement (dB re 1 nm)')
title(D(1).ExpName)

% Third figure: Equivalent SPL vs. RMS displacement
S = SuppressionEffectiveLevel([Dflat D Dref],0);
S = S([S.is_supp_data]);
figure
for ii = 1:numel(S)
    cS = S(ii);
    if ~cS.Fsupp_eq_CF, mrkr = '.'; else mrkr = '^'; end
    xplot(cS.SignDispPower+120,cS.EqSPL,'Color',local_color(cS.Fsupp),'MarkerFaceColor',local_color(cS.Fsupp),'Marker',mrkr)
end
xlabel('RMS displacement (dB re 1 nm)')
ylabel('Equivalent level (dB SPL)')
title(D(1).ExpName)

function C = local_color(value)
% function to get color for each line according to associated frequency

allFreq = 50:30e3;
nFreq = numel(allFreq);
clrs = jet(nFreq);
indx = local_findnearest(allFreq,value);
C = clrs(indx,:);

function indx = local_findnearest(V,val)
% finds the element of a vector closest to a specified value

% Search within V for smallest difference
W = V-val;
[junk indx] = min(abs(W));

function local_colorbar(uFsupp,CF)
% function teo generate colorbar

hcb = colorbar('Location','East','Position',[.936 .11 .016 .815]);
title(hcb,'F_{supp} (kHz)')

nFsupp = numel(uFsupp);
set(hcb,'YLim',[0 nFsupp],'YTick',(1:nFsupp)-.5,'TickLength',[0 0],'YTickLabel',arrayfun(@(F) num2str(F/1e3,3),uFsupp,'UniformOutput',false))
cla(hcb)
for ii = 1:nFsupp
    if uFsupp(ii) ~= CF, lwdth = 1; else lwdth = 3; end
    xplot(hcb,[0 1],[ii ii]-.5,'Color',local_color(uFsupp(ii)),'LineWidth',lwdth*2)
end
