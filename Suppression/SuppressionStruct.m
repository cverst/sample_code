function D = SuppressionStruct(DSflat,DSsupp,DSref,DSrefX,ME,CF)
% SuppressionStruct - processes suppression data for further use
%   D = SuppressionStruct(DSflat,DSsupp,DSref,DSrefX,ME,CF) processes the
%   data for further analysis.
%
%   INPUTS
%       DSflat  dataset array containing flat zwuis
%       DSsupp  struct with dataset arrays containing suppression data*
%       DSref   dataset array containing reference data (single references)*
%       DSrefX  dataset array containing reference data (multiple reference)
%       ME      middle ear reference
%       CF      characteristic frequency in kHz
%
%   * The order of the dataset arrays in the fields of DSsupp should match
%   the order of the dataset array DSref.
%
%
%  OUTPUT D has fields
%
%           ExpName: experiment name
%              irec: (pooled) recording number
%          irec_ref: (pooled) recording number of associated reference (multiple references)
%         irec_refX: (pooled) recording number of associated reference (single reference)
%                CF: characteristic frequency
%      is_flat_data: logical whether data is flat zwuis data
%      is_supp_data: logical whether data is suppression zwuis data
%       is_ref_data: logical whether data is reference data
%                '': ''
%             Fprim: stimulus components
%             Fsupp: suppressor frequency
%               iCF: index to characteristic frequency
%             isupp: index to suppressor frequency
%       Fsupp_eq_CF: logical whether Fsupp == CF
%           baseSPL: base sound intensity
%           SPLjump: jump in sound intensity
%             Lsupp: total sound intensity of suppressor
%            allSPL: vector with sound intensity per stimulus component
%                '': ''
%              Gain: gain per stimulus component re stimulus
%             Phase: phase per stimulus component re stimulus
%           Gain_ME: gain per stimulus component re middle ear
%          Phase_ME: phase per stimulus component re middle ear
%             Alpha: Rayleigh significance per stimulus component
%             Rcrit: used Rayleigh criterion
%           Rsignif: logical to Rayleigh-significant stimulus components
%               Vel: BM velocity per stimulus component
%              Disp: BM displacement per stimulus component
%      SignVelPower: total BM velocity
%     SignDispPower: total BM displacement
%                '': ''
%          Gain_ref: gain per stimulus component re reference (single)
%         Phase_ref: phase per stimulus component re reference (single)
%         Alpha_ref: Rayleigh significance per stimulus component of the reference (single)
%       Rsignif_ref: logical to Rayleigh-significant stimulus components of the reference (single)
%         Gain_refX: gain per stimulus component re reference (multiple)
%        Phase_refX: phase per stimulus component re reference (multiple)
%        Alpha_refX: Rayleigh significance per stimulus component of the reference (multiple)
%      Rsignif_refX: logical to Rayleigh-significant stimulus components of the reference (multiple)
%                '': ''
%      PlateauStart: first element of Fprim that has response in plateau
%          qPlateau: logical index to plateau response
%
%   See also SuppressionDisplacement, SuppressionProfile, SuppressionSensitivity,
%   SuppressionSimilarity, SuppressionSPL, SuppressionTuning.
%

% Parameters
Rcrit = .001;

% Analyze flat data
F = local_flat(DSflat,ME,CF,Rcrit);
F = sortAccord(F, [F.baseSPL]);

% Analyze suppression and reference data
[S R] = local_supp_ref(DSsupp,DSref,DSrefX,ME,CF,Rcrit);
S = sortAccord(S, [S.Fsupp] + 0.1*[S.baseSPL]+0.001*[S.Lsupp]);
R = sortAccord(R, [R.baseSPL]);

% Concatenate structs
D = [F S R];

% Indicate HF-plateau
indx = SuppressionPlateau(D);
for ii = 1:numel(D)
    D(ii).d_______________ = '________________';
    D(ii).PlateauStart = indx(ii);
    if ~isnan(indx(ii))
        D(ii).qPlateau = [false(1,indx(ii)-1) true(1,numel(D(ii).Fprim)-indx(ii)+1)];
    else
        D(ii).qPlateau = false(1,numel(D(ii).Fprim));
    end
end


function S = local_flat(DS,ME,CF,Rcrit)
% local_flat - analyzes flat zwuis data

% Apple data
QQ = apple(DS,2,0);

% Build struct
for ii = 1:numel(QQ)
    S(ii).ExpName = QQ(ii).ExpName;
    S(ii).irec = QQ(ii).irec_pool;
    S(ii).irec_ref = {};
    S(ii).irec_refX = {};
    S(ii).CF = CF;
    S(ii).is_flat_data = true;
    S(ii).is_supp_data = false;
    S(ii).is_ref_data = false;
    %
    S(ii).a_______________ = '________________';
    S(ii).Fprim = QQ(ii).Fprim;
    S(ii).Fsupp = [];
    S(ii).iCF = local_findfreq(QQ(ii).Fprim/1e3,CF);
    S(ii).isupp = [];
    S(ii).Fsupp_eq_CF = [];
    S(ii).baseSPL = QQ(ii).baseSPL;
    S(ii).SPLjump = 0;
    S(ii).Lsupp = [];
    S(ii).allSPL = ones(size(S(ii).Fprim))*S(ii).baseSPL;
    %
    S(ii).b_______________ = '________________';
    S(ii).Gain = QQ(ii).Gain + 94;
    S(ii).Phase = delayPhase(QQ(ii).Phase, QQ(ii).Fprim, 0, 2); % unwrap & set start to zero
    S(ii).Gain_ME = QQ(ii).Gain - a2db(abs(eval(ME,QQ(ii).Fprim))) + 94; % reference to ME
    S(ii).Phase_ME = delayPhase(QQ(ii).Phase - cangle(eval(ME,QQ(ii).Fprim)), QQ(ii).Fprim, 0, 2); % reference to ME, unwrap & set start to zero
    S(ii).Alpha = QQ(ii).Alpha;
    S(ii).Rcrit = Rcrit;
    S(ii).Rsignif = QQ(ii).Alpha <= Rcrit;
    S(ii).Vel = QQ(ii).Gain + S(ii).allSPL;
    S(ii).Disp = QQ(ii).Gain + S(ii).allSPL - a2db(2*pi) - a2db(QQ(ii).Fprim);
    S(ii).SignVelPower = p2db(nansum(db2p(S(ii).Vel) + pmask(S(ii).Rsignif)));
    S(ii).SignDispPower = p2db(nansum(db2p(S(ii).Disp) + pmask(S(ii).Rsignif)));
    %
    S(ii).c_______________ = '________________';
    S(ii).Gain_ref = [];
    S(ii).Phase_ref = [];
    S(ii).Alpha_ref = [];
    S(ii).Rsignif_ref = [];
    S(ii).Gain_refX = [];
    S(ii).Phase_refX = [];
    S(ii).Alpha_refX = [];
    S(ii).Rsignif_refX = [];
end


function [S D] = local_supp_ref(DSsupp,DSref,DSrefX,ME,CF,Rcrit)
% local_supp_ref - analyzes flat zwuis data

% Use recursion to have appropriate link between DSsupp and DSref
nDSref = numel(DSref);
S = [];
if isstruct(DSsupp)
    FN = fieldnames(DSsupp);
    nFN = numel(FN);
    for ii = 1:nFN
        [SS D] = local_supp_ref(DSsupp.(FN{ii}),DSref,DSrefX,ME,CF,Rcrit); % directly returns D; otherwise DSref keeps "repeating" itself
        S = [S SS];
    end
    return
elseif nDSref > 1
    D = [];
    for ii = 1:nDSref
        [SS DD] = local_supp_ref(DSsupp(ii),DSref(ii),DSrefX,ME,CF,Rcrit);
        S = [S SS];
        D = [D DD];
        if ii > 1 % remove doublures
            if D(end).irec{1} == D(end-2).irec{1}
                D(end-2) = [];
            end
        end
    end
    return
end

% Apple data
QQ = apple(DSsupp,2,0);
RR = apple(DSref,2,0);
RX = apple(DSrefX,2,0);

% Build suppression struct
for ii = 1:numel(QQ)
    S(ii).ExpName = QQ(ii).ExpName;
    S(ii).irec = QQ(ii).irec_pool;
    S(ii).irec_ref = RR.irec_pool;
    S(ii).irec_refX = RX.irec_pool;
    S(ii).CF = CF;
    S(ii).is_flat_data = false;
    S(ii).is_supp_data = true;
    S(ii).is_ref_data = false;
    %
    S(ii).a_______________ = '________________';
    S(ii).Fprim = QQ(ii).Fprim;
    S(ii).Fsupp = QQ(ii).Fprofile;
    S(ii).iCF = local_findfreq(QQ(ii).Fprim/1e3,CF);
    S(ii).isupp = local_findfreq(QQ(ii).Fprim,QQ(ii).Fprofile);
    S(ii).Fsupp_eq_CF = S(ii).isupp == S(ii).iCF;
    S(ii).baseSPL = QQ(ii).baseSPL;
    S(ii).SPLjump = QQ(ii).SPLjump;
    S(ii).Lsupp = QQ(ii).baseSPL + QQ(ii).SPLjump;
    S(ii).allSPL = ones(size(S(ii).Fprim))*S(ii).baseSPL; S(ii).allSPL(S(ii).isupp) = S(ii).Lsupp;
    %
    S(ii).b_______________ = '________________';
    S(ii).Gain = QQ(ii).Gain + 94;
    S(ii).Phase = delayPhase(QQ(ii).Phase, QQ(ii).Fprim, 0, 2); % unwrap & set start to zero
    S(ii).Gain_ME = QQ(ii).Gain - a2db(abs(eval(ME,QQ(ii).Fprim))) + 94; % reference to ME
    S(ii).Phase_ME = delayPhase(QQ(ii).Phase - cangle(eval(ME,QQ(ii).Fprim)), QQ(ii).Fprim, 0, 2); % reference to ME, unwrap & set start to zero
    S(ii).Alpha = QQ(ii).Alpha;
    S(ii).Rcrit = Rcrit;
    S(ii).Rsignif = QQ(ii).Alpha <= Rcrit;
    S(ii).Vel = QQ(ii).Gain + S(ii).allSPL;
    S(ii).Disp = QQ(ii).Gain + S(ii).allSPL - a2db(2*pi) - a2db(QQ(ii).Fprim);
    S(ii).SignVelPower = p2db(nansum(db2p(S(ii).Vel) + pmask(S(ii).Rsignif)));
    S(ii).SignDispPower = p2db(nansum(db2p(S(ii).Disp) + pmask(S(ii).Rsignif)));
    %
    S(ii).c_______________ = '________________';
    S(ii).Gain_ref = QQ(ii).Gain - RR.Gain;
    S(ii).Phase_ref = delayPhase(QQ(ii).Phase - RR.Phase, QQ(ii).Fprim, 0, 2);
    S(ii).Alpha_ref = RR.Alpha;
    S(ii).Rsignif_ref = RR.Alpha <= Rcrit;
    S(ii).Gain_refX = QQ(ii).Gain - RX.Gain;
    S(ii).Phase_refX = delayPhase(QQ(ii).Phase - RX.Phase, QQ(ii).Fprim, 0, 2);
    S(ii).Alpha_refX = RX.Alpha;
    S(ii).Rsignif_refX = RX.Alpha <= Rcrit;
end

% Build reference struct
D.ExpName = RR.ExpName;
D.irec = RR.irec_pool;
D.irec_ref = [];
D.irec_refX = [];
D.CF = CF;
D.is_flat_data = false;
D.is_supp_data = false;
D.is_ref_data = true;
%
D.a_______________ = '________________';
D.Fprim = RR.Fprim;
D.Fsupp = [];
D.iCF = local_findfreq(RR.Fprim/1e3,CF);
D.isupp = [];
D.Fsupp_eq_CF = [];
D.baseSPL = RR.baseSPL;
D.SPLjump = 0;
D.Lsupp = [];
D.allSPL = ones(size(D.Fprim))*D.baseSPL;
%
D.b_______________ = '________________';
D.Gain = RR.Gain + 94;
D.Phase = delayPhase(RR.Phase, RR.Fprim, 0, 2); % unwrap & set start to zero
D.Gain_ME = RR.Gain - a2db(abs(eval(ME,RR.Fprim))) + 94; % reference to ME
D.Phase_ME = delayPhase(RR.Phase - cangle(eval(ME,RR.Fprim)), RR.Fprim, 0, 2); % reference to ME, unwrap & set start to zero
D.Alpha = RR.Alpha;
D.Rcrit = Rcrit;
D.Rsignif = RR.Alpha <= Rcrit;
D.Vel = RR.Gain + D.allSPL;
D.Disp = RR.Gain + D.allSPL - a2db(2*pi) - a2db(RR.Fprim);
D.SignVelPower = p2db(nansum(db2p(D.Vel) + pmask(D.Rsignif)));
D.SignDispPower = p2db(nansum(db2p(D.Disp) + pmask(D.Rsignif)));
%
D.c_______________ = '________________';
D.Gain_ref = [];
D.Phase_ref = [];
D.Alpha_ref = [];
D.Rsignif_ref = [];
D.Gain_refX = [];
D.Phase_refX = [];
D.Alpha_refX = [];
D.Rsignif_refX = [];

% Build referenceX struct
D(2).ExpName = RX.ExpName;
D(2).irec = RX.irec_pool;
D(2).irec_ref = [];
D(2).irec_refX = [];
D(2).CF = CF;
D(2).is_flat_data = false;
D(2).is_supp_data = false;
D(2).is_ref_data = true;
%
D(2).a_______________ = '________________';
D(2).Fprim = RX.Fprim;
D(2).Fsupp = [];
D(2).iCF = local_findfreq(RX.Fprim/1e3,CF);
D(2).isupp = [];
D(2).Fsupp_eq_CF = [];
D(2).baseSPL = RX.baseSPL;
D(2).SPLjump = 0;
D(2).Lsupp = [];
D(2).allSPL = ones(size(D(2).Fprim))*D(2).baseSPL;
%
D(2).b_______________ = '________________';
D(2).Gain = RX.Gain + 94;
D(2).Phase = delayPhase(RX.Phase, RX.Fprim, 0, 2); % unwrap & set start to zero
D(2).Gain_ME = RX.Gain - a2db(abs(eval(ME,RX.Fprim))) + 94; % reference to ME
D(2).Phase_ME = delayPhase(RX.Phase - cangle(eval(ME,RX.Fprim)), RX.Fprim, 0, 2); % reference to ME, unwrap & set start to zero
D(2).Alpha = RX.Alpha;
D(2).Rcrit = Rcrit;
D(2).Rsignif = RX.Alpha <= Rcrit;
D(2).Vel = RX.Gain + D(2).allSPL;
D(2).Disp = RX.Gain + D(2).allSPL - a2db(2*pi) - a2db(RX.Fprim);
D(2).SignVelPower = p2db(nansum(db2p(D(2).Vel) + pmask(D(2).Rsignif)));
D(2).SignDispPower = p2db(nansum(db2p(D(2).Disp) + pmask(D(2).Rsignif)));
%
D(2).c_______________ = '________________';
D(2).Gain_ref = [];
D(2).Phase_ref = [];
D(2).Alpha_ref = [];
D(2).Rsignif_ref = [];
D(2).Gain_refX = [];
D(2).Phase_refX = [];
D(2).Alpha_refX = [];
D(2).Rsignif_refX = [];


function ifreq = local_findfreq(allfreqs,freq)
% local_findfreq - finds index to specific frequency

[junk ifreq] = min(abs(allfreqs/freq-1));
