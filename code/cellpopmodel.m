function [Tmat, Emat, CD20samps, CD16samps, CPXvec, LDHvec, perfvec] = cellpopmodel(CD20dist0, CD16dist0,tvec, T0, E0,...
    adcxpars,CPXpars)
% This function is going to call the adcx function for an individual value
% of CD20 and CD16 that is sampled from the observed distribution of a
% patient. Within the adcx function, the reaction_ss function will be
% called, which, for a given CD20, CD16, R, kd20 and kd16, outputs the
% number of complexes that are formed. In order to determine the rate at
% which effector cells (E) kill tumor cells (T) and produce depleted E
% cells (Estar) and dead cells (D which produce perf and LDH), we assume
% this is a function of the number of complexes formed (i.e. the output of
% reactios_ss), linearly proportional to some parameter gamma. Then, for a
% given CD20, CD16 level, we get a number of E cells in time, number of
% tumor cells in time, and number of E star in time. 
%
%The CD16 level of the E star cells is:
% CD16Estar = CD16E- CPX
% giving us two locations on our histogram for CD16 from each forward run
% The CD20 level fo the T cells is still the same, but over time the number
% of tumor cells changes.

% OUTPUTS:
% Tmat: matrix of length t and width of the sampled CD20 values containing
% number of tumor cells T in time for each CD20 value
% Emat: matrix of length t and width of the sampled CD16 values containing
% number of tumor cells E in time for each CD16 value
% Cpxvec: row vector of width of sampled CD16/CD20values containing
% the resulting number of complexes in time for each CD16/CD20 value
% CD20samps: these are the values of CD20 that were sampled that correspond 
% to the columns in your Tmat
% CD16samps: these are the values of CD16 that were sampled that correspond
% the columns in your E mat
% LDHvec: column vector of length t containing the sum of the number of
% dead cells at each time for all CD16/CD20 combos, times some LDH
% production factor
% perfvec: column vector of length t containing the sum of the number of
% effector cells at each time for all CD16/CD20 combos, times some perf
% production factor

% INPUTS:
% CD20dist0: initial matlab distribution object where we should be able to either
% randomly sample from it, or sample in steps according to the probability
% obtained from data
% CD16dist0: initial matlab distribution object also obtained from data
% tvec: input time vector that we want to model (assuming this is in hours)
% T0: initial number of tumor cells (T0= T+D)
% E0: initial number of effector cells (E0= E+E*) assume effector cell
% total is constant
% adcvpars: R, kd20, kd16 (RTX concentration, rate of CD20 binding, rate of
% CD16binding)
% CPXpars: g, r, kexp (tumor growth rate, stoichiometry of E cells needed
% to kill 1 tumor cell, rate of return of E* cells to E cells (replenishing
% of CD16 after shedding)

nsamps = 1000;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);

for i = 1:nsamps
    % generate one column of your sampled data by calling the adcx model
    % which calles the reaction_ss function
    [Ti, Ei, Estari, LDHi, perfi, CPXi] = adcx(tvec, T0,E0, adcxpars,...
        CD20samps(i), CD16samps(i),CPXpars);
    Tmat(:,i) = Ti;
    Emat(:,i) = Estari;
    LDHmat(:,i) = LDHi;
    perfmat(:,i) =perfi;
    delCD16(i,1) = CPXi;
    
end
% For a given CD16, CD20 distribution, here are the expected LDH and perf
% trajectories over time
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);


end

