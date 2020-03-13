function [Tmat, Emat, CD20mat, CD16mat, delCD16, LDHmat, perfmat] = cellpopmodel(CD20samps, CD16samps,...
     nsamps, lambda, pstruct, expt)
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
% CD20mat: these are the values of CD20 that were sampled that correspond 
% to the columns in your Tmat
% CD16mat: these are the values of CD16 that were sampled that correspond
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
% tvec: time parameters that adcx takes in
% T0: initial number of tumor cells (T0= T+D)
% E0: initial number of effector cells (E0= E+E*) assume effector cell
% total is constant
% adcvpars: R, kd20, kd16 (RTX concentration, rate of CD20 binding, rate of
% CD16binding)
% CPXpars: g, r, kexp (tumor growth rate, stoichiometry of E cells needed
% to kill 1 tumor cell, rate of return of E* cells to E cells (replenishing
% of CD16 after shedding)

% Assign all the parameters needed for adcx from the pstruct
switch expt
    case 'Z138'
        g = pstruct.g_Z138;
        kon16 = pstruct.kon16_V158; % TOTALLY GUESSED ON THIS
    case 'SUDHL4'
        g = pstruct.g_SUDHL4;
        kon16 = pstruct.kon16_F158; % ALSO A GUESS
end






% Allocate the input parameters
% Padcx = num2cell(adcxpars);
% [g,r,kexp,gamma] = deal(Padcx{:});
% Pcpx = num2cell(CPXpars);
% [RTX,kon20,koff20,kon16,koff16,gamma_perf] = deal(Pcpx{:});
% Ptvec = num2cell(tvec);
% [tf_mol,tf_et,nr_t_mol,nr_t_et] = deal(Ptvec{:});
% E0toT0 = E0/T0;

Tmat       = zeros(pstruct.nr_t_et,nsamps);
Emat       = zeros(pstruct.nr_t_et,nsamps);
LDHmat     = zeros(pstruct.nr_t_et,nsamps);
perfmat    = zeros(pstruct.nr_t_et,nsamps);
delCD16    = zeros(nsamps,1);


for i = 1:nsamps
    %i
%     [T,E,Estar,LDH,perf,CPX] = adcx(tf_mol,tf_et,nr_t_mol,nr_t_et,T0,E0,g,r,kexp,gamma,...
%     CD20,CD16,RTX,kon20,koff20,kon16,koff16,gamma_perf);
    % generate one column of your sampled data by calling the adcx model
    % which calles the reaction_ss function
    [Ti, Ei, Estari, LDHi, perfi, CPXi] = adcx(pstruct.tf_mol,pstruct.tf_et,...
        pstruct.nr_t_mol,pstruct.nr_t_et,pstruct.T0,pstruct.E0toT0,...
        pstruct.Estar0,g,pstruct.r,pstruct.kexp,pstruct.gamma,...
        (CD20samps(i)),(CD16samps(i)),pstruct.RTX,...
        pstruct.kon20,pstruct.koff20,kon16,pstruct.koff16,...
        pstruct.h,pstruct.gamma_perf);  
	% whoever wrote this, you might want to consider callling adcx_wrapper
	% instead in the previous like, so that you don't have to change this code everytime we
	% introduce a new parameter... oh wait maybe that won't work because of
	% CD20samps and CD16samps...
    
    Tmat(:,i)     = Ti;
    Emat(:,i)     = Ei;
    Estarmat(:,i) = Estari;
    LDHmat(:,i)   = LDHi;
    perfmat(:,i)  = perfi;
    delCD16(i,1)  = CPXi;
   
 
    
end

% Make a matrix that is CD20 values * number of tumor cells at each time
CD20mat = zeros(size(Tmat));
    for i = 1:nsamps
    CD20mat(:,i) = Tmat(:,i).*CD20samps(i); 
    end

% Make the same matrix for CD16 using E and Estar mat


    CD16mat = zeros(size(Emat));
    for i = 1:nsamps
    CD16mat(:,i) = Emat(:,i).*CD16samps(i);
    CD16matEstar(:,i) = Estarmat(:,i).*(CD16samps(i)-delCD16(i));
    end


end

