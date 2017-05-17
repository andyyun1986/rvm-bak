% SB2_INITIALISATION  Initialise everything for the SPARSEBAYES algorithm
%
% [LIKELIHOOD, BASIS, BASISSCALES, ALPHA, BETA, MU, PHI, USED] = ...
%    SB2_INITIALISATION(TYPE, BASIS, TARGETS, SETTINGS, OPTIONS)
%
% OUTPUT ARGUMENTS:
% 
%	LIKELIHOOD	Likelihood structure (from SB2_LIKELIHOODS)
%	BASIS		Pre-processed full basis matrix
%	BASISSCALES	Scaling factors from full basis pre-processing
%	ALPHA		Initial hyperparameter alpha values
%	BETA		Initial noise level (Gaussian)
%	MU			Initial weight values
%	PHI			Initial "relevant" basis
%	USED		Indices of initial basis vectors (columns of BASIS)
% 
% INPUT ARGUMENTS:
% 
%	TYPE		Likelihood: one of 'Gaussian', 'Bernoulli', 'Poisson'
%
%	BASIS		NxM matrix of all possible basis vectors 
%				(one column per basis function)
%
%	TARGETS		N-vector with target output values
% 
%	SETTINGS	Initialisation structure for main parameter values via
%				SB2_PARAMETERSETTINGS 
% 
%	OPTIONS		User options structure from SB2_USEROPTIONS
%	
% NOTES: 
% 
% This function initialises all necessary model parameters appropriately
% before entering the main loop of the SPARSEBAYES inference algorithm.
%
% This function is intended for internal use by SPARSEBAYES only.
%

%
% Copyright 2009, Vector Anomaly Ltd
%
% This file is part of the SPARSEBAYES library for Matlab (V2.0).
%
% SPARSEBAYES is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% SPARSEBAYES is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with SPARSEBAYES in the accompanying file "licence.txt"; if not, write to
% the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
% MA 02110-1301 USA
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function [BASIS, BasisScales, Alpha, beta, Mu, PHI, Used] = ...
    SB2_Initialisation(BASIS, Targets)

% A "reasonable" initial value for the noise in the Gaussian case
GAUSSIAN_SNR_INIT	= 0.1;

% "Reasonable" initial alpha bounds
INIT_ALPHA_MAX	= 1e3;

% 
% BASIS PREPROCESSING: 
% 
% Scale basis vectors to unit norm. This eases some calculations and 
% will improve numerical robustness later.
% 
[BASIS, BasisScales] = SB2_PreProcessBasis(BASIS);

%% Noise model considerations
%
% In the Gaussian case, initialise "sensibly" to a proportion of the
% signal level (e.g. 10 percent)
% 

  % Noise initialisation

    % catch the pathological case where all outputs are zero
	% (although we're probably doomed anyway if that's true)
    stdt	= max([1e-6 std(Targets)]);
    
    % Initialise automatically approximately according to "SNR"
    beta	= 1/(stdt*GAUSSIAN_SNR_INIT)^2;



% Initialise basis (PHI), mu and alpha
%
% Either as specified by the SETTINGS structure, or heuristically
% 

% First, compute 'linearised' output for use in heuristic initialisation
%
TargetsPseudoLinear	= Targets; % standard linear case

% 1) the starting basis, PHI

%
% At this point Used will contain both pre-specified relevant basis
% functions (SB2_PARAMETERSETTINGS) and any "free" ones (SB2_USEROPTIONS).
% If neither are active, Used will be empty.
% 

  % Set initial basis to be the largest projection with the targets
  proj			= (BASIS'*TargetsPseudoLinear);
  [~, Used]	= max(abs(proj));
 
PHI	= BASIS(:,Used);

% 2) the most probable weights

Mu		= [];
 
% 3) the hyperparameters


    % Exact for single basis function case (diag irrelevant), 
    % heuristic in the multiple case
    p		= diag(PHI'*PHI)*beta;
    q		= (PHI'*Targets)*beta;
    Alpha	= p.^2./(q.^2-p); %
	% The main algorithm will handle these automatically shortly
	% (i.e. prune them)
    Alpha(Alpha<0) = INIT_ALPHA_MAX; 
    %
