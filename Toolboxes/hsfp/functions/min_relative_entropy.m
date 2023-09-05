%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_post, DLagrang] = min_relative_entropy(prior, inequalityMatrix, inequalityVector, equalityMatrix, equalityVector, options)
% Calculates posterior probabilites from a supplied prior using the dual
% Lagrangian multiplier and the Hessian. 
%
%% INPUTS: 
% prior - desired prior probability distribution (usually exp smoothed Prs)
% (Type: array double, [1 x T])
%
% inequalityMatrix - matrix defining inequality views 
% Here for z^2
% (Type: matrix double, [1 x T])
%
% inequalityVector - vector for the defined inequality views 
% Here: (cr_mu^2)+cr_variance
% (Type: double) 
%
% equalityMatrix - equality constrains
% Here: % Pr sum to 1
% (Type: matrix, [number of EConstraints x T])
%
% equalityVector - equality views 
% Here: Constrain the 1st moments
% (Type: array double, [number of EConstraints x 1])
%
% options : options that can be specified for optimization
% Note: specify using "optimoptions"


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $


numVariables = size(inequalityMatrix, 1);
numEqualityConstraints = size(equalityMatrix, 1);

if nargin < 6 || isempty(options)
    if isempty(inequalityMatrix)
        options = optimset('GradObj', 'on', 'Hessian', 'on', ...
            'MaxIter', 10^6, 'TolFun', 1e-16, 'MaxFunEvals', ...
            10^6, 'Display', 'off');
    else
        options = optimset('GradObj', 'on', 'Hessian', ...
            'user-Supplied', ...
            'HessFcn', @(lv, lambda) nestedmHessian(lv, lambda, ...
            prior, inequalityMatrix, equalityMatrix, numVariables, numEqualityConstraints), ...
            'MaxIter', 10^6, 'TolFun', 1e-16, 'MaxFunEvals', ...
            10^6, 'Display', 'off');
    end
end

%% Optimization initialisation:

lv0 = zeros(numVariables + numEqualityConstraints, 1);

% Apply Lagrangian based on type of constraints:
% 1. Lagrangian multiplier of view on mu, specify constraints: l_ != 0
% 2. Lagrangian multiplier of view on sigma may = 0

% If EQALITY:
if ~numVariables
    v_ = fminunc(@(v) nestedmDLagrang_eq(v, prior, equalityMatrix, equalityVector, numEqualityConstraints), ...%
        lv0, options);
    p_post = exp(log(prior) - 1 - v_' * equalityMatrix);
    DLagrang = v_;
    
    % If EQUALITY & INEQUALITY
else
    alpha = -eye(numVariables + numEqualityConstraints);
    alpha(numVariables + 1:end, :) = [];
    beta = zeros(numVariables, 1);
    DLagrang = fmincon(@(lv) nestedmDLagrang(lv, prior, inequalityMatrix, inequalityVector, equalityMatrix , equalityVector, numVariables), ...
        lv0, alpha, beta, [], [], [], [], [], options);
    numEqualityConstraints = DLagrang(1:numVariables);
    v_ = DLagrang(numVariables + 1:end);
    p_post = exp(log(prior) - 1 - numEqualityConstraints' * inequalityMatrix - v_' * equalityMatrix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested function 1: Dual Langrang for equality constraints
% Dual Langrangian and Hessina with aeq & beq
    function [mh, mgrad, mHess] = nestedmDLagrang_eq(v, prior, equalityMatrix, equalityVector, numEqualityConstraints)

        p = exp(log(prior) - 1 - v' * equalityMatrix);
        p = max(p, 10^(-32));
        p_trans = p';
        h = (log(p) - log(prior)) * p_trans + v' * (equalityMatrix * p_trans - equalityVector);

        mh = - h; %min
        mgrad = [];
        mHess = [];

        if nargout > 1
            mgrad = equalityVector - equalityMatrix * p_trans; 
        end

        if nargout > 2
            % h = a * diag(p) * a'
            mHess = equalityMatrix * diag(p) * equalityMatrix'; 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested function 2: Dual Langrang for equality & inequality constraints
% Dual Langrangian and Hessina with a & b
    function [mh, mgrad] = nestedmDLagrang(lv, prior, inequalityMatrix, inequalityVector, equalityMatrix , equalityVector, numVariables)
        l = lv(1 : numVariables);
        v = lv(numVariables + 1 : end);

        p = exp(log(prior) - 1 - l' * inequalityMatrix - v' * equalityMatrix);
        p = max(p, 10^(-32));
        p_trans = p';

        %  DLagrang
        h = (log(p) - log(prior)) * p_trans + l' * (inequalityMatrix * p_trans - inequalityVector) + v' * (equalityMatrix * p_trans - equalityVector);
        mh = - h; % min

        if nargout > 1
            mgrad = [inequalityVector; equalityVector] - [inequalityMatrix; equalityMatrix] * p_trans; 
        end
    end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function mHess = nestedmHessian(lv, lambda, prior, inequalityMatrix, equalityMatrix, numVariables, numEqualityConstraints)
        % Hessian of opposite dual Lagrangian called by
        % calculations for dual Lagrangian for inequality & equality constraints
        l = lv(1 : numVariables);
        v = lv(numVariables + 1 : end);

        p = exp(log(prior) - 1 - l' * inequalityMatrix - v' * equalityMatrix);
        p = max(p, 10^(-32));
        % h = [c; a] * diag(p) * [c; a]'
        mHess = ([inequalityMatrix; equalityMatrix] .* (ones(numVariables + numEqualityConstraints, 1) * p )) * [inequalityMatrix', equalityMatrix'];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


