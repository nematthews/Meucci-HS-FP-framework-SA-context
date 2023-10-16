% HSFPparameters class definition
classdef HSFPparameters
    properties
        RollWindow          (1,1) double {mustBePositive,mustBeInteger} = 12
        Tau                 (1,1) double    {mustBePositive} = 12
        Alpha               (1,1) double    {mustBePositive} = 0.1
        Z_target            (1,1) string {mustBeMember(Z_target,{'latest',...
                                'mean'})} = 'latest'
        Gamma               (1,1) double    {mustBePositive} = 2
        h                   (1,1) double     {mustBePositive} = 0.5
        Tau_prior           (1,1) double    {mustBePositive} = 12
        Ensemble_wt_method  (1,1) string {mustBeMember(Ensemble_wt_method,...
                             {'simple','log-linear'})} = 'log-linear'
    end
    
    methods
    %%%%% CONSTRUCTOR FUNCTION %%%%%%%%
    function HSFPparameters = HSFPparameters(RollWindow,Tau,Alpha,Z_target,Gamma,h,Tau_prior,Ensemble_wt_method)
            if  (nargin == 8)
                HSFPparameters.RollWindow = RollWindow;
                HSFPparameters.Tau = Tau;  
                HSFPparameters.Alpha = Alpha;
                HSFPparameters.Z_target = Z_target;
                HSFPparameters.Gamma = Gamma;
                HSFPparameters.h = h;
                HSFPparameters.Tau_prior = Tau_prior;
                HSFPparameters.Ensemble_wt_method = Ensemble_wt_method;
            elseif (nargin > 0)
                error("You need to provide all 8 required inputs (or leave empty for default values): " + ...
                    "RollWindow, Tau, Alpha, Z_target, Gamma, h, Tau_prior, Ensemble_wt_method ")
            end
    end 
    end
end