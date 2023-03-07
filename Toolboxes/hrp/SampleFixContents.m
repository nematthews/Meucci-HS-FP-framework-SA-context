% Hierarchical Risk Parity (HRP) Toolbox 
% Version 1.0 (R2022b) 07-March-2023 
%
% Author: 	     Nina Matthews
% Supervisor: 	 A/Prof. Tim Gebbie
% Project: 		 Time- and Market-Conditioned Flexible Probabilities: 
%		  		 in the Context of the South African Bond Market.
% Last edit: 	 25/02/2023
% Resources:     The MathWorks, Inc. (2019)  &  Marcos Lopez De Prado (2016) 
%			     (See reference list below)		
%
%   This program was developed to be utilised in my MSc Dissertation. 
%   The HRP package is used for portfolio optimisation for a variety of 
%   portfolios consisting of different combinations of asset classes. 
%   The functions come WITHOUT ANY WARRANTY.
%
% NO WARRANTY
%
% General information.
%   This toolbox houses sub functions needed for each for the 3 Steps of 
%   De Prado's HRP algorithm as defined in his 2016 paper  Building 
%   Diversified Portfolios That Outperform Out-Of-Sample :
%
%	Step 1: Tree Clustering: Computing distance matrix & Building the 
%           clustering tree structure
%	Step 2: Quasi-Diagonalization: Creates quasidiagonalized covar matrix 
%           based on linkage matrix
%	Step 3: Recursive Bisection: Allocating the weights with bisection 
%           algorithm.
%
% The hrpestimate() function wraps all other sub function files needed for
% the HRP weights estimation into one implementation, executing all steps:
%
% The order of the files listed below correspond to their call order and 
% their nesting structure, see the individual function files for further 
% details on inputs and outputs of each.
% 1. hrpestimate(cov):
%   1.1 qdiag(link)
%        1.1.1 leafnodes(rootGroupNodeId, link)
%   1.2 hrpbisect(cov, sortedIdx)
%            1.2.1 ivp(cov(sortedIdx, sortedIdx))

% General Functions:
%   hrpestimate  - HRP weight estimations 
%   qdiag        - sorted asset index list for quasidiagonalized covar matrix 
%   leafnodes    - leaf nodes for a given group node id
%   hrpbisect    - recursively bisects assets and allocates weights using vip
%   ivp		   - performs simple inverse-variance allocation
% Test Code
%   hrp_weights_001 - generates the HRP weights for a portfolio of 10 assets
%					  given a dummy dataset developed by De Prado. 
%                     Results are consistent. 
%
% Test Data
%   Kipnis_covMat.csv
% 
% Obsolete functions.
%   None at this time.
%
% Others
%   None at this time
%
% GUI Utilities.
%   None at this time.

% Copyright(c) 2023 Nina Matthews
%
%
% 11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
% FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
% OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
% PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
% OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
% TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
% PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
% REPAIR OR CORRECTION.
%
% 12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
% WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
% REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
% INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
% OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
% TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
% YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
% PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGES.
%






