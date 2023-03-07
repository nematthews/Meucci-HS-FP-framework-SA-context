% Hierarchical Risk Parity (HRP) Toolbox 
% Version 1.0 (R2022b) 07-March-2023 
%
% Author:  	     Nina Matthews
% Supervisor: 	 A/Prof. Tim Gebbie
% Project: 		 Time- and Market-Conditioned Flexible Probabilities: 
%		  		 in the Context of the South African Bond Market.
% Last edit: 	 25/02/2023
% Resources:     The MathWorks, Inc. (2019)  &  Marcos Lopez De Prado (2016) 
%				 (See reference list below)		
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
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
% It will call on sub function files within /Source when executing the 
% function. The Fn_HRP_estimate() function wraps all other sub function 
% files into one implementation, executing all steps:
%
% The order of the files listed below correspond to their call order and 
% their nesting structure, see the individual function files for further 
% details on inputs and outputs of each. 
% 1. Fn_HRP_estimate(assetCovar):
%   1.1 quasiDiagSort(link)
%        1.1.1 getLeafNodesInGroup(rootGroupNodeId, link)
%   1.2 helperBisectHRP(assetCovar, sortedIdx)
%            1.2.1 allocByInverseVariance(assetCovar(sortedIdx, sortedIdx))

% General Functions:
%   besaaip      - Allinprice for South African bonds using BESA specification
%   besainfaip   - Allinprice for inflation linked bonds
%   besatenor    - Bond tenor
%   besaconv     - Bond Convexity
%   besamoddur   - Bond Modified Duration
%   besaimpytm   - Implied yield-to-maturity 
%   besafwdprc   - Bond forward price
%   besaoption   - Bond option using Blacks model as per BESA specification
%   ncdproceeds  - NCD proceeds
%   ncdmvalue    - NCD maturity value
%   eqvalue      - Equivalent value of continuous rate
%   settledate   - Settlement date rule
%   valuedate    - Valuation date
%   cpiratio     - Compute the CPI ratio for inflation link bonds
%   cpirefdates  - Compute the reference dates for CPI data
%   holidays     - Holidays and non-trading days in South Africa and the US
%       
% Test Code
%   besa_gilts_001     - test allinprice computations
%   besa_options_002   - test options pricing 
%
% Test Data
%   None at this time.
% 
% Obsolete functions.
%   None at this time.
%
% Others
%   None at this time
%
% GUI Utilities.
%   None at this time.

% Copyright(c) 2004 Tim Gebbie
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%
% NO WARRANTY
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

% $Revision: 1.4 $ $Date: 2006-03-30 15:39:31+02 $



