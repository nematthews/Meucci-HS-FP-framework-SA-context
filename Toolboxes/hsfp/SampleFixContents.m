% Historical Simulation with Flexible Probabilities (HS-FP) Toolbox 
% Version 1.0 (R2022b) 07-March-2023 
%
% Author: 	     Nina Matthews
% Supervisor: 	 A/Prof. Tim Gebbie
% Project: 		 Time- and Market-Conditioned Flexible Probabilities: 
%		  		 in the Context of the South African Bond Market.
% Last edit: 	 17/05/2023
% Resources:     The MathWorks, Inc. (2019)  &  Attilio Meucci 
%			     (See Meucci reference list below)		
%
%   This program was developed to be utilised in my MSc Dissertation. 
%   The HSFP toolbox is used in portfolio optimisations for a variety of 
%   portfolios consisting of different combinations of asset classes. 
%   The functions come WITHOUT ANY WARRANTY.
%
% NO WARRANTY
%
% General information.
%   This toolbox houses funtions used to generate flexible probabilites
%   in the HSFP framework developed by Meucci (see references below). 
%
% General Functions:
%   xxxMainFNxx - xxxxxxx 
%   rw_probs    -  facilitates the calculation of the Time-conditioned 
%                      probabilities used in the HS-FP framework.
%                      It uses window lengths specified in number of mnths
%                      to generate the sequence of normalised Prs. 
%                      Default window length = 60 mnths (5 yrs).
%   xxxx    - 
%   xxxx    - 
%                  
%   xxx		     - 
%
% Test Code
%   hsfp_probs_1 - xxxxx
%
% The hrpestimate() function wraps all other sub function files needed for
% the HRP weights estimation into one implementation, executing all steps:
%
%
% Test Data
%   xxxxx
% 
% Obsolete functions.
%   None at this time.
%
% Others
%   None at this time.
%
% GUI Utilities.
%   None at this time.
%
% References:
% Meucci, A., 2008, Fully flexible views: Theory and practice, Risk 21, 
% 97—102 Article and code available at http://symmys.com/node/158.
%
% ______, 2010, Historical scenarios with Fully Flexible Probabilities. 
% http://symmys.com/node/150. Working Paper.
%
% ______, 2012, Effective number of scenarios with Fully Flexible Probabilities.
% http://symmys.com/node/362. GARP Risk Professional, 45-46.
%
% ______, 2013, Estimation and Stress-Testing via Time- and Market-Conditional
% Flexible Probabilities. http://symmys.com/node/600
%
% Copyright(c) 2023 Nina Matthews
%
%
% 1. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
% FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
% OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
% PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
% OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
% TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
% PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
% REPAIR OR CORRECTION.
%
% 2. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
% WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
% REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
% INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
% OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
% TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
% YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
% PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGES.
%








