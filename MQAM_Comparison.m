%--------------------------------------------------------------------------
% MQAM_Comparison.m -- ETSIT-UPM LTDS 2010-2011
%
% Analisys of several QAM modulation schemas returning the Bit Error Rate
% from each one for a given SNR.
%
% Authors:
%   Luis Antonio �beda Medina (lubeme23@gmail.com)
%   H�ctor Veiga Ortiz (hveiga@hawk.iit.edu)
%
% Input:
%  SNR - Signal Noise Ratio for the AWGN Channel in dB. 
%
%  Example: MQAM_Comparison(5)
%
%
% Copyright 2010 H�ctor Veiga Ortiz and Luis Antonio �beda Medina
% 
%
% This file is part of MQAM Simulator.
% 
% MQAM Simulator is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MQAM Simulator is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MQAM Simulator.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function MQAM_Comparison(SNR)

% Modulation deepness. Analysing until 4096-QAM.
% It would be possible to analyze bigger modulations adding
% integers to depth.
depth = [1 2 3 4 5 6];
% Amount of random bits
longitud = 2000; 

BER=[];
% Analysing BER for each modulation. Showing results on screen.
disp('---------- M-QAM MODULATION COMPARISON ---------------------');
for n=1:length(depth)
    nivel = 4^depth(n);
    BER(n) = MqamMed(longitud,depth(n),SNR);
    disp(['Modulation: ' num2str(nivel) '-QAM BER = ' num2str(BER(n))]);
end
disp('------------------------------------------------------------');



end


