%--------------------------------------------------------------------------
% MqamMed.m -- ETSIT-UPM LTDS 2010-2011
%
% Support function used by M-QAM_Election.m and M-QAM_Comparison.m
%
% Authors:
%   Luis Antonio Úbeda Medina (lubeme23@gmail.com)
%   Héctor Veiga Ortiz (hveiga@hawk.iit.edu)
%
% Input:
%  longitud - Amount of random bits to be tested. 
%  M - Deepness of modulation: 1 = 4-QAM, 2 = 16-QAM, 3 = 64-QAM.
%  SNR - SNR - Signal Noise Ratio for the AWGN Channel in dB.
%
%  Example: MqamMed(1000,2,5)
%
%
% Copyright 2010 Héctor Veiga Ortiz and Luis Antonio Úbeda Medina
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

function [BER] = MqamMed(longitud,M,SNR)
close all;
resto = rem(longitud,2*M);
bin = randint(longitud+(2*M-resto),1,2);
L=length(bin); 

k=100;
f=10;
bit1=ones(1,k);
bit0=0*bit1;
Nivel=M;
symbol=ones(1,2*Nivel*k);
mbit=[];mx=[];my=[];
bin = fliplr(bin);

x=0;
y=0;

for n=0:2*Nivel:L-2*Nivel
    bit=[];
    xi=0;
    yi=0;

    for m= 1:2:2*Nivel
        if bin(n+m)==0      
            xi=xi+(2^((m-1)/2));
            bit=[bit bit0];
        else
            xi=xi-(2^((m-1)/2));
            bit=[bit bit1];
        end
        if bin(n+m+1)==0
            yi=yi+(2^((m-1)/2));
            bit=[bit bit0];
        else
            yi=yi-(2^((m-1)/2));
            bit=[bit bit1];
        end
    end
    
    x=xi*symbol;
    y=yi*symbol;

    mx=[mx x];
    my=[my y];
    
    mbit=[mbit bit];
end

v=0:2*pi/k:2*pi*L-2*pi/k;
msync = mx + my*j;
qam=real(msync).*cos(f*v)-imag(msync).*sin(f*v); 

Vn=awgn(qam,SNR,'measured');

Vnx=Vn.*cos(f*v);
Vny=Vn.*-1.*sin(f*v);

[b,a]=butter(2,0.04);
Hx=2.*filter(b,a,Vnx);
Hy=2.*filter(b,a,Vny);
ML=length(Hx);

mdeb=[];
for m=k:2*Nivel*k:ML
    sym=[];
    thx=0;thy=0;
    fx=0;fy=0;
    for n=1:Nivel
        if Hy(m) > thy
            sym = [sym bit0];
 
            fy=1;
        else
            sym = [sym bit1];
            fy=-1;
        end

        if Hx(m) > thx
            sym = [sym bit0];
            fx=1;
        else
            sym = [sym bit1];
            fx=-1;
        end
        thy = thy + fy*(2^(Nivel-n));
        thx = thx + fx*(2^(Nivel-n));
    end
   
    mdeb=[mdeb fliplr(sym)];
end

err = 0;
bin2 = [];

for n=1:length(bin)
    if(bin(n) == 1) 
        bin2 = [bin2 bit1];
    elseif bin(n) == 0
        bin2 = [bin2 bit0];
    end
end

for n=1:length(bin2)
    if(bin2(n) ~= mdeb(n)) 
        err = err + 1;
    end
end

error_bits_total = round(err/k);

ber = error_bits_total/length(bin);
BER = ber;