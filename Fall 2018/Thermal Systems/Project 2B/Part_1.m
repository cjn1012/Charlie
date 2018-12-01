% A study of the effects of compressor, turbine and nozzle efficiencity to the overall cycle

% Sum of efficiencies cannot exceed 260%
    % Nozzle can not exceed 98%
    % Turbine can not exceed 90%
    % Compressor can not exceed 90%
    
% Choose the best percentages and graph a T-s diagram

clear all
close all

CompressorEff = linspace(.70,.90,21);
TurbineEff = linspace(.78,.98,21);
NozzleEff = linspace(.70,.90,21);

EffMatrix = zeros(200,4);
l = 1;

for i = CompressorEff
    for j = TurbineEff
        for k = NozzleEff
            if i+j+k == 2.60
                EffMatrix(l,1) = i;
                EffMatrix(l,2) = j;
                EffMatrix(l,3) = k;
                EffMatrix(l,4) = Optimization(i,j,k);
                l = l + 1;
            end
        end
    end
end

                
Sorted = sort(EffMatrix,4);                

Eff_Comp = Sorted(1,1);
Eff_Turb = Sorted(1,2);
Eff_Nozz = Sorted(1,3);

[EntropyIdeal,TemperatureIdeal,EntropyActual,TemperatureActual] = T_s_Data(Eff_Comp,Eff_Turb,Eff_Nozz);

plot(EntropyIdeal,TemperatureIdeal,'b')
hold on
plot(EntropyActual,TemperatureActual,'r')






































