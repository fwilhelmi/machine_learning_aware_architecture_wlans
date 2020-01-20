%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% SaveMatrixToFile.m --> 
%-------------------------------------------------------------------------

function [] = SaveMatrixToFile(m)

filename = 'output_stas.csv';

if isfile(filename)
    dlmwrite(filename, m, 'delimiter', ',', '-append');
else
    csvwrite(filename, m);
end

end
