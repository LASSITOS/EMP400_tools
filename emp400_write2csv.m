function emp400_write2csv(emi, outfile)
% A function to take an EMP-400 Matlab data structure and
% write the principal data fields to a comma-separate ASCII file
% Specifically, the following fields are output:
% - time / date
% - lat
% - lon
% - inphase component for each frequency
% - quadrature component for each frequency
% - conductivity for each frequency
% - total snow+ice thickness (TT) for each frequency (if present in emi)
%
% USAGE: emp400_write2csv(emi, outfile)
%
% INPUTS:
%      emi: EMP-400 data structure
%  outfile: fully-qualified filename for CSV file to be written
%
%
% Andy Mahoney      March 2023
% -----------------------------------------------------------------------

% Open specified output file for writing
outfid = fopen(outfile, 'wt');


% Check what fields are present in data structure and prepare
% header line and format string accordingly
headline = 'time';
outfmt = {'%19s'};
outfields = {'t'};

if isfield(emi, 'lat') && isfield(emi, 'lon')
    headline = [headline ',lat,lon'];
    outfmt = [outfmt, '%10.6f', '%11.6f'];
    outfields = [outfields, 'lat', 'lon'];
end

for f=1:3
    fn = num2str(f,'%d');
    if isfield(emi, ['I' fn])
        headline = [headline ',' ['I' fn]];
        outfmt = [outfmt, '%7.1f'];
        outfields = [outfields, ['I' fn]];
    end
    if isfield(emi, ['Q' fn])
        headline = [headline ',' ['Q' fn]];
        outfmt = [outfmt, '%7.1f'];
        outfields = [outfields, ['Q' fn]];
    end
    if isfield(emi, ['C' fn])
        headline = [headline ',' ['C' fn]];
        outfmt = [outfmt, '%10.4f'];
        outfields = [outfields, ['C' fn]];
    end
    if isfield(emi, ['TT' fn])
        headline = [headline ',' ['TT' fn]];
        outfmt = [outfmt, '%5.2f'];
        outfields = [outfields, ['TT' fn]];
    end
end

% Write the header line
fprintf(outfid, '%s\n', headline);
 
% Write each line of data to file, field-by-field    
Np = numel(emi.t);
outdatefmt = 'yyyy-mm-ddTHH:MM:SS.fff';
for p=1:Np
    fprintf(outfid, outfmt{1}, datestr(emi.t(p), outdatefmt));
    for f = 2:numel(outfields)
        fprintf(outfid, [',' outfmt{f}], emi.(outfields{f})(p));
    end
    fprintf(outfid, '\n');
end
        
end
