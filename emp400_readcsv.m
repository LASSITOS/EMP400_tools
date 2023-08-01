function em = emp400_readcsv(infile)
% A simple function to read abbreviated EMP 400 data from a comma-separated
% ASCII file written by emp400_write2csv

% USAGE:
%  em = emp400_readcsv(infile)
%
% INPUTS:
%  infile: fully-qualified filename for CSV file to be read
%
% OUTPUTS:
%      em: data structure containing abbreviated EMP400 data
%          - time: time of each observations as date string
%          -    t: time of each observations as 
%          -  lat: latitude of each observations (if GPS working)
%          -  lon: longitude of each observations (if GPS working)
%          -   Ix: in-phase component for frequency x
%          -   Qx: quadrature component for frequency x
%          -   Cx: apparent conductivity for frequency x
%          -  TTx: derived total snow+ice thickness for frequency x
%
% NOTES:
%   - lat and lon won't be present if GPS wasn't functioning during survey
%   - data for up to 3 frequencies may be present
%   - TT only present if computed (by emp400_Q2TT)

% Andy Mahoney      March 2023
% -----------------------------------------------------------------------


% Open specified output file for writing
infid = fopen(infile, 'rt');

% Get header fields
header = fgetl(infid);
headfields = regexp(header, ',', 'split');
Nf = numel(headfields);

% Construct format string
fmtstr = '';
for f=1:Nf
    if regexp(headfields{f}, 'time')
        fmtstr = [fmtstr '%s'];
    else
        fmtstr = [fmtstr '%f'];
    end
    
    if f < Nf
        fmtstr = [fmtstr ' '];
    else
        fmtstr = [fmtstr '\n'];
    end
end

% Read in data using text scan
indata = textscan(infid, fmtstr, 'delimiter', ',');

% Assign data columns to fields
for f = 1:Nf
    em.(headfields{f}) = indata{f};
end

em.t = datenum(em.time, 'yyyy-mm-ddTHH:MM:SS.fff');


end