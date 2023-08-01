function emp400_emi2csv(varargin)
% A simple function to read data files produced by the EMP-400 
% and return the following data structure containing EM data 
% and any GPS information 

% USAGE: emp400_emi2csv
%        emp400_emi2csv('infiles', infiles);
%        emp400_emi2csv(___, 'outpref', outpref);
%        emp400_emi2csv(___, 'outsuff', outsuff);
%        emp400_emi2csv(___, 'fnamedate', fnamedate);
%        emp400_emi2csv(___, 'outpath', outpath);
%        emp400_emi2csv(___, 'surveydates', surveydates);
%        emp400_emi2csv(___, 'logtoffsets', logtoffsets);


% Inputs:
%     infiles: list of filename for EMP-400 files (.EMI)
%              Default: if no files specified, user is prompted to select
%     outpref: prefix to use at beginning of output filename
%              Default: input filename excluding extension
%     outsuff: sufffix to use at end of output filename
%              Default: empty string
%   fnamedate: if set to "on" adds date-time string to filename based on
%              timestamp of first record on EMI file
%              default: "off"
%     outpath: Path for folder to use for output files. Folder is created
%              if it doesn't already exist
%              Default: same as input file
% logtoffsets: array of time offsets (in days) to correct for error in
%              logger timestamps. There should either be one offset for
%              each input file, or a single offset to be applied to all
%              files
% surveydates: cell array of date strings specifying survey dates  to 
%              correct for errors in logger timestamps. There should either
%              be one date for each input file, or a single date to be 
%              applied to all files


% Outputs
% An comma-separated ASCII output file is created for each input file
% with a filename of the format:
%    [outpath]/[outpref]_[fnamedate][outsuff].csv


% Check for optional arguments
a = 1;
while a < numel(varargin)
    switch varargin{a}
        case 'infiles'
            infiles = varargin{a+1};
            a = a + 2;
        case 'outpref'
            outpref = varargin{a+1};
            fixedoutpref = 1;
            a = a + 2;
        case 'outsuff'
            outsuff = varargin{a+1};
            a = a + 2;
        case 'fnamedate'
            fnamedate = varargin{a+1};
            a = a + 2;
        case 'outpath'
            outpath = varargin{a+1};
            a = a + 2;
        case 'logtoffsets'
            logtoffsets = varargin{a+1};
            a = a + 2;
        case 'surveydates'
            surveydates = varargin{a+1};
            a = a + 2;
        otherwise
            disp('Argument not recognized: ');
            disp(varargin{a});
            a = a + 1;
    end
end

% Check to see if input file list was specified in function call
if ~exist('infiles', 'var')
    [flist, inpath] = uigetfile('~/*.EMI', 'Select EMI files', ...
                                  'multiselect', 'on');
    
    % If only one file selected, convert to 1-element cell array
    if ~iscell(flist)
        temp = flist;
        flist = {temp};
        clear('temp');
    end                          
                              
    % Add path to each file name
    Nf = numel(flist);
    infiles = cell(Nf,1);
    for f=1:Nf
        infiles{f} = [inpath filesep flist{f}];
    end
end

% If input file not specified as a cell array, convert to 1-element cell array
if ~iscell(infiles)
    temp = infiles;
    infiles = {temp};
    clear('temp');
end
Nf = numel(infiles);

% Set fnamedate to "off" if not specified in function call
if ~exist('fnamedate', 'var')
    fnamedate = 'off';
end

% Set default values for optional time-related arguments
if ~exist('logtoffsets', 'var')
    logtoffsets = 0;
end
if numel(logtoffsets) == 1, logtoffsets = repmat(logtoffsets, Nf,1); end

if ~exist('surveydates', 'var')
    surveydates = {'not specified'};
end
if numel(surveydates) == 1, surveydates = repmat(surveydates, Nf,1); end

if ~exist('fixedoutpref', 'var')
    fixedoutpref = 0;
end

% Go through each file in file list
for f=1:Nf
    
    % Read EMI file
    emi = emp400_reademifile(infiles{f}, 'logtoffset', logtoffsets(f), ...
                                           'surveydate', surveydates{f});

    % Extract output prefix from filename if not specified in function call
    if ~fixedoutpref
        [~,outpref,~] = fileparts(infiles{f});
    end
    
    % Assign date stamp to filename if specified in function call
    if strcmp(fnamedate, 'on')
        fdstr = ['_' datestr(emi.t(1), 'yyyymmdd_HHMMSS')];
    else
        fdstr = '';
    end
    
    % Assign blank string if output suffix not specified in function call
    if ~exist('outsuff', 'var')
        outsuff = '';
    end
    
    % Extract output path from filename if not specified in function call
    if ~exist('outpath', 'var')
        [outpath,~,~] = fileparts(infiles{f});
    end
    
    % Make output folder if it doesn't exist
    if ~exist(outpath, 'file')
        mkdir(outpath);
    end
    
    % Construct output filename including path
    outfile = [outpath filesep outpref fdstr outsuff '.csv'];

    % Write EMI data to csv file
    emp400_write2csv(emi, outfile);
    
    disp(['Written: ' outpref fdstr outsuff '.csv']);
end

end
    
    