function emi = emp400_reademifile(infile, varargin)
% A simple function to read data files produced by the EMP-400 
% and return the following data structure containing EM data 
% and any GPS information 
%
% USAGE: 
%    emi  = readEMP400_emifile(infile)
%    emi  = readEMP400_emifile(___, 'logtoffset', logtoffset)
%    emi  = readEMP400_emifile(___, 'surveydate', surveydate)
%
% INPUTS
%       infile: fully-qualified filename for EMP-400 file
%   logtoffset: time offset in days to correct for logger timestamp error
%               Default: 0
%   surveydate: datestring specifying survey date in to correct for logger
%               timestamp error. Specify in same time zone as logger clock
%               Default: 'not specified'
%
% OUTPUTS
% emi: a data structure with the following fields
%               info: a data structure with numerous fields including 
%                     time, date, and configuration parameters            
%             record: a serial number for each measurement
%             xcoord: x-coord for a pre-assigned survey grid
%             ycoord: y-coord for a pre-assigned survey grid
%               time: time of each measurement as string
%                  t: time of each measurement as Matlab Julian day
%                     (corrected as necessary by toffset)
%                 I1: inphase component for 1st frequency
%                 Q1: quadrature component for 1st frequency
%                 C1: derived conductivity for 1st frequency
%                 I2: inphase component for 2nd frequency
%                 Q2: quadrature component for 2nd frequency
%                 C2: derived conductivity for 2nd frequency
%                 I3: inphase component for 3rd frequency
%                 Q3: quadrature component for 3rd frequency
%                 C3: derived conductivity for 2nd frequency
%             remark: string value of any remarks entered during survey
%               mark: 
%                lat: latitude for each measurement (if GPS enabled)
%                lon: longitude for each measurement (if GPS enabled)
%                alt: altitude for each measurement (if GPS enabled)
%               tilt: instrument tilt during each measurement
%             errors: error flags
%                gps: structure containing GPS NMEA info (if GPS enabled)
%                     NOTE: NMEA data only provides time of day and
%                           relies on logger time to calculate date
%         logtoffset: time offset correction specified in function call
%                     (see toffset input argument above)
%         surveydate: date string for survey date correction
%                     (see toffset input argument above)
%
%
% NOTES:
%     - if correcting problem with logger time stamps, user should use
%       EITHER the logtoffset or surveydate arguments and NOT BOTH
%
%       Andy Mahoney      March 2023
% -----------------------------------------------------------------------


% Check for optional input arguments
a = 1;
while (a < numel(varargin))
    switch varargin{a}
        case 'logtoffset'
            logtoffset = varargin{a+1};
            a = a + 2;
        case 'surveydate'
            surveydate = varargin{a+1};
            a = a + 2;
        otherwise
            disp('Argument not recognized:');
            disp(varargin{a});
            a = a + 1;
    end
end

% Apply default values as necessary
if ~exist('logtoffset', 'var')
    logtoffset = 0;
end
if ~exist('surveydate', 'var')
    surveydate = 'not specified';
end


% Check whether input file exists
if ~exist(infile, 'file')
    disp(['readEMP400_emifile: file not found: ', infile]);
    disp('Quitting.');
    return
end

infid = fopen(infile, 'rt');

% Assign string and float types to specific information fields
stringinfofields = {'gssi_profiler_ver', 'gssi_profiler_firmware_ver', ...
               'project_name', 'date', 'units', 'gps_status', ...
               'gps_validity_criteria', 'collection_mode', ...
               'grid_type', 'instrument_orientation', 'data_smoothing', ...
               'user_notes', 'data_file_recovery', ...
               'tlt_bit', 'bat_bit', 'adc_bit', 'tns_bit', 'nef_bit', 'rce_bit'};

floatinfofields = {'line_spacing', 'station_spacing', 'mark_spacing', ...
                'xmin', 'ymin', 'xmax', 'ymax', 'inphase_zero_levels', ...
                'quad_zero_levels', 'number_of_stacks', ...
                'calibration_height', 'power_line_frequency'};

            
% Read info from file header
% assign to variables depending on type of information
line = fgetl(infid);
while ~strcmp(line, '$$$')

    if ~any(line)
        line = fgetl(infid);
        continue;
    end
    
    linesplit = regexp(line, ',', 'split');
    infofield = regexprep(regexprep(linesplit{1},'\s','_'), '\W', '');
    infofield = lower(infofield);

    switch true
        case any(strcmp(infofield, stringinfofields))
            emi.info.(infofield) = linesplit{2};
        case any(strcmp(infofield, floatinfofields))
            emi.info.(infofield) = str2double(linesplit{2});
        case strcmp(infofield, 'frequencies')
            emi.info.(infofield) = str2double(linesplit(2:end-1));
        case strcmp(infofield, 'collection_mode')
            emi.info.(infofield) = linesplit{2};
            if ~strcmp(emi.info.collection_mode, 'Stationary')
                emi.info.sample_interval = str2double(linesplit{3});
                emi.info.sample_unit = linesplit{4};
            end
        case numel(regexp(infofield, '.*bit.*')) > 0
            infosplit = regexp(infofield, '__', 'split');
            emi.info.(infosplit{1}).val = infosplit{2};
            emi.info.(infosplit{1}).str = infosplit{3};
        otherwise
            disp(['Unrecognized header field: ' line]);
    end
    line = fgetl(infid);
end

emi.info.t0 = datenum(emi.info.date, 'mm/dd/yyyy HH:MM:SS.fff') ...
              + logtoffset;
if ~strcmp(surveydate, 'not specified')
    emi.info.t0 = datenum(surveydate) + mod(emi.info.t0, 1);
end

% Count number of frequencies measured
Nfreq = numel(emi.info.frequencies);

% Assign specific variables to be string variables
stringheadfields = {'time', 'remark', 'mark', 'tilt', 'errors'};

% Read data header for field names
headline = fgetl(infid);

% Rename frequencies by numbers (e.g. f1, f2, f3)
for f=1:Nfreq
    freqstr = num2str(emi.info.frequencies(f),'%0.0f');
    fstr = num2str(f, '_f%d');
    headline = regexprep(headline, ['\[' freqstr '\]'], fstr);
    notalot = 0;
end

% Extract data field names and convert to lower case
datafields = lower(regexp(headline, ',', 'split'));
datafields = regexprep(datafields, '[\W ]', '');

Nfields = numel(datafields);

% Concatenate format string for reading data
infstr = '';
for f=1:Nfields
    if any(strcmp(stringheadfields, datafields{f}))
        infstr = [infstr '%s '];
    else
        infstr = [infstr '%f '];
    end
end
infstr = infstr(1:end-1);

% Read data
indata = textscan(infid, infstr, 'delimiter', ',');

% Close file
fclose(infid);

% Assign to data structure
for f = 1:Nfields
    emi.(datafields{f}) = indata{f};
end

% Convert timestamps into date numbers
emi.t = floor(emi.info.t0) + mod(datenum(emi.time, 'HH:MM:SS.fff'),1);
   
                       
% Do a little QC
if isfield(emi, 'lat') && isfield(emi, 'long')
    badlatlon = (emi.lat < -90) | (emi.lat > 90) | ...
                (emi.long < -180) | (emi.long > 180);
    emi.lat(badlatlon) = NaN;
    emi.long(badlatlon) = NaN;
end                       
   
% Check to see if associated GPS file exists
[inpath, infstem, inext] = fileparts(infile);
gpsfile = [inpath filesep infstem '.GPS'];
gpsfmt = '%s %s %f %s %f %s %f %f %f %f %s %f %s %s %s %s';
if exist(gpsfile, 'file')
    % Read NMEA data lines from .GPS file
    infid2 = fopen(gpsfile, 'rt');
    nmea = textscan(infid2, gpsfmt, 'delimiter', ',');
    fclose(infid2);
    
    % Check to see if valid NMEA data was found
    nmeavalid = ~any(cellfun('isempty', nmea));
    if nmeavalid
    
        % Convert NMEA data to useful numerical values
        emi.gps.tod = mod(datenum(nmea{2}, 'HHMMSS'),1);       % Time of day

        latdeg = floor(nmea{3}/100);                           % Latitude (ddmmm.mm)
        latmin = mod(nmea{3}, 100);
        lathem = strcmp(nmea{4}, 'N') - strcmp(nmea{4}, 'S');
        emi.gps.lat = lathem.*(latdeg + latmin/60);
        londeg = floor(nmea{5}/100);                           % Longitude (ddmmm.mm)
        lonmin = mod(nmea{5}, 100);
        lonhem = strcmp(nmea{6}, 'E') - strcmp(nmea{6}, 'W');
        emi.gps.lon = lonhem.*(londeg + lonmin/60);

        emi.gps.fix = nmea{7};                                 % GPS fix quality (0 = none; 1 = GPS; 2 = DGPS)
        emi.gps.nsat = nmea{8};                                % # satellites
        emi.gps.hdop = nmea{9};                                % HDOP
        emi.gps.alt = nmea{10};                                % GPS altitude
        emi.gps.altunit = nmea{11};
        emi.gps.geoidht = nmea{12};                            % Geoid height
        emi.gps.geoidhtunit = nmea{13};
        emi.gps.timesinceDGPS = nmea{14};
        emi.gps.checksum = nmea{15};                            % Check sum
        emi.gps.logtime = datenum(nmea{16}, 'mm/dd/yy HH:MM:SS.fff');  % Logger time stamp

        % Correct for time zone or clock offset
        t_off = emi.gps.tod - mod(emi.gps.logtime,1);
        t_off_mean = mean(t_off);
        emi.gps.meanclockoffset = t_off_mean + 1*(t_off_mean < -0.5) ...
                                             - 1*(t_off_mean >= 0.5);
        emi.t = emi.t + emi.gps.meanclockoffset; 

        notalot = 0;
    end    
end

% Rename some length fields to make them shorter and more usable
fieldorder = fieldnames(emi);
fields2rename = {'inphase_f1', 'quad_f1', 'conductivity_f1', ...
                 'inphase_f2', 'quad_f2', 'conductivity_f2', ...
                 'inphase_f3', 'quad_f3', 'conductivity_f3', ...
                 'long'};
newfieldnames = {'I1','Q1','C1','I2','Q2','C2','I3','Q3','C3','lon'};
for f=1:numel(fields2rename)
    if isfield(emi, fields2rename{f})
        emi.(newfieldnames{f}) = emi.(fields2rename{f});
        emi = rmfield(emi, fields2rename{f});
        fieldi = ~cellfun(@isempty, regexp(fieldorder, fields2rename{f}));
        fieldorder{fieldi} = newfieldnames{f};
    end
end
emi = orderfields(emi, fieldorder);


end

