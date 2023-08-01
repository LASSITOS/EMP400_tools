function cal = emp400_extractcalpoints(em,outfile,calspd,calseglen)
% Function to identify observations along one or more EMP-400 survey tracks 
% where the instrument was stationary for an extended period of time. 

% These observations can be used to help match EMP400 data with 
% drill-hole data for calibration

% Andy Mahoney - 2023

% USAGE:
%     cal = emp400_extractcalpoints(em)
%     cal = emp400_extractcalpoints(em, outfile)
%     cal = emp400_extractcalpoints(em, outfile, calspd)
%     cal = emp400_extractcalpoints(em, outfile, calspd, calseglen)

% INPUT: 
%        em: data structure containing abbreviated EMP400 data 
%            or an array of such structures
%            (see emp400_write2csv.m or emp400_readcsv.m)
%   outfile: name of output file for saving calibration data points
%            in a comma-separated ASCII file
%    calspd: cut-off speed for identifiying stationary moments
%            default: 0.5 m/s
% calseglen: minimum duration of stationary moment
%            default: 10 s


% OUTPUT:
%      cal: data structure containing the following fields
%         - time: time at which EMP400 was stopped
%         -  lat: mean latitude for stopped period
%         -  lon: mean longitude for stopped period
%         -  spd: mean speed during stopped period
%         -   I1: mean inphase component during stopped period
%                 for frequency 1 
%         -   Q1: mean quadrature component during stopped period
%                 for frequency 1 
%         -   Q1: mean conductivity during stopped period
%                 for frequency 1
%         -   I2: mean inphase component during stopped period
%                 for frequency 2 
%         -   Q2: mean quadrature component during stopped period
%                 for frequency 2 
%         -   Q2: mean conductivity during stopped period
%                 for frequency 2
%         -   I3: mean inphase component during stopped period
%                 for frequency 3 
%         -   Q3: mean quadrature component during stopped period
%                 for frequency 3 
%         -   Q3: mean conductivity during stopped period
%                 for frequency 3


% Prompt user for input files if no emi data structures specified
if ~exist('em', 'var')
    [emfile, empath] = uigetfile('~/*.csv');
    
    if ~iscell(emfile)
        xxx = emfile;
        emfile = cell(1);
        emfile{1} = xxx;
        clear('xxx');
    end
    
    Nem = numel(emfile);
    for f = 1:Nem
        emfile{f} = [empath filesep emfile{f}];
    end
    
    em = [];
    for f=1:Nem
        em_in = emp400_reademifile(emfile{f});
        em = [em; em_in];
    end
    
end

  
% Number of EMP400 surveys
Nem = numel(em);

% Set default calspd if none specified
if ~exist('calspd','var')
    calspd = 0.5; % m/s
end

% Set default calseglen if none specified
if ~exist('calseglen','var')
    calseglen = 10; % seconds
end

% Define reference Ellispoid for distance/speed calculations
ellipsoid = referenceEllipsoid('wgs84', 'm');

% Go through each input file
for i = 1:Nem 
    
    % Exclude bad GPS data
    badi = isnan(em.lat) | isnan(em.lon);
    if any(badi)
        fn = fieldnames(em);
        for f=1:numel(fn)
            em.(fn{f}) = em.(fn{f})(~badi);
        end
    end
    
    % Calculate time interval between each EM measurment
    dt = (em(i).t(2:end)-em(i).t(1:end-1))*86400;

    % Calculate running mean lat/lon over calseglen period
    kw = round(calseglen/mean(dt));
    smkern = ones(kw,1)/kw;
    %lat_sm = nanconv(em(i).lat, smkern, 'same');
    %lon_sm = nanconv(em(i).lon, smkern, 'same');
    lat_sm = smooth1_noedge(em(i).lat, kw);
    lon_sm = smooth1_noedge(em(i).lon, kw);

    % Calculate distance moved between smoothed EM positions
    d = distance(lat_sm(1:end-1), lon_sm(1:end-1), ...
                 lat_sm(2:end), lon_sm(2:end), ellipsoid);

    % Set beginning and end margins to NaN
    % (to exclude edge effects from convolution)
    d(1:ceil(kw/2)) = NaN;
    d(end-ceil(kw/2):end) = NaN;

    % Calulate resulting speed         
    spd = d./dt;

    % Find occasions when EM speed was below threshold
    lowspd = spd < calspd;
    em_stop = diff(lowspd) > 0;
    em_go = diff(lowspd) < 0;
    seg0 = find([lowspd(1); em_stop] == 1);
    seg1 = find([em_go; lowspd(end)] == 1);
    if numel(seg0) ~= numel(seg1)
        wtf = 1;
    end

    % Identify stops that lasted longer than minimum time
    seglen = (em(i).t(seg1)-em(i).t(seg0))*86400;
    cseg_i = seglen > calseglen;
    cseg0 = seg0(cseg_i);
    cseg1 = seg1(cseg_i);
    Ncal = numel(cseg0);

    % Calculate mean measurements during calibration stop
    cal0.Q1 = zeros(Ncal,1);
    cal0.I1 = zeros(Ncal,1);
    cal0.C1 = zeros(Ncal,1);
    cal0.Q2 = zeros(Ncal,1);
    cal0.I2 = zeros(Ncal,1);
    cal0.C2 = zeros(Ncal,1);
    cal0.Q3 = zeros(Ncal,1);
    cal0.I3 = zeros(Ncal,1);
    cal0.C3 = zeros(Ncal,1);
    cal0.lat = zeros(Ncal,1);
    cal0.lon = zeros(Ncal,1);
    cal0.time = zeros(Ncal,1);
    cal0.spd = zeros(Ncal,1);
    for c=1:Ncal
        segi = cseg0(c):cseg1(c);
        cal0.I1(c) = mean(em(i).I1(segi), 'omitnan');
        cal0.Q1(c) = mean(em(i).Q1(segi), 'omitnan');
        cal0.C1(c) = mean(em(i).C1(segi), 'omitnan');
        cal0.I2(c) = mean(em(i).I2(segi), 'omitnan');
        cal0.Q2(c) = mean(em(i).Q2(segi), 'omitnan');
        cal0.C2(c) = mean(em(i).C2(segi), 'omitnan');
        cal0.I3(c) = mean(em(i).I3(segi), 'omitnan');
        cal0.Q3(c) = mean(em(i).Q3(segi), 'omitnan');
        cal0.C3(c) = mean(em(i).C3(segi), 'omitnan');
        cal0.lat(c) = mean(em(i).lat(segi), 'omitnan');
        cal0.lon(c) = mean(em(i).lon(segi), 'omitnan');
        cal0.time(c) = mean(em(i).t(segi), 'omitnan');
        cal0.spd(c) = mean(spd(segi), 'omitnan');
    end


    
    % Concatentate calibration data from each file into output data
    % structure
    calfields = fieldnames(cal0);
    if i==1
        for f = 1:numel(calfields)
            cal.(calfields{f}) = cal0.(calfields{f});
        end
    else
        for f = f:numel(calfields)
            cal.(calfields{f}) = [cal.(calfields{f}); ...
                                   cal0.(calfields{f})];
        end
    end    
    
    % % % % ** Diagnostics
    Nc = numel(cal0.lon);
    clab = num2str((1:Nc)','%02.0f');
    fig = figure;
    ax = wellspacedaxes(fig, 'Nrows', 1, 'Ncols', 2);
    line(em(i).lon, em(i).lat, 'parent', ax(1), 'color', 'blue');
    line(cal0.lon, cal0.lat, 'parent', ax(1), ...
                             'color', 'red', ...
                             'linestyle', 'none', ...
                             'marker', 'x');
    text(cal0.lon, cal0.lat, clab, 'parent', ax(1), ...
                                   'horizontalal', 'center', ...
                                   'verticalal', 'bottom');
                                        
    line(em(i).t(2:end), spd*100, 'parent', ax(2), ...
                                   'color', 'blue', ...
                                   'linestyle', '-', ...
                                   'marker', 'none');
    line(cal0.time, cal0.spd*100, 'parent', ax(2), ...
                                  'color', 'black', ...
                                  'linestyle', 'none', ...
                                  'marker', 'x');
    
    line(em(i).t, em(i).Q1, 'parent', ax(2), ...
                                           'color', 'red', ...
                                           'linestyle', '-', ...
                                           'marker', 'none');
    line(cal0.time, cal0.Q1, 'parent', ax(2), ...
                             'color', 'black', ...
                             'linestyle', 'none', ...
                             'marker', 'x');
    text(cal0.time, cal0.spd*100, clab, 'parent', ax(2), ...
                                        'horizontalal', 'center', ...
                                        'verticalal', 'bottom');
    text(cal0.time, cal0.Q1, clab, 'parent', ax(2), ...
                                   'horizontalal', 'center', ...
                                   'verticalal', 'bottom');
    


     notalot = 0;
                               
end


% Create ASCII file for calibration points corresponding to llQ file 
if ~exist('outfile', 'var')
    outfile = 'EMP400_calpoints.csv';
end
outfid = fopen(outfile,'wt');
fprintf(outfid, 'time,lat,lon,spd,Q1,I1,C1,Q2,I2,C2,Q3,I3,C3,\n');
outfstr = ['%s,%10.6f,%11.6f,%4.2f,%6.2f,%6.2f,%6.2f,' ...
                                  '%6.2f,%6.2f,%6.2f,' ...
                                  '%6.2f,%6.2f,%6.2f\n'];
for c=1:numel(cal.time)
    fprintf(outfid, outfstr, datestr(cal0.time(c), 'yyyy-mm-ddTHH:MM:SS'), ...
                             cal0.lat(c), ...
                             cal0.lon(c), ...
                             cal0.spd(c), ...
                             cal0.Q1(c), ...
                             cal0.I1(c), ...
                             cal0.C1(c), ...
                             cal0.Q2(c), ...
                             cal0.I2(c), ...
                             cal0.C2(c), ...
                             cal0.Q3(c), ...
                             cal0.I3(c), ...
                             cal0.C3(c));
end
fclose(outfid);
disp(['Written: ' outfile]);


end

