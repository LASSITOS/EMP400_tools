function plot_LASSITOSemp400profiles

empath = ['/Users/amahoney/Google Drive/My Drive/SyncFolders/LASSITOS' ...
          filesep 'FlightTestUTQ2023/EMP400data/PloverPtTransects'];

emilist = filedigger(empath, '*.EMI');
Nf = numel(emilist);


freq2plot = (1:16)*1000;
Nfreq = numel(freq2plot);

lat0 = 71.36366802;
lon0 = -156.377742;
lat1 = 71.36831101;
lon1 = -156.372794;
mstruct = createmstruct_AKAlbNAD83();
[x0,y0] = projfwd(mstruct, lat0,lon0);
[x1,y1] = projfwd(mstruct, lat1,lon1);
laz = atan2(y1-y0, x1-x0);



Q = cell(Nfreq,1);
C = cell(Nfreq,1);
d = cell(Nfreq,1);
 
for f = 1:Nf
    em = emp400_reademifile([empath filesep emilist{f}]);

    fprintf([emilist{f} ': %dkHz, %dkHz, %dkHz\n'],  em.info.frequencies/1000);
    
    for i=1:3
        freqi = em.info.frequencies(i);
        ff = find(freq2plot == freqi);
        Qfield = ['Q' num2str(i, '%d')];
        Cfield = ['C' num2str(i, '%d')];
        fprintf('%dkHz: ', freqi);
        if isempty(Q{ff})
            Q{ff} = em.(Qfield);
            C{ff} = em.(Cfield);
            [x, y] = projfwd(mstruct, em.lat, em.lon);
            d{ff} = (x-x0)*cos(laz) + (y-y0)*sin(laz);
            fprintf('new\n');
        else
            [x, y] = projfwd(mstruct, em.lat, em.lon);
            dnew = (x-x0)*cos(laz) + (y-y0)*sin(laz);
            dgap0 = max(diff(d{ff}));
            dgap1 = max(diff(dnew));
            fprintf('old - ');
            if dgap0 > dgap1
                fprintf('replaced with new\n');
                Q{ff} = em.(Qfield);
                C{ff} = em.(Cfield);
                [x, y] = projfwd(mstruct, em.lat, em.lon);
                d{ff} = (x-x0)*cos(laz) + (y-y0)*sin(laz);
            else
                fprintf('kept old\n');
            end
        end
    end
end

fig = figure;
ax = axes('parent', fig, 'pos', [0.1, 0.1, 0.85, 0.85]);

lcols = generatecolorramp({'red','cya','blu'}, Nfreq);
for ff = 1:Nfreq
    freqname = sprintf('%0.0f kHz', freq2plot(ff)/1000);
    line(d{ff}, C{ff}, 'parent', ax, 'displayname', freqname, ...
                                     'color', lcols(ff,:), ...
                                     'linewidth', 2);
end

set(ax, 'xlim', [-50, 550], 'xtick', -50:50:550, ...
        'ylim', [-10,900], 'ytick', 0:100:900, ...
        'xgrid', 'on', 'ygrid', 'on', 'fontsize', 20);
xlabel(ax, 'Distance from tide crack (m)', 'fontsize', 24);
ylabel(ax, 'Quadrature', 'fontsize', 24);

legend(ax, 'location', 'EastOutside');

end

