% AUTHOR: Bruce C. Bowler
function mbd2matlab(glider,mission,type)
global x_sci_salinity;
global x_sci_sigmat;
global eastbound;
global westbound;
global x_corrected_lat;
global x_corrected_lon;
global x_yo_number;
global x_solar_az
global x_solar_el
global x_corrected_chl

global eastYoRange;
global westYoRange;

if nargin == 2
    type = 'mbd';
end

switch type
    case 'sbd'
        flight = 's';
        science = 't';
    case 'mbd'
        flight = 'm';
        science = 'n';
    case 'dbd'
        flight = 'd';
        science = 'e';
end

if isunix
    lscomm = 'ls';
elseif ispc
    lscomm = 'dir/b';
end

tic;
cd ([vdrive() '/glider/' glider '/' mission '/' flight 'bd'])
% cd ([vdrive() '/glider/' glider '/' mission '/preflight']);% flight 'bd'])

missionNumber = sscanf(mission,'mission%f'); %converts string mission# to float

if missionNumber == 10 && strcmp(glider,'grampus')
    disp ('skipping merging, it was done manually due to issues with the mbdlist on the glider')
else
    if missionNumber < 8 && strcmp(glider,'henry')
        command = [lscomm ' *.' flight 'bd | dbd2asc -s | dba2_orig_matlab > ' mission '.mfile'];
        [~,~]=system(command);
        t = toc;
        fprintf(' (%1.2f seconds) \n%s',t,'data converted to .mat')
    else %sends unix or pc specific commands to command prompt, dbd2asc command
        command1 = [lscomm ' *.' flight 'bd | dbd2asc -s > ' mission '.' flight 'asc'];
        command2 = [lscomm ' *.' science 'bd | dbd2asc -s > ' mission '.' science 'asc'];
        system(command1);
        t = toc;
        fprintf(' (%1.2f seconds) \n%s',t,[flight 'bd processed'])
        system(command2);
        t = toc;
        fprintf(' (%1.2f seconds) \n%s',t,[science 'bd processed'])
        command3 = ['dba_merge ' mission '.' flight 'asc ' mission '.' science 'asc | dba2_orig_matlab > ' mission '.mfile'];
        system(command3);
        t = toc;
        fprintf(' (%1.2f seconds) \n%s',t,'data merged and converted to .mat')
    end
end

fid = fopen([mission '.mfile'],'r');
mfile = [pwd '/' fgetl(fid)];
fclose(fid);
run (mfile)
clear comm command stat mfile

lastElement = size(data,2);
x_sci_salinity  = lastElement + 1;
x_sci_sigmat    = lastElement + 2;
x_corrected_lat = lastElement + 3;
x_corrected_lon = lastElement + 4;
x_yo_number     = lastElement + 5;
x_solar_az      = lastElement + 6;
x_solar_el      = lastElement + 7;
x_corrected_chl = lastElement + 8;

%
% make sure data is sorted based on m_present_time
%
pt = data(:,m_present_time); %#ok<*NODEF>
[~,ix]=sort(pt);
data=data(ix,:);
%
% the bb2fl reports beta.  We really want bb.  This bit of code
% converts beta to bb and stores the result back into the data
% file
%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Converting beta to bb')
beta = data(:,sci_bbfl2s_bb_scaled);
if missionNumber == 12
    beta = beta * 10;   % there was an issue with the calibration data for bb during mission 12
end
bb = 2*pi*1.1*beta;
data(:,sci_bbfl2s_bb_scaled) = bb;
clear beta bb;
%
% 'corrrect' chl data based on bbfl2s chl vs discrete chls taken at either
% end of a mission.  (v:/gnats/glider_discrete.xlsx)
%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Correcting chl')
newChl = 1.8242*(data(:,sci_bbfl2s_chlor_scaled).^0.6146);
data(:,x_corrected_chl) = newChl;
clear newChl
%
% cvt mg nitrate to uM nitrate for grampus mission 5 only
%
if strcmp(glider,'grampus') && strcmp(mission,'mission5')
    t = toc;
    fprintf(' (%1,2f seconds) \n %s',t,'Converting mg nitrate to uM nitrrate')
    nitrate = data(:,sci_suna_nitrate_mg)/14.0067e-3;
    data(:,sci_suna_nitrate_mg) = nitrate;
    clear nitrate;
end
%
% calculate salinity and sigma T for inclusion in the data file
%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Calculating salinity and sigma t')
press = data(:,sci_water_pressure);
cond  = data(:,sci_water_cond);
temp  = data(:,sci_water_temp);
salt  = sw_salt(cond*10/sw_c3515,temp,press*10);
[~,sig] = swstate(salt,temp,press*10);
sig(sig<1) = nan;

data(:,x_sci_salinity) = salt;
data(:,x_sci_sigmat)   = sig;
clear press cond temp salt junk sig
%
% eliminate any bogus GPS fixes (lat or lon = 69696969
%
t=toc;
fprintf(' (%1.2f seconds) \n%s',t,'Eliminating bad fixes')
if exist('m_gps_lat','var')
    for i = [m_gps_lat m_gps_lon]
        ix = data(:,i)==69696969;
        data(ix,:)=[];
    end
end
%
% put m_lat, m_lon, m_gps_lat and m_gps_lon into a "computer friendly"
% format (ie degree and decimal degree)
%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Converting latitude and longitude to decimal')
for i=[m_lat m_lon]
    x = data(:,i);
    deg = fix(x/100);
    minutes = x-deg*100;
    x = deg+minutes/60;
    data(:,i) = x;
end
if exist('m_gps_lat','var')
    for i=[m_gps_lat m_gps_lon]
        x = data(:,i);
        deg = fix(x/100);
        minutes = x-deg*100;
        x = deg+minutes/60;
        data(:,i) = x;
    end
end
clear x deg minutes i

%
t = toc;
if exist('m_gps_lat','var')
    fprintf(' (%1.2f seconds) \n%s',t,'Correcting latitude and longitude for set/drift')
    
    lats = data(:,m_lat); lons = data(:,m_lon); mtime = data(:,m_present_time);
    glats = data(:,m_gps_lat); % glons = data(:,m_gps_lon);
    
    newlats = zeros(numel(lats),1);
    newlons = newlats;
    
    ix = isfinite(lats) & isfinite(glats);
    iy = find(ix == 1);
    bb = [1; diff(iy)];
    iz = find(bb>10);
    
    slope = zeros(numel(iz),2);
    incpt = zeros(numel(iz),2);
    range = zeros(numel(iz),3);
    
    initialGpsIndex = [1; iy(iz)];
    
    for i = 2:numel(initialGpsIndex)
        ix1  = initialGpsIndex(i-1);
        ix3  = initialGpsIndex(i);
        cidx = find(isfinite(lats));
        ivec = setdiff(cidx,ix3);
        ivec = ivec(ivec<ix3);
        ix2  = ivec(end);
        
        lat1 = lats(ix2);
        lat2 = lats(ix3);
        lon1 = lons(ix2);
        lon2 = lons(ix3);
        tm1  = mtime(ix1);
        tm2  = mtime(ix2);
        deltalat = lat2-lat1;
        deltalon = lon2-lon1;
        deltatm  = tm2 - tm1;
        %     if deltatm == 0; disp([i ix1 ix2 ix3]); end
        range(i-1,1) = tm1;
        range(i-1,2) = tm2;
        range(i-1,3) = tm2-tm1;
        slope(i-1,1) = deltalat/deltatm;
        slope(i-1,2) = deltalon/deltatm;
        incpt(i-1,1) = 0;
        incpt(i-1,2) = 0;
    end
    clear bb cidx deltalat deltalon lat1 lat2 lon1 lon2 tm1 tm2 deltatm
    clear initialGpsIndex ix ix1 ix2 ix3 iy iz ivec glats
    
    
    for i = 1:numel(lats)
        thisTime = mtime(i);
        for j = 1:numel(range)/3
            if range(j,1) <= thisTime && thisTime < range(j,2)
                k = j;
                ldt  = thisTime - range(k,1);
                dLat = slope(k,1)*ldt + incpt(k,1);
                dLon = slope(k,2)*ldt + incpt(k,2);
                newlats(i) = lats(i) + dLat;
                newlons(i) = lons(i) + dLon;
            end
        end
    end
    clear i thisTime j mtime range k ldt dLat dLon lats lons slope incpt
    newlats(newlats==0) = nan;
    newlons(newlons==0) = nan;
else
    newlats = data(:,m_lat);
    newlons = data(:,m_lon);
end

newlats = naninterp(newlats);
newlons = naninterp(newlons);

if exist('x_measured_depth','var')
    depth = data(:,x_measured_depth);
    depth = naninterp(depth);
    data(:,x_measured_depth) = depth;
else
    depth = data(:,m_depth);
    depth = naninterp(depth);
    data(:,m_depth) = depth;
end

data(:,x_corrected_lat) = newlats;
data(:,x_corrected_lon) = newlons;
clear newlats newlons depth

%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Tag each upcast with number')

endPoints = findRuns(data(:,sci_bbfl2s_bb_scaled),10);
endPoints(:,1) = endPoints(:,1)-1;
endPoints(:,2) = endPoints(:,2)+1;
numPoints = numel(endPoints);
ep1 = reshape(endPoints',numPoints,1);
ep1 = ep1(2:end-1);
numPoints = numPoints - 2;
endPoints = reshape(ep1,2,numPoints/2)';
numPoints = size(endPoints,1);

for i = 1:numPoints
    data((endPoints(i,1):endPoints(i,2)),x_yo_number) = i;
end
clear numPoints endPoints i ep1

t=toc;
fprintf(' (%1.2f seconds) \n%s',t,'Add solar azimuth and elevation data')

matlabDate = datetime(data(:,m_present_time)/86400+datenum(1970,1,1),'convertfrom','datenum');
% for i = 1:size(data,1)
%     if fix(i/500)*500 == i
%         disp(num2str(i))
%     end
%     [az,el]=calcSolarData(matlabDate(i),data(i,x_corrected_lat),data(i,x_corrected_lon));
%     data(i,x_solar_az)  = az;
%     data(i,x_solar_el) = el;
% end

[az, el]= SolarAzEl(matlabDate,data(:,x_corrected_lat),data(:,x_corrected_lon),0);
az(az>180) = az(az>180) - 360;
data(:,x_solar_az) = az;
data(:,x_solar_el) = el;
clear matlabDate az el %i

%
% split data into eastbound and westbound segments, based on
% farthest east gps fix
%
t = toc;
fprintf(' (%1.2f seconds) \n%s',t,'Splitting data into eastbound and westbound legs')
if exist('m_gps_lon','var')
    [~,id]=max(data(:,m_gps_lon));
else
    [~,id]=max(data(:,m_lon));
end
eastbound = data(1:id,:);
westbound = data(id:end,:);
clear dists md id junk

eastYoRange = unique(eastbound(:,x_yo_number));
westYoRange = unique(westbound(:,x_yo_number));
eastYoRange = eastYoRange(2:end);
westYoRange = westYoRange(2:end);

t = toc;
fprintf(' (%1.2f seconds)\n',t)
clear t pt command1 command2 command3 fid
save(['../matlabdata/' mission '.mat'])

function endPoints = findRuns(array,runLen)
%
% a is a vector of data with nans in it.  The goal is to identify long
% stretches of nans as these will separate one data segment from the next
% "long" is relative and depends on circumstances.
%
badvalue = -9999; % known bad value
array(~isfinite(array)) = badvalue; % change it to a known "impossible" value

p = find([true;diff(array)~=0;true]);
q = find(diff(p)>=runLen);
%
%             +-------start of run of "bad values"
%             |      +end of run of "bad values"
%             v      v
endPoints = [p(q) p(q+1)-1];