clear all
% clf

% fpath = '/Users/lnlammers/Desktop/Active Projects/CEC project/Calculations/Trace_gas_flux';
outfile = '/home/greg/data/rawdata/MojaveCarbon/FLV_N2Oflux_inverse_commonatmboundary.csv'

fid = fopen('/home/greg/data/gdrive_berkeley/MojaveCarbon/Data/TG_analysis_Laura/trace_gas_N2O_2.txt');

% Catm is the atmospheric boundary condition - some fluxes are VERY
% sensitive to this
%Catm = 1.834;%2.163; %CH4
%molecular_mass = 12.01 + 4*1.008; % CH4

%Catm = 398; %CO2 382
%molecular_mass = 12.01 + 2*15.998; % CO2

Catm = .328;%0.396; %N2O
molecular_mass = 44.013; % N2O

data = textscan(fid,'%s %s %f %f %f %f %f %f %f');
soil = data{:,1}; % soil ID
date = data{:,2}; % collection date
depth_all = data{:,3};% depth (cm)
c_all = data{:,4}; %concentration (ppm)
D_all = data{:,7}; %diffusivity Moldrop 1999 eq8 cm^2/s
% D_all = data{:,6}; %diffusivity Moldrop 1999 eq7 cm^2/s
% D_all = data{:,8}; %diffusivity Penman 1940 cm^2/s
% D_all = data{:,9}; %diffusivity Millington 1959 cm^2/s

soils = {'A1', 'A2','A3','B','C','D'};
dates = {'0513','0813','1113','0414','0514'};
net_consumption_kgHayr = zeros(length(soils),length(dates));
d1_consumption = zeros(length(soils),length(dates));
d2_consumption = zeros(length(soils),length(dates));
d3_consumption = zeros(length(soils),length(dates));
d4_consumption = zeros(length(soils),length(dates));
surface_flux = zeros(length(soils),length(dates));
net_consumption_ron = zeros(length(soils),length(dates));
surface_flux_ron = zeros(length(soils),length(dates));
j = 1;
while j <= length(soil)
    for k = 1:5
        if (j+k)>length(soil)
            increment = k-1;
        elseif (depth_all(j+k) < depth_all(j+(k-1)))
            increment = k-1;         
        end
    end
    if soil{j} == soils{1}
        m = 1;
    elseif soil{j} == soils{2}
        m = 2;
    elseif soil{j} == soils{3}
        m = 3;
    elseif soil{j} == soils{4}
        m = 4;
    elseif soil{j} == soils{5}
        m = 5;
    elseif soil{j} == soils{6}
        m = 6;
    end
    if date{j} == dates{1}
        n = 1;
    elseif date{j} == dates{2}
        n = 2;
    elseif date{j} == dates{3}
        n = 3;
    elseif date{j} == dates{4}
        n = 4;
    elseif date{j} == dates{5}
        n = 5;
    end
    
    c = c_all(j:j+increment)*(molecular_mass/24.45)*1000/(100^3)/(10^6); % g/cm^3
    
    % To fix the surface boundary condition to the atmospheric value
    c(1) = Catm*(molecular_mass/24.45)*1000/(100^3)/(10^6); % g/cm^3 

    D = D_all(j:j+increment);
    depth = depth_all(j:j+increment);
      % discretize grid
    dz = 0.1; % cm
    z = min(depth):dz:(max(depth)+dz);
    c_CH4 = zeros(1,length(z));
    c_function = fit(depth,c,'smoothingspline','SmoothingParam',.1);
    D_function = fit(depth,D,'smoothingspline','SmoothingParam',.1);
    c_CH4 = c_function(z); % Now gas concentration is distributed along
                           % the new depth grid
                           % Diffusivity could be too, but a step function
                           % is used instead (below)

% %     % linear segments
%     for i = 1:length(depth)-1
%         ytemp = polyfit([depth(i),depth(i+1)],[c(i),c(i+1)],1);
%         for k = 1:length(z)
%             if z(k)<= depth(i+1) && z(k)>=depth(i)
%                 c_CH4(k) = polyval(ytemp,z(k));
%             end
% %             if(z(k)==depth(i))
% %                 c_CH4(k) = c(i);
% %             end
%             
%         end
%     end
    
    % Ron's calculation
    clear fluxron consumptionron net_consumptionron
    fluxron = zeros(1,length(depth));
    for i = 2:length(depth)
        fluxron(i-1) = -D(i)*(c(i)-c(i-1))/(depth(i)-depth(i-1));
    end
    for i = 2:length(fluxron)
        consumptionron(i-1) = (fluxron(i)-fluxron(i-1))/(depth(i)-depth(i-1));
        net_consumptionron(i-1) = -fluxron(i)+fluxron(i-1);
    end
        net_consumption_ron(m,n) = sum(net_consumptionron)/1000*100^2*10000*60*60*24*365; %kg/Ha/yr
        surface_flux_ron(m,n) = fluxron(1)/1000*100^2*10000*60*60*24*365; % %kg/Ha/yr

      
    figure(1)
%     subplot(2,1,1)
    hold on
%     plot(c_function,depth,c)
    plot(c_CH4(1:end-1),z(1:end-1))
    plot(c,depth,'k*')
    xlabel('Concentration, g/cm^3')
    ylabel('depth, cm')
    set(gca,'YDir','reverse')
  
%     D_CH4 = D_function(z);
    consumption = zeros(1,length(z)-2);
    flux = zeros(1,length(z)-2);
    D_CH4 = zeros(1,length(z)-2);
    for i = 1:length(z)-2
            if z(i) < depth(2) %Diffusivity step-function
                D_CH4(i) = D(2);
            elseif z(i)>=depth(2) && z(i)<depth(3)
                 D_CH4(i) = D(3);
            elseif z(i)>=depth(3) && z(i)<depth(4)
                D_CH4(i) = D(4);
            elseif length(depth)>4 && z(i)>=depth(4) && z(i)<depth(5)
                D_CH4(i) = D(5);
            end
    end
    
% Plot Diffusivity vs. Depth
%     figure(1)
%     subplot(2,1,2)
%     hold on
%     plot(D_CH4,z(1:end-2))
%     xlabel('Diffusivity, cm^2/s')
%     ylabel('depth, cm')
    % For each depth in the grid calculated the gradient
    for i = 2:length(z)-2
        % Gradient for each interval in profile (First deriv)
        dc_dz = (c_CH4(i+1)-c_CH4(i))/dz;
        % Second derivative limit calculation
        d2c_dz2 = (c_CH4(i+1) - 2*c_CH4(i) + c_CH4(i-1))/dz^2;
        % Calculate flux and consumption (flux * diffusivity)
        flux(i) = -D_CH4(i)*dc_dz; %flux into the soil
        consumption(i) = -D_CH4(i)*d2c_dz2; % g/cm^3/s
    end
    
    % Surface flux is the flux across the uppermost layer 
    surface_flux(m,n) = flux(2)/1000*100^2*10000*60*60*24*365; % %kg/Ha/yr
    % Net consumption is the sum of consumption along the entire profile
    net_consumption = -sum(consumption*dz);% g/cm^2/s
    d1_consumption(m,n) = -(consumption(depth(2)/dz)*dz)/1000*100^2*10000*60*60*24*365;
    d2_consumption(m,n) = -(consumption(depth(3)/dz)*dz)/1000*100^2*10000*60*60*24*365;
    d3_consumption(m,n) = -(consumption(depth(4)/dz)*dz)/1000*100^2*10000*60*60*24*365;
    if length(depth) > 4
        d4_consumption(m,n) = -(consumption(depth(5)/dz)*dz)/1000*100^2*10000*60*60*24*365;
    else
        d4_consumption(m,n) = NaN;
    end
    net_consumption_kgHayr(m,n) = net_consumption/1000*100^2*10000*60*60*24*365; %kg/Ha/yr    % 10000 m^2 = 1 Ha

    figure(2)
    hold on
    plot(depth(2:end),net_consumptionron)
    
    j = j+increment+1;
end
figure(3)
hold on
plot(surface_flux,net_consumption_kgHayr,'r*')
plot(surface_flux_ron,net_consumption_ron,'b*')
xlabel('Surface flux, kg/Ha/yr')
ylabel('Flux by profile integration, kg/Ha/yr')
fprintf('%f\t%f\t%f\t%f\t%f\n',surface_flux')
fprintf('\n \n')
fprintf('%f\t%f\t%f\t%f\t%f\n',net_consumption_kgHayr')


sf_table = array2table(surface_flux', 'VariableNames', soils)
sf_table.date = datestr(datenum(dates', 'mmyy'))
sf_table.flux = {'surface', 'surface', 'surface', 'surface', 'surface'}'
cs_table = array2table(net_consumption_kgHayr', 'VariableNames', soils)
cs_table.date = datestr(datenum(dates', 'mmyy'))
cs_table.flux = {'netconsump', 'netconsump', 'netconsump', 'netconsump', 'netconsump'}'
d1_table = array2table(d1_consumption', 'VariableNames', soils)
d1_table.date = datestr(datenum(dates', 'mmyy'))
d1_table.flux = {'d1consump', 'd1consump', 'd1consump', 'd1consump', 'd1consump'}'
d2_table = array2table(d2_consumption', 'VariableNames', soils)
d2_table.date = datestr(datenum(dates', 'mmyy'))
d2_table.flux = {'d2consump', 'd2consump', 'd2consump', 'd2consump', 'd2consump'}'
d3_table = array2table(d3_consumption', 'VariableNames', soils)
d3_table.date = datestr(datenum(dates', 'mmyy'))
d3_table.flux = {'d3consump', 'd3consump', 'd3consump', 'd3consump', 'd3consump'}'
d4_table = array2table(d4_consumption', 'VariableNames', soils)
d4_table.date = datestr(datenum(dates', 'mmyy'))
d4_table.flux = {'d4consump', 'd4consump', 'd4consump', 'd4consump', 'd4consump'}'

out = [sf_table; cs_table; d1_table; d2_table; d3_table; d4_table]

writetable(out, outfile ,'Delimiter', ',') 
