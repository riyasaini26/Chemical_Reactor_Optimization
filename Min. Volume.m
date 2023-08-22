% Prompt user for input
%num_elements_CAo = input('Enter the number of elements in the array(CAo): ');

% Preallocate array
%my_array_CAo = zeros(1, num_elements_CAo);
my_array_CAo =[2,5,6,6,11,14,16,24];
% Prompt user for values and store in array
% for i = 1:num_elements_CAo
%     my_array_CAo(i) = input(['Enter value for element(CAo) ' num2str(i) ': ']);
% end
% Preallocate array
%my_array_Ca = zeros(1, num_elements_CAo);
my_array_Ca =[0.5,3,1,2,6,10,8,4];
% Prompt user for values and store in array
% for i = 1:num_elements_CAo
%     my_array_Ca(i) = input(['Enter value for element(Ca) ' num2str(i) ': ']);
% end

% Preallocate array
%my_array_tow = zeros(1, num_elements_CAo);
my_array_tow =[30,1,50,8,4,20,20,4];
% Prompt user for values and store in array
% for i = 1:num_elements_CAo
%     my_array_tow(i) = input(['Enter value for element(tow) ' num2str(i) ': ']);
% end
my_array_rate = zeros(1, 8);
for i = 1:8
    my_array_rate(i) = my_array_tow(i)/(my_array_CAo(i)-my_array_Ca(i));
end
%n = 3; % Degree of the polynomial
%coefficients = polyfit(my_array_Ca, my_array_rate, n); % Compute the polynomial coefficients
x_interp = linspace(min(my_array_Ca), max(my_array_Ca), 10000); % Generate 1000 evenly spaced x-values for interpolation
y_interp=spline(my_array_Ca, my_array_rate,x_interp);
%y_interp = polyval(coefficients, x_interp)+1.8; % Evaluate the polynomial at the x-values
plot(x_interp, y_interp);

% Input parameters
flowrate = input('Enter flow rate (m^3/h): ');
C0 = input('Enter initial concentration of organics (mg/L): ');
Ct = input('Enter required concentration in treated stream (mg/L): ');
reactor_system_type = input('Enter reactor system type (1-5): ');

% % Constants
% Vmin = linspace(1,Vmax,100); % Minimum volume
% Vmax = 1000; % Maximum volume
tolerance = 1e-6; % Tolerance for convergence

% Function to calculate concentration in treated stream for different reactor systems
f = @(V) reactor_system(flowrate, C0, Ct, V, reactor_system_type);
area = zeros(1, 1000);
% Calculate minimum volume for different reactor systems
 min_y = min(y_interp); % find minimum value of y
        min_y_indices = find(y_interp == min_y); % find indices of elements in y that are equal to min_y
        min_x = x_interp(min_y_indices(1)); % find the corresponding value of x using the first index in min_y_indices
switch reactor_system_type
    case 1 % Single PFR
        area_under_plot = trapz(x_interp, y_interp);
        Vmin=area_under_plot*flowrate;
%         Vmin = fminbnd(f, Vmin(1), Vmax);
    case 2 % Single CSTR
        Vmin = (C0-Ct)*flowrate*(polyval(coefficients, Ct));
    case 3 % Two stirred tanks of any size
        ci=1;
        for i=1:1000
            area1=(ci-Ct)*(polyval(coefficients, Ct));
            area2=(C0-ci)*(polyval(coefficients, ci));
            area(i)=area1+area2;
            ci=ci+0.1;
        end
        min_area=min(area);
        Vmin=min_area*flowrate;
        %Vmin = fsolve(f, [Vmin(1), Vmin(1)]);
    case 4 % Combination of a PFR and a MFR
        min_y = min(y_interp); % find minimum value of y
        min_y_indices = find(y_interp == min_y); % find indices of elements in y that are equal to min_y
        min_x = x_interp(min_y_indices(1)); % find the corresponding value of x using the first index in min_y_indices
        area1=(C0-min_x)*y_interp(min_y_indices);
        x_start=Ct;
        x_end=min_x;
        [temp,idx_start]=min(abs(x_interp-x_start));
        [temp,idx_end]=min(abs(x_interp-x_end));
%         idx_start = find(x_interp == x_start);
%         idx_end = find(x_interp == x_end);

        % Extract the corresponding y values for the specified range
        y_range = y_interp(idx_start:idx_end);

        % Calculate the area under the curve using the trapz function
        area2 = trapz(x_interp(idx_start:idx_end), y_range);
        Vmin=(area1+area2)*flowrate;
        %Vmin = fsolve(f, [Vmin(1), Vmin(1)]);
    case 5 % PFR with recycle
         
         [temp,C0_indices]=min(abs(x_interp-C0));
         [temp,Ct_indices]=min(abs(x_interp-Ct));
         %indx=find(yy==min(yy));
%         i=find(xx==Cai);
%             [temp1, i] = min(abs(xx - Cai));
%             optimum = 10;
            for k=min_y_indices+3:C0_indices
%                 k = 500;
                h=y_interp(k);
                
                [temp,m]=min(abs(y_interp(Ct_indices:min_y_indices)-h));
                a1=trapz(x_interp(Ct_indices:m),y_interp(Ct_indices:m)-h);
                a2=h*(x_interp(k)-x_interp(m))-trapz(x_interp(m:k),y_interp(m:k));
                if abs(a1-a2)<0.1
                    %optimum = k;
                    break
                end
            end
            Vmin=x_interp(k);
%         % for i=min_y_indices+3:C0_indices
%         i=550;
%              h=y_interp(i);
%             % [temp,m]=min(abs(y_interp(Ct_indices:i-3)-h));
%              m=find(y_interp(Ct_indices:i-3)==h);
%              % Extract the corresponding y values for the specified range
%             y_range = y_interp(m:i);
% 
%             % Calculate the area under the curve using the trapz function
%             area1 = trapz(x_interp(m:i), y_range);
%             area1=h*(i-m)-area1;
%             y_range = y_interp(Ct_indices:m);
% 
%             % Calculate the area under the curve using the trapz function
%             area2 = trapz(x_interp(Ct_indices:m), y_range);
%             area2=area2-h*(m-Ct_indices);
%             %fprintf(abs(area1-area2));
%             
%             if((abs(area1-area2))<1)
%                % break;
%             end
%         % end
%          
        % Define the y value for which you want to find the corresponding x values
%         y_value = polyval(coefficients,2); % specify the y value
% 
%         % Find the indices of x values that correspond to the specified y value
%         idx = find(y_interp == y_value);
% 
%         % Extract the corresponding x values using logical indexing
%         x_values = x_interp(idx);
        %Vmin = fsolve(f, [Vmin(1), Vmin(1)]);
    otherwise
        error('Invalid reactor system type.');
end
% Display results
fprintf('Minimum volume required for reactor system %d: %.2f m^3\n', reactor_system_type, Vmin);