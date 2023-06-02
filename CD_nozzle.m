%% Constants

g = 1.4; % gamma constant for air

%% Given Values

pB = 14.7; %back pressure in psi
R1 = 0.49; %inlet radius in in
Rt = 0.2; %throat radius in in
R2 = 0.4; %exit radius in in
Lt = 1; %length from inlet to throat in in
L = 4; %length in in
pstar_p0 = 0.5283; %ratio of p* to p0
Tstar_T0 = 0.8333; %ratio of T* to T0

%% Calculation Area of the Nozzle

%create an equation to find area at any location
m1 = -.29/1; %slope of nozzle in converging
m2 = .2/3; %slope of nozzle in diverging
%A = pi*r^2

%% Choked Flow Conditions

R = 287; % ideal gas constant
P0 = 675686; %Po value in Pa
T0 = 300; %T0 value in K
Mt = 1; %Mach # at choked flow condition

%% Given Table 1 Values

points = [1 2 3 4 5 6 7 8 9 10 11];
location = [0 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];
A = [0 0 0 0 0 0 0 0 0 0 0];
M = [0 0 0 0 0 0 0 0 0 0 0];
r = [0 0 0 0 0 0 0 0 0 0 0];

%% Calculating Table 1 Values

%calculating area values for each point
r(1) = R1;
r(2) = R1 + (m1*.25);
r(3) = R1 + (m1*.5);
r(4) = R1 + (m1*.75);
r(5) = Rt;
r(6) = Rt+(m2*(location(6)-1));
r(7) = Rt+(m2*(location(7)-1));
r(8) = Rt+(m2*(location(8)-1));
r(9) = Rt+(m2*(location(9)-1));
r(10) = Rt+(m2*(location(10)-1));
r(11) = R2;

%% Conversions

%convert all inch values to mm
R1 = R1*25.4;
Rt = Rt*25.4;
R2 = R2*25.4;
Lt = Lt*25.4;
L = L*25.4;
r = r.*25.4;
location = location.*25.4;

%% Area

A0 = pi*R1^2; %area at inlet
At = pi*Rt^2; %area at throat
Ae = pi*R2^2; %area at exit, equal to A*

A = pi*(r.^2); %area at each location

%% Mach Number
%calculating area value at location 3 using choked flow conditions
A(5) = At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*Mt^2)^((g+1)/(2*(g-1))))/Mt);

%solving for Mach numbers

syms M
assume(M,'real')
M1 = double(vpasolve(A(1)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[0,1]));
M2 = double(vpasolve(A(2)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[0,1]));
M3 = double(vpasolve(A(3)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[0,1]));
M4 = double(vpasolve(A(4)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[0,1]));
M5 = 1;
M6 = double(vpasolve(A(6)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));
M7 = double(vpasolve(A(7)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));
M8 = double(vpasolve(A(8)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));
M9 = double(vpasolve(A(9)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));
M10 = double(vpasolve(A(10)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));
M11 = double(vpasolve(A(11)==At*(((g+1)/2)^(-(g+1)/(2*(g-1))))*(((1+((g-1)/2)*M^2)^((g+1)/(2*(g-1))))/M),M,[1,Inf]));

M = [M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11];

%% Other Values
%Calculating the pressure values

for i = 1:11
    p(i) = P0*(1+((g-1)/2)*M(i).^2)^(-g/(g-1));
end

%Calculating the temperature values

for i = 1:11
    T(i) = T0*(1+((g-1)/2)*M(i).^2)^(-1);
end

%Calculating the density (rho) values

%First we need to find rho0
% PV = mRT is ideal gas law
%P = rhoRT with density
%rho = P/RT
rho0 = P0/(R*T0);

%Now we can calculate the rho values
for i = 1:11
    rho(i) = rho0*(1+((g-1)/2)*M(i).^2)^(-1/(g-1));
end


%Calculating the velocity values

%first calculate the speed of sound
%then plug in and calculate the velocity
for i = 1:11
    a(i) = sqrt(g*R*T(i));
    V(i) = M(i)*a(i);
end


%Calculating mass flow rate

for i = 1:11
    mdot(i) = rho(i)*A(i)*V(i);
end

%% Organizing the Table Values

table = [points; location; A; p; T; rho; V; M; mdot]';


%% Quick Pressure Conversion

%converting Pa to psi
p = p.*0.000145038;

%printing the result
table2 = [location; p];

%% Experiment Data

%first we needed to average each column in each spreadsheet
%now we can import the data

%% Importing the 1_0 Data

Data_10 = readmatrix('1_0.xlsx');
T01_10 = Data_10(519,1); %Tank T01 (F)
T02_10 = Data_10(519,2); %Pitot T02 (F)
P01_10 = Data_10(519,3); %Tank P01 (Volts)
P02_10 = Data_10(519,4); %Pitot P02 (Volts)
Ps_10 = Data_10(519,5); %Pitot Ps (Volts)

%% Importing the 1_25 Data

Data_125 = readmatrix('1_25.xlsx');
T01_125 = Data_125(259,1); %Tank T01 (F)
T02_125 = Data_125(259,2); %Pitot T02 (F)
P01_125 = Data_125(259,3); %Tank P01 (Volts)
P02_125 = Data_125(259,4); %Pitot P02 (Volts)
Ps_125 = Data_125(259,5); %Pitot Ps (Volts)

%% Importing the 1_5 Data

Data_15 = readmatrix('1_5.xlsx');
T01_15 = Data_15(247,1); %Tank T01 (F)
T02_15 = Data_15(247,2); %Pitot T02 (F)
P01_15 = Data_15(247,3); %Tank P01 (Volts)
P02_15 = Data_15(247,4); %Pitot P02 (Volts)
Ps_15 = Data_15(247,5); %Pitot Ps (Volts)

%% Importing the 1_75 Data

Data_175 = readmatrix('1_75.xlsx');
T01_175 = Data_175(222,1); %Tank T01 (F)
T02_175 = Data_175(222,2); %Pitot T02 (F)
P01_175 = Data_175(222,3); %Tank P01 (Volts)
P02_175 = Data_175(222,4); %Pitot P02 (Volts)
Ps_175 = Data_175(222,5); %Pitot Ps (Volts)

%% Importing the 2_0 Data

Data_20 = readmatrix('2_0.xlsx');
T01_20 = Data_20(225,1); %Tank T01 (F)
T02_20 = Data_20(225,2); %Pitot T02 (F)
P01_20 = Data_20(225,3); %Tank P01 (Volts)
P02_20 = Data_20(225,4); %Pitot P02 (Volts)
Ps_20 = Data_20(225,5); %Pitot Ps (Volts)

%% Importing the 2_5 Data

Data_25 = readmatrix('2_5.xlsx');
T01_25 = Data_25(307,1); %Tank T01 (F)
T02_25 = Data_25(307,2); %Pitot T02 (F)
P01_25 = Data_25(307,3); %Tank P01 (Volts)
P02_25 = Data_25(307,4); %Pitot P02 (Volts)
Ps_25 = Data_25(307,5); %Pitot Ps (Volts)

%% Importing the 3_0 Data

Data_30 = readmatrix('3_0.xlsx');
T01_30 = Data_30(212,1); %Tank T01 (F)
T02_30 = Data_30(212,2); %Pitot T02 (F)
P01_30 = Data_30(212,3); %Tank P01 (Volts)
P02_30 = Data_30(212,4); %Pitot P02 (Volts)
Ps_30 = Data_30(212,5); %Pitot Ps (Volts)

%% Importing the 3_5 Data

Data_35 = readmatrix('3_5.xlsx');
T01_35 = Data_35(168,1); %Tank T01 (F)
T02_35 = Data_35(168,2); %Pitot T02 (F)
P01_35 = Data_35(168,3); %Tank P01 (Volts)
P02_35 = Data_35(168,4); %Pitot P02 (Volts)
Ps_35 = Data_35(168,5); %Pitot Ps (Volts)

%% Importing the 4_0 Data

Data_40 = readmatrix('4_0.xlsx');
T01_40 = Data_40(201,1); %Tank T01 (F)
T02_40 = Data_40(201,2); %Pitot T02 (F)
P01_40 = Data_40(201,3); %Tank P01 (Volts)
P02_40 = Data_40(201,4); %Pitot P02 (Volts)
Ps_40 = Data_40(201,5); %Pitot Ps (Volts)

%% Importing the 4_5 Data

Data_45 = readmatrix('4_5.xlsx');
T01_45 = Data_45(183,1); %Tank T01 (F)
T02_45 = Data_45(183,2); %Pitot T02 (F)
P01_45 = Data_45(183,3); %Tank P01 (Volts)
P02_45 = Data_45(183,4); %Pitot P02 (Volts)
Ps_45 = Data_45(183,5); %Pitot Ps (Volts)

%% Importing the 5_0 Data

Data_50 = readmatrix('5_0.xlsx');
T01_50 = Data_50(762,1); %Tank T01 (F)
T02_50 = Data_50(762,2); %Pitot T02 (F)
P01_50 = Data_50(762,3); %Tank P01 (Volts)
P02_50 = Data_50(762,4); %Pitot P02 (Volts)
Ps_50 = Data_50(762,5); %Pitot Ps (Volts)

%% Importing the CD Throat Data

Data_cdt = readmatrix('CD_throat.xlsx');
T01_cdt = Data_cdt(2:2397,1); %Tank T01 (F)
T02_cdt = Data_cdt(2:2397,2); %Pitot T02 (F)
P01_cdt = Data_cdt(2:2397,3); %Tank P01 (Volts)
P02_cdt = Data_cdt(2:2397,4); %Pitot P02 (Volts)
Ps_cdt = Data_cdt(2:2397,5); %Pitot Ps (Volts)

%calculate Mach number at throat
M_cdt = sqrt((((Ps_cdt./P02_cdt).^(-(g-1)/g))-1)*(2/(g-1)));

%% Collecting all Experiment Data

%concatenating all values of each variable into one vector
T01 = [T01_10 T01_125 T01_15 T01_175 T01_20 T01_25 T01_30 T01_35 T01_40 T01_45 T01_50];
T02 = [T02_10 T02_125 T02_15 T02_175 T02_20 T02_25 T02_30 T02_35 T02_40 T02_45 T02_50];
P01 = [P01_10 P01_125 P01_15 P01_175 P01_20 P01_25 P01_30 P01_35 P01_40 P01_45 P01_50];
P02 = [P02_10 P02_125 P02_15 P02_175 P02_20 P02_25 P02_30 P02_35 P02_40 P02_45 P02_50];
Ps = [Ps_10 Ps_125 Ps_15 Ps_175 Ps_20 Ps_25 Ps_30 Ps_35 Ps_40 Ps_45 Ps_50];

%% Converting Values to SI Units

%converting temperatures from Fahrenheit to Kelvin
T01 = ((T01-32).*(5/9))+273.15;
T02 = ((T02-32).*(5/9))+273.15;
T01_cdt = ((T01_cdt-32).*(5/9))+273.15;
T02_cdt = ((T02_cdt-32).*(5/9))+273.15;

%converting pressures from volts to psi to Pa
P01 = (P01*14.9912-0.0728).*6894.76;
P02 = (P02*14.9981-0.1234).*6894.76;
Ps = (Ps*14.99908+0.0279).*6894.76;
P01_cdt = (P01_cdt*14.9912-0.0728).*6894.76;
P02_cdt = (P02_cdt*14.9981-0.1234).*6894.76;
Ps_cdt = (Ps_cdt*14.99908+0.0279).*6894.76;

%% Calculating Mach number from pitot tube

%calculating the p/p02 values
ps_p02 = Ps./P02;

%solving for Mach number
%we will denote Mach number as N here
syms N
assume(N,'real')

N1 = double(vpasolve(ps_p02(1) == (1+((g-1)/2)*N^2)^(-g/(g-1)),N,[0,1]));
N2 = double(vpasolve(ps_p02(2) == (1+((g-1)/2)*N^2)^(-g/(g-1)),N,[0,1]));
N3 = double(vpasolve(ps_p02(3) == (1+((g-1)/2)*N^2)^(-g/(g-1)),N,[0,1]));
N4 = double(vpasolve(ps_p02(4) == (1+((g-1)/2)*N^2)^(-g/(g-1)),N,[0,1]));
%N5 = double(vpasolve(ps_p02(5) == (1+((g-1)/2)*N^2)^(-g/(g-1)),N,[1,Inf]));

N4 = double(vpasolve(ps_p02(4) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N5 = double(vpasolve(ps_p02(5) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N6 = double(vpasolve(ps_p02(6) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N7 = double(vpasolve(ps_p02(7) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N8 = double(vpasolve(ps_p02(8) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N9 = double(vpasolve(ps_p02(9) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N10 = double(vpasolve(ps_p02(10) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));
N11 = double(vpasolve(ps_p02(11) == (((g+1)/(1-g+(2*g*N^2)))*(((((g+1)^2)*N^2)/((4*g*N^2)-(2*(g-1))))^(-g/(g-1)))),N,[1,Inf]));

N1 = .1893;

M_exp = [N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11];

%% Normalizing Locations

location_norm = location/101.6; %normalized in mm

location2 = [0 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];
location2_norm = location2/4;

%% Plot #1

%creating new set of points with the values for the experiment values
location2 = [0 12.5 25 37.5 50 62.5 75 87.5 100 112.5 125];

figure(1);
plot(location_norm,M,'r-o','LineWidth',2,'MarkerSize', 5)
hold on;
plot(location2_norm,M_exp,'b-o','LineWidth',2,'MarkerSize', 5)
grid on;
title('M vs. Normalized Nozzle Distance');
xlabel('Normalized Nozzle Distance Along Centerline');
ylabel('Mach Number (M)');
legend('Analytical Values','Experimental Values');

%% Plot #2

figure(2);
plot(location2_norm,T01,'r','LineWidth',2,'MarkerSize', 15)
hold on;
plot(location2_norm,T02,'b','LineWidth',2,'MarkerSize', 15)
grid on;
title('T_{01} & T_{02} vs. Nozzle Distance');
xlabel('Normalized Nozzle Distance Along Centerline');
ylabel('Temperature (K)');
legend('T_{01}','T_{02}');

%% Plot #3

%Calculating pcdt/p01
p_ratio = Ps_cdt./P01_cdt;

figure(3);
plot(M_cdt,p_ratio,'r','LineWidth',2,'MarkerSize', 15)
hold on;
grid on;
title('p_{throat}/p_{01} vs. Mach Number at Nozzle Throat');
xlabel('Mach Number (M)');
ylabel('p_{throat}/p_{01}');

%Creating a line of best fit
coefficients = polyfit(M_cdt, p_ratio, 3);
xFit = linspace(min(M_cdt), max(M_cdt), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
equation = sprintf('y = %.2fx^3 + %.2fx^2 + %.2fx + %.2f', coefficients);
text(min(M_cdt)+0.32, max(p_ratio)-0.01, equation,'FontSize',10);

%Determining M value for pressure ratio
M_pratio = sqrt((((pstar_p0).^(-(g-1)/g))-1)*(2/(g-1)));

%Identifying critical pressure ratio on graph
yline(pstar_p0,'r--', 'LineWidth', 1.5);
plot(M_pratio,pstar_p0,'o-','MarkerSize', 10,'MarkerFaceColor','k','MarkerEdgeColor','k');
text(M_pratio+.01,pstar_p0+.02,'Critical P Ratio = 0.5283');

%% Plot #4

figure(4);
plot(M_cdt,T01_cdt,'r','LineWidth',2,'MarkerSize', 15)
hold on;
plot(M_cdt,T02_cdt,'b','LineWidth',2,'MarkerSize', 15)
grid on;
title('T_{01} & T_{02} vs. Mach Number at Nozzle Throat');
xlabel('Mach Number (M)');
ylabel('Temperature (K)');
legend('T_{01}','T_{02}');

%% Plot #5

figure(5);
plot(M_cdt,P01_cdt,'r','LineWidth',2,'MarkerSize', 15)
hold on;
plot(M_cdt,P02_cdt,'b','LineWidth',2,'MarkerSize', 15)
hold on;
plot(M_cdt,Ps_cdt,'g','LineWidth',2,'MarkerSize', 15)
grid on;
title('P_{01} P_{02} & P vs. Mach Number at Nozzle Throat');
xlabel('Mach Number (M)');
ylabel('Pressure (Pa)');
legend('P_{01}','P_{02}','P');