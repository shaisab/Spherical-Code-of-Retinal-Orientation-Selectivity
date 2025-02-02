% Question 1

%1.a

clear all; 

%1.b

doc which

%1.c
path_to_doc = cd ('C:\Program Files\MATLAB\R2019b\toolbox\matlab\helptools\doc.m');

%1.d
current_path = pwd;

%1.e
whos

% Question 2

%2.a

%ans =

%   1.0000e-08

%2.b

%5e-5 = 0.00005
%3e11 = 300000000000

% Question 3

Initial_investment = 100000; % Initial amount of money Dorin invested

time5 = 5; % time period elapsed(in years)
time40 = 40;

growth = 1.05; % the growth of the market each year

% 3.a
Expected_return_5th_Year = Initial_investment .* growth .^ time5;

% 3.b
Expected_return_40th_Year = Initial_investment .* growth .^ time40;

% 3.c
Management_Fee = 0.005;

Adjusted_growth = growth - Management_Fee;

Adjusted_Expected_return_5th_Year = Initial_investment .* Adjusted_growth .^ time5;
Adjusted_Expected_return_40th_Year = Initial_investment .* Adjusted_growth .^ time40;

% 3.d
New_Management_Fee = 0.004;
New_Adjusted_growth = growth - New_Management_Fee; 

New_Expected_return_40th_Year = Initial_investment .* New_Adjusted_growth .^ time40;
Gain = New_Expected_return_40th_Year - Adjusted_Expected_return_40th_Year;

% Question 4
a = 3;
b = 5;
c = 7;

ans4a = a + b < c;
ans4b = ((a > b) & (a <= c));
ans4c = a | b;
ans4d = b ~= c;

clear all; 