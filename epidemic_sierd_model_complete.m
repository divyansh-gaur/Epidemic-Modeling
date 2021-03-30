% This script simulates the spread of an epidemic using the SEIRD or Susceptible-Exposed-Infected-Recovered-Dead model used in epidemiology

%% Parameters

incubation = 6;  % Time taken for exposed individual to show symptoms
rec_time = 14; % Time taken for infected individual to recover

r0 = 2.5; % The reproduction rate of the pathogen. One infected person infects r0 susceptible people.
pop = 5e6; % Total population of the region

sigma = r0/(pop*rec_time);

period = 365; % The total number of days for the simulation
int_day = 120; % The day on which preventive measures are put into place
vacc = 0; % Number of people receiving vaccines everyday

%% Modelling through the Differential equations

t_span = 0:1:period; % Duration of the model and time units

y_init = [pop-165, 3, 150, 11, 1]; % Initial Conditions [Susceptible, Exposed, Infected, Recovered, Dead]

[t,y]=ode45(@(t,y) diff_eq_solver(t,y,sigma, vacc), t_span, y_init); % Solves the system of ordinary differential eqations

%% Plot
figure(1)
title('Spread of the epidemic without intervention policies')    
xlabel('Days since initial exposure')
ylabel('Population')
hold on

susceptible = animatedline('Color','b', 'LineWidth', 1.6, 'MarkerSize', 12);
exposed = animatedline('Color','y', 'LineWidth', 1.6, 'MarkerSize', 12);
infectious = animatedline('Color','r', 'LineWidth', 1.6, 'MarkerSize', 12);
recovered = animatedline('Color','g', 'LineWidth', 1.6, 'MarkerSize', 12);
death = animatedline('Color','k', 'LineWidth', 1.6, 'MarkerSize', 12);

legend('Susceptible','Exposed','Infectious','Recovered','Death', 'Location', 'East')
grid minor;
set(gca, 'FontSize', 12)

for j = 1:period
    set(gca, 'XLim', [0 length(y)],'YLim', [0 pop])
    addpoints(susceptible,j,y(j,1));
    addpoints(exposed,j,y(j,2));
    addpoints(infectious,j,y(j,3));
    addpoints(recovered,j,y(j,4));
    addpoints(death,j,y(j,5));
    hold on
    drawnow limitrate
    pause(0.0001)
end

 hold off
 pause


%% Prevention Modelling

% This model is adjusted to show the effects of preventive measures

t_span1 = 0:1:int_day; % Time before the introduction of preventive measures
y_init = [pop-165, 3, 150, 11, 1]; % Initial Conditions [Susceptible, Exposed, Infected, Recovered, Dead]

[t,y]=ode45(@(t,y) diff_eq_solver(t,y,sigma, vacc), t_span1, y_init); % Solves the system of ordinary differential equations

y_2nd = [y(end,1), y(end,2), y(end,3), y(end,4), y(end,5)];

vacc = 100; % Vaccines are introduced and 'vacc' number of people are receiving them per day thereby being moved from susceptible to recovered

t_span2 = int_day+1:1:period; % Time after the introduction of preventive measures
[t,y2]=ode45(@(t,y) diff_eq_solver(t, y, sigma/2, vacc), t_span2, y_2nd); % Solves the system of ordinary differential eqations with halved r0 as a result of preventive measures

y_c = [y;y2];

%% Plot
figure(2)
title('Spread of the epidemic with intervention policies')    
xlabel('Days since initial exposure')
ylabel('Population')
hold on

susceptible_int = animatedline('Color','b', 'LineWidth', 1.6, 'MarkerSize', 12);
exposed_int = animatedline('Color','y', 'LineWidth', 1.6, 'MarkerSize', 12);
infectious_int = animatedline('Color','r', 'LineWidth', 1.6, 'MarkerSize', 12);
recovered_int = animatedline('Color','g', 'LineWidth', 1.6, 'MarkerSize', 12);
death_int = animatedline('Color','k', 'LineWidth', 1.6, 'MarkerSize', 12);

legend('Susceptible','Exposed','Infectious','Recovered','Death', 'Location', 'West')
grid minor;
set(gca, 'FontSize', 12)

for i = 1:period
    set(gca, 'XLim', [0 length(y_c)],'YLim', [0 pop])
    addpoints(susceptible_int,i,y_c(i,1));
    addpoints(exposed_int,i,y_c(i,2));
    addpoints(infectious_int,i,y_c(i,3));
    addpoints(recovered_int,i,y_c(i,4));
    addpoints(death_int,i,y_c(i,5));
    hold on
    drawnow limitrate
    pause(0.0001)
end

hold off

 
%% Function for solving the differential equations
 
 function dydt = diff_eq_solver(t, y,sigma, vacc)

% Parameters
drate = 0.034;  % Case fatality rate
incubation = 6;  % Time taken for exposed individual to show symptoms
rec_time = 14; % Time taken for infected individual to recover

% Initial Conditions
S = y(1);
E= y(2);
I = y(3);

% System of Differential Equations
dS = -sigma*I.*S-vacc;
dE = sigma*I.*S -  (1/incubation).*E;
dI = (1/incubation)*E - (1/rec_time)*I;
dR = (1/rec_time)*(1-drate)*I+vacc;
dD = (drate)*(1/rec_time)*I;

dydt = [dS; dE; dI; dR; dD];

end
