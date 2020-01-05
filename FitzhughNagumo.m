% Change the Iins as an array of numbers between 0 and 1
Iins = [0.2 0.5 0.6 .75 0.8];
visualize(Iins);

% Visualize the trajectory of some initial values of the FitzHugh-Nagamo Model
function visualize(Iins)
    for n = 1:length(Iins)
        Iin = Iins(n);
        
        % Sets up figure
        figure('Name',['Initial=' num2str(Iin)],'NumberTitle','off');
        y1Lim = [-5,+5]; y2Lim = [0,0.1];
        clear tmp;
        subplot(2,3,[1,2,4,5]);


        % Then, plots trajectory
        [tout,yout] = ode15s(@(t,y)FHode(y,Iin),(0:0.01:625),[-1.0;0.05]);
        subplot(2,3,[1,2,4,5]);
        hold on;
        plot(yout(1,1),yout(1,2),'ko');
        plot(yout(:,1),yout(:,2),'k-');
        hold off;
        xlim(y1Lim); xlabel('x variable');
        ylim(y2Lim); ylabel('y variable');
        
        % Then, plots the x- and y-value as two functions of t
        axis square;
        subplot(2,3,3);plot(tout,yout(:,1),'k-'); xlim([min(tout),max(tout)]); xlabel('time');ylabel('x variable');
        subplot(2,3,6);plot(tout,yout(:,2),'k-'); xlim([min(tout),max(tout)]); xlabel('time');ylabel('y variable');
    end
end

% Sets up the ODE for FitzHugh-Nagomo Model
function rhs = FHode(y,Iinput)
    v = y(1,:); w = y(2,:);
    rhs = [v - 2 - v.^3/3 - (40*w-7/2) + Iinput ; (v/4+7/16)/200-1/40*w];
end


