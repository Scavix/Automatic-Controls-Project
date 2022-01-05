%%animation script
T=t_step;
Y=y_step;
fhand=figure(20);
clf;
max_T=max(T);
min_Y=min(Y);
max_Y=max(Y);
Y_ss=Y(end);
axis([-10 10  0 5])

theta = rad2deg(W)
for i=1:10:1000
figure(20);
clf;
    axis([0 5 0 5])
    % Plot one frame...

    x_dron = [1 4 4 1]
    y_dron = [2 2 1 1]
    dron =  patch(x_dron, y_dron,'g');

    rotate(dron,[0 0 1], Y(i) * theta)
    hold on;
  
    title('Risposta al gradino W 1(t)')
 pause(0.1);
    if ishandle(fhand)==false; break;
    end
end