%%% plots the shape functions used in the operator aproximation
addpath('..\..\matlabLib');
x=linspace(0,1,100);

% 1stOrder "hat functions"
y1st1 = (1-x);
y1st2 = [x,(1-x)] ;
y1st3 = x;
figure; plot(x,y1st1,[x,1+x],y1st2,1+x, y1st3, 'LineWidth', 1.5); axis('tight');
v=axis;
axis([v(1) v(2) v(3) 1.1]);
printEPS('','1st')


% 2ndOrder Lagrange
x=linspace(0,1,100);
y2nd1 = (1-x).*(2*(1-x)-1);
y2nd2 = 2*(1-x).*2.*(x);
y2nd3 = x.*(2*x-1);
x=linspace(0,2,100);
figure; plot(x,y2nd1,x,y2nd2,x,y2nd3, 'LineWidth', 1.5); axis('tight');
v=axis;
axis([v(1) v(2) v(3) 1.1]);
printEPS('','2nd')

% 3rdOrder Lagrange
x=linspace(0,1,100);
y3rd1 = (1-x).*(2-3*x).*(1-3*x)./2;
y3rd2 = 9*x.*(1-x).*(2-3*x)./2;
y3rd3 = 9*x.*(1-x).*(3*x-1)./2;
y3rd4 = x.*(2-3*x).*(1-3*x)./2;
x=linspace(0,3,100);
figure; plot(x,y3rd1,x,y3rd2,x,y3rd3,x,y3rd4, 'LineWidth', 1.5); axis('tight');
v=axis;
axis([v(1) v(2) v(3) 1.1]);
printEPS('','3rd')

%% PWS
x = linspace(0,1,100);
k=pi/2;
xi = (x(length(x))-x)/(x(length(x))-x(1));
yPWS1 = sin(k*(1-xi))./sin(k*(xi(1)));
yPWS2 = sin(k*(xi))./sin(k*(1-xi(length(xi))));
figure; plot(x,yPWS1,x,yPWS2, 'LineWidth', 1.5); axis('tight');
v=axis;
axis([v(1) v(2) v(3) 1.1]);
printEPS('','PWS')

%%
% x=linspace(1,10,100);
% 
% xi = (x(length(x))-x)/(x(length(x))-x(1));
% 
% y2nd1 = (1-xi).*(2*(1-xi)-1);
% y2nd2 = 2*(1-xi).*2.*(xi);
% y2nd3 = xi.*(2*xi-1);
% % x=linspace(0,2,100);
% figure; plot(x,y2nd1,x,y2nd2,x,y2nd3, 'LineWidth', 1.5); axis('tight');
% v=axis;
% axis([v(1) v(2) v(3) 1.1]);
% % printEPS('','2nd')



