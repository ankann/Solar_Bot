clc
lat_k=deg2rad(22.57);
lat_d = deg2rad(28.25);
long_k = deg2rad(88.364);
v_x=2;
v_y=2;
e=0.0167;
omega=30*10^3;
w=(108/60)*10^6;
D_p=3;
e0 = 23.4*2*pi/360;
D_w =4;
m = [0 0 0 0 1 0];
s = [0 0 0 0 0 1];

phi = long_k;
%for t=1:5
c=clock;
for t=0:10
long_k=long_k+v_x*t;
lat_k=lat_k+v_y*t;
c(6)=c(4)+mod(t,60);
c(5)=c(5)+floor(t/60);
c(4)=c(5)+floor(t/3600);
c(3)=c(5)+floor(t/(3600*24));
c(2)=c(5)+floor(t/(3600*24*30));
for i=1:t+1
    c=double(c);
    s1=c(4);
    c(4)=s1+(i-1);
    s1=0;
    lambda(i) = days(c(3),c(2),c(1),D_p)/365.25;
    delta(i) = asin(sin(lambda(i))*sin(e0));
    %delta(size(lambda))=asin(sin(lambda(size(lambda)))*sin(e0));
    t_ecc=2*e*w*cos(2*days(c(3),c(2),c(1),D_p)/365.25)/omega;
    t_obli=(1/4*omega)*e0^2*w*vpa(cos((4*pi*days(c(3),c(2),c(1),D_w))/365.25));
    s_c = c + m.*(4*(lat_k - lat_d)) + s.*(t_ecc +t_obli);
    v(i) = (s_c(4) + s_c(5)/60 + s_c(6)/3600 - 12)*15*2*pi/360;
    %disp(Z(1,1));
    syms x y z
    x = 0:0.2:2*pi;
    y = 0:0.1:pi;
    [X, Y] = meshgrid(x,y);
    %X= cat(3,X,X1);
    A = ones(numel(x),numel(y));
    Z=-1.*A -sin(Y).*sin(X).*(cos(delta(i))+1).*sin(v(i))+sin(Y).*cos(X).*sin(delta(i)).*cos(phi)+cos(Y).*(cos(delta(i))*cos(phi)*cos(v(i))+sin(delta(i))*sin(phi));
    % surfht(x,y,Z);
    %Z=cat(3,Z,[(1:t));
    b=[0 0];
    %disp(i);
    %X1=X;
    contour(X,Y,Z,b);
    hold on;
end
end

function d=days(day,month,yr,d_p)
    if(mod(yr,4)==0&&mod(yr,400)~=0)
        mo=[31 29 31 30 31 31 30 31 30 31 30 31];
    else
        mo=[31 28 31 30 31 31 30 31 30 31 30 31];
    end
    if(d_p==3)
        day1=0;
        for i=0:month
            day1=day1+mo(mod(i,12)+1);
        end
        d=day1+day-3;
    end
    if(d_p==4)
        day1=0;m1=0;
        for i=0:month
            day1=day1+mo(mod(i,12)+1);
        end
        for j=1:6
            m1=m1+mo(j);
        end
        d=abs(day1+day-4-m1);
    end
end
