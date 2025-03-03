function Bpq=BpqxyErFun(x, y, p, q)
%Expansion in terms of the error function
global r s2 c2 xi tp t0 zeta zeta2
r= p+q; s2= q/r; c2= p/r; xi= x*y/2/r; tp= 1/y;
t0= t0proc(xi,c2);
zeta2= 2*(phit(s2,xi,tp)-phit(s2,xi,t0));
zeta= real(sqrt(zeta2));
if tp < t0
  zeta= -zeta;
end
if abs(zeta)>0.1
  Bpq=Bpqxy630sub1(x, y, p, q);
else
  Bpq=Bpqxy630sub2(x, y, p, q);
end
end

function Bpq=Bpqxy630sub1(x, y, p, q)
%Expansion for abs(zeta) not small
global r s2 c2 xi tp t0 zeta zeta2
gk=g2kcoeff(y, t0, c2, zeta);
Sr=gk(1)-gk(2)/r+gk(3)*3/r^2-gk(4)*15/r^3;
Reta=exp(-r*zeta2/2)*Sr/sqrt(2*pi*r);
Bpq=0.5*erfc(zeta*sqrt(r/2))+ Reta;
end

function Bpq=Bpqxy630sub2(x, y, p, q)
% Expansion for small abs(zeta)
global r s2 c2 xi tp t0 zeta zeta2
r2=r*r;
r3=r2*r;
gjk=gjkcoeff(t0,c2);
for k=0:3
  g= gjk(1,2*k+1);
  for j=1:6-2*k
    g= g+gjk(j+1,2*k+1)*zeta^j;
  end
  gk(k+1)= g;
end
Sr=gk(1)-gk(2)/r+gk(3)*3/r2-gk(4)*15/r3;
Reta=exp(-r*zeta2/2)*Sr/sqrt(2*pi*r);
Bpq=0.5*erfc(zeta*sqrt(r/2))+ Reta;
end

function phi=phit(s2,xi,t)
phi=log(t)-s2*log(t-1)+xi*t;
end

function t0=t0proc(xi,c2)
if xi >= c2
  t0= (xi-c2+sqrt((xi-c2)^2+4*xi))/2/xi;
else
  t0= 2/(c2-xi+sqrt((xi-c2)^2+4*xi));
end
end

function tk=tkcoeff(t0, c2)
%output tk, coefficients of the expansion.
% fk0(k) is the k-th derivative of phi(t) at t=t0
global t02 t03 t04 t05 t06 t07 t08 t09 t010 t011 t012
t0m1=t0-1;
t0m12=t0m1*t0m1;
t0m13=t0m12*t0m1;
t0m14=t0m13*t0m1;
t0m15=t0m14*t0m1;
t0m16=t0m15*t0m1;
t0m17=t0m16*t0m1;
t0m18=t0m17*t0m1;
t0m19=t0m18*t0m1;
t02=t0*t0;
t03=t02*t0;
t04=t03*t0;
t05=t04*t0;
t06=t05*t0;
t07=t06*t0;
t08=t07*t0;
t09=t08*t0;
t010=t09*t0;
t011=t010*t0;
t012=t011*t0;
fk0(2)= -(c2*t02-2*t0+1)/t02/t0m12;
fk0(3)= 2*(c2*t03-3*t02+3*t0-1)/t03/t0m13;
fk0(4)= -6*(c2*t04-4*t03+6*t02-4*t0+1)/t04/t0m14;
fk0(5)= 24*(c2*t05-5*t04+10*t03-10*t02+5*t0-1)/t05/t0m15;
fk0(6)= -120*(c2*t06-6*t05+15*t04-20*t03+15*t02-6*t0+1)/t06/t0m16;
fk0(7)= 720*(c2*t07-7*t06+21*t05-35*t04+35*t03-21*t02+7*t0-1)/t07/t0m17;
fk0(8)= -5040*(c2*t08-8*t07+28*t06-56*t05+70*t04-56*t03+28*t02-8*t0+1)/t08/t0m18;
fk0(9)= 40320*(c2*t09-9*t08+36*t07-84*t06+126*t05-126*t04+84*t03-...
    36*t02+9*t0-1)/t09/t0m19;
tk(1)= 1/sqrt(fk0(2));
tk(2)= -1/6*fk0(3)*tk(1)^2/fk0(2);
tk(3)= -1/72*tk(1)^3*(3*fk0(2)*fk0(4)-5*fk0(3)^2)/fk0(2)^2;
tk(4)= -1/1080*tk(1)^4*(9*fk0(2)^2*fk0(5)-45*fk0(2)*fk0(3)*fk0(4)+40*fk0(3)^3)/fk0(2)^3;
tk(5)= -1/17280*tk(1)^5*(24*fk0(2)^3*fk0(6)-168*fk0(2)^2*fk0(3)*fk0(5)-105*fk0(2)^2*fk0(4)^2+...
    630*fk0(2)*fk0(3)^2*fk0(4)-385*fk0(3)^4)/fk0(2)^4;
tk(6)= -1/136080*tk(1)^6*(27*fk0(2)^4*fk0(7)-252*fk0(2)^3*fk0(3)*fk0(6)-...
    378*fk0(2)^3*fk0(4)*fk0(5)+1260*fk0(2)^2*fk0(3)^2*fk0(5)+1575*fk0(2)^2*fk0(3)*fk0(4)^2-...
    4200*fk0(2)*fk0(3)^3*fk0(4)+1960*fk0(3)^5)/fk0(2)^5;
tk(7)= -1/43545600*tk(1)^7*(1080*fk0(2)^5*fk0(8)-12960*fk0(2)^4*fk0(3)*fk0(7)-...
    22680*fk0(2)^4*fk0(4)*fk0(6)-13608*fk0(2)^4*fk0(5)^2+83160*fk0(2)^3*fk0(3)^2*fk0(6)+...
    249480*fk0(2)^3*fk0(3)*fk0(4)*fk0(5)+51975*fk0(2)^3*fk0(4)^3-...
    360360*fk0(2)^2*fk0(3)^3*fk0(5)-675675*fk0(2)^2*fk0(3)^2*fk0(4)^2+...
    1126125*fk0(2)*fk0(3)^4*fk0(4)-425425*fk0(3)^6)/fk0(2)^6;
tk(8)= -1/3265920*tk(1)^8*(9*fk0(2)^6*fk0(9)-135*fk0(2)^5*fk0(3)*fk0(8)-...
    270*fk0(2)^5*fk0(4)*fk0(7)-378*fk0(2)^5*fk0(5)*fk0(6)+...
    1080*fk0(2)^4*fk0(3)^2*fk0(7)+3780*fk0(2)^4*fk0(3)*fk0(4)*fk0(6)+...
    2268*fk0(2)^4*fk0(3)*fk0(5)^2+2835*fk0(2)^4*fk0(4)^2*fk0(5)-...
    5880*fk0(2)^3*fk0(3)^3*fk0(6)-26460*fk0(2)^3*fk0(3)^2*fk0(4)*fk0(5)-...
    11025*fk0(2)^3*fk0(3)*fk0(4)^3+23520*fk0(2)^2*fk0(3)^4*fk0(5)+...
    58800*fk0(2)^2*fk0(3)^3*fk0(4)^2-70560*fk0(2)*fk0(3)^5*fk0(4)+...
    22400*fk0(3)^7)/fk0(2)^7;
end

function gk=g2kcoeff(y, t0, c2, zeta)
%Output the coefficents g2k 
global t02 t03 t04 t05 t06 t07 t08 t09 t010 t011 t012
zeta2=zeta*zeta;
zeta3=zeta*zeta2;
zeta5=zeta3*zeta2;
zeta7=zeta5*zeta2;
y2=y*y;
y3=y2*y;
y4=y3*y;
y5=y4*y;
y6=y5*y;
tk=tkcoeff(t0, c2);
gk(1)= -tk(1)/t0/(t0*y-1)-1/zeta;
gk(2)= -(3*t04*y2*tk(3)-6*t03*y2*tk(1)*tk(2)+3*t02*y2*tk(1)^3-...
    6*t03*y*tk(3)+9*t02*y*tk(1)*tk(2)-3*t0*y*tk(1)^3+...
    3*t02*tk(3)-3*t0*tk(1)*tk(2)+tk(1)^3)/t03/(t0*y-1)^3-1/zeta3;
gk(3)= -(5*t08*y4*tk(5)-10*t07*y4*tk(1)*tk(4)-10*t07*y4*tk(2)*tk(3)+...
    15*t06*y4*tk(1)^2*tk(3)+15*t06*y4*tk(1)*tk(2)^2-...
    20*t05*y4*tk(1)^3*tk(2)+5*t04*y4*tk(1)^5-20*t07*y3*tk(5)+...
    35*t06*y3*tk(1)*tk(4)+35*t06*y3*tk(2)*tk(3)-...
    45*t05*y3*tk(1)^2*tk(3)-45*t05*y3*tk(1)*tk(2)^2+...
    50*t04*y3*tk(1)^3*tk(2)-10*t03*y3*tk(1)^5+30*t06*y2*tk(5)-...
    45*t05*y2*tk(1)*tk(4)-45*t05*y2*tk(2)*tk(3)+50*t04*y2*tk(1)^2*tk(3)+...
    50*t04*y2*tk(1)*tk(2)^2-50*t03*y2*tk(1)^3*tk(2)+10*t02*y2*tk(1)^5-...
    20*t05*y*tk(5)+25*t04*y*tk(1)*tk(4)+25*t04*y*tk(2)*tk(3)-25*t03*y*tk(1)^2*tk(3)-...
    25*t03*y*tk(1)*tk(2)^2+25*t02*y*tk(1)^3*tk(2)-5*t0*y*tk(1)^5+5*t04*tk(5)-...
    5*t03*tk(1)*tk(4)-5*t03*tk(2)*tk(3)+5*t02*tk(1)^2*tk(3)+...
    5*t02*tk(1)*tk(2)^2-5*t0*tk(1)^3*tk(2)+tk(1)^5)/t05/(t0*y-1)^5-1/zeta5;
gk(4)= -(7*t012*y6*tk(7)-14*t011*y6*tk(1)*tk(6)-14*t011*y6*tk(2)*tk(5)-...
    14*t011*y6*tk(3)*tk(4)+21*t010*y6*tk(1)^2*tk(5)+...
    42*t010*y6*tk(1)*tk(2)*tk(4)+21*t010*y6*tk(1)*tk(3)^2+...
    21*t010*y6*tk(2)^2*tk(3)-28*t09*y6*tk(1)^3*tk(4)-84*t09*y6*tk(1)^2*tk(2)*tk(3)-...
    28*t09*y6*tk(1)*tk(2)^3+35*t08*y6*tk(1)^4*tk(3)+...
    70*t08*y6*tk(1)^3*tk(2)^2-42*t07*y6*tk(1)^5*tk(2)+...
    7*t06*y6*tk(1)^7-42*t011*y5*tk(7)+77*t010*y5*tk(1)*tk(6)+...
    77*t010*y5*tk(2)*tk(5)+77*t010*y5*tk(3)*tk(4)-...
    105*t09*y5*tk(1)^2*tk(5)-210*t09*y5*tk(1)*tk(2)*tk(4)-...
    105*t09*y5*tk(1)*tk(3)^2-105*t09*y5*tk(2)^2*tk(3)+...
    126*t08*y5*tk(1)^3*tk(4)+378*t08*y5*tk(1)^2*tk(2)*tk(3)+...
    126*t08*y5*tk(1)*tk(2)^3-140*t07*y5*tk(1)^4*tk(3)-...
    280*t07*y5*tk(1)^3*tk(2)^2+147*t06*y5*tk(1)^5*tk(2)-...
    21*t05*y5*tk(1)^7+105*t010*y4*tk(7)-175*t09*y4*tk(1)*tk(6)-...
    175*t09*y4*tk(2)*tk(5)-175*t09*y4*tk(3)*tk(4)+...
    217*t08*y4*tk(1)^2*tk(5)+434*t08*y4*tk(1)*tk(2)*tk(4)+...
    217*t08*y4*tk(1)*tk(3)^2+217*t08*y4*tk(2)^2*tk(3)-...
    238*t07*y4*tk(1)^3*tk(4)-714*t07*y4*tk(1)^2*tk(2)*tk(3)-...
    238*t07*y4*tk(1)*tk(2)^3+245*t06*y4*tk(1)^4*tk(3)+...
    490*t06*y4*tk(1)^3*tk(2)^2-245*t05*y4*tk(1)^5*tk(2)+...
    35*t04*y4*tk(1)^7-140*t09*y3*tk(7)+210*t08*y3*tk(1)*tk(6)+...
    210*t08*y3*tk(2)*tk(5)+210*t08*y3*tk(3)*tk(4)-...
    238*t07*y3*tk(1)^2*tk(5)-476*t07*y3*tk(1)*tk(2)*tk(4)-...
    238*t07*y3*tk(1)*tk(3)^2-238*t07*y3*tk(2)^2*tk(3)+...
    245*t06*y3*tk(1)^3*tk(4)+735*t06*y3*tk(1)^2*tk(2)*tk(3)+...
    245*t06*y3*tk(1)*tk(2)^3-245*t05*y3*tk(1)^4*tk(3)-...
    490*t05*y3*tk(1)^3*tk(2)^2+245*t04*y3*tk(1)^5*tk(2)-...
    35*t03*y3*tk(1)^7+105*t08*y2*tk(7)-140*t07*y2*tk(1)*tk(6)-...
    140*t07*y2*tk(2)*tk(5)-140*t07*y2*tk(3)*tk(4)+...
    147*t06*y2*tk(1)^2*tk(5)+294*t06*y2*tk(1)*tk(2)*tk(4)+...
    147*t06*y2*tk(1)*tk(3)^2+147*t06*y2*tk(2)^2*tk(3)-...
    147*t05*y2*tk(1)^3*tk(4)-441*t05*y2*tk(1)^2*tk(2)*tk(3)-...
    147*t05*y2*tk(1)*tk(2)^3+147*t04*y2*tk(1)^4*tk(3)+...
    294*t04*y2*tk(1)^3*tk(2)^2-147*t03*y2*tk(1)^5*tk(2)+...
    21*t02*y2*tk(1)^7-42*t07*y*tk(7)+49*t06*y*tk(1)*tk(6)+...
    49*t06*y*tk(2)*tk(5)+49*t06*y*tk(3)*tk(4)-49*t05*y*tk(1)^2*tk(5)-...
    98*t05*y*tk(1)*tk(2)*tk(4)-49*t05*y*tk(1)*tk(3)^2-...
    49*t05*y*tk(2)^2*tk(3)+49*t04*y*tk(1)^3*tk(4)+...
    147*t04*y*tk(1)^2*tk(2)*tk(3)+49*t04*y*tk(1)*tk(2)^3-...
    49*t03*y*tk(1)^4*tk(3)-98*t03*y*tk(1)^3*tk(2)^2+...
    49*t02*y*tk(1)^5*tk(2)-7*t0*y*tk(1)^7+7*t06*tk(7)-...
    7*t05*tk(1)*tk(6)-7*t05*tk(2)*tk(5)-7*t05*tk(3)*tk(4)+...
    7*t04*tk(1)^2*tk(5)+14*t04*tk(1)*tk(2)*tk(4)+7*t04*tk(1)*tk(3)^2+...
    7*t04*tk(2)^2*tk(3)-7*t03*tk(1)^3*tk(4)-21*t03*tk(1)^2*tk(2)*tk(3)-...
    7*t03*tk(1)*tk(2)^3+7*t02*tk(1)^4*tk(3)+14*t02*tk(1)^3*tk(2)^2-...
    7*t0*tk(1)^5*tk(2)+tk(1)^7)/t07/(t0*y-1)^7-1/zeta7;
end
function gjk=gjkcoeff(t0, c2)
%Output the coefficents gjk
global t02 t03 t04 t05 t06 t07 t08 t09 t010 t011 t012
tk=tkcoeff(t0, c2);
gjk(1,1)= -(t0*tk(2)-tk(1)^2)/t0/tk(1);
gjk(2,1)= -(tk(1)*tk(3)-tk(2)^2)/tk(1)^2;
gjk(3,1)= -(tk(1)^2*tk(4)-2*tk(1)*tk(2)*tk(3)+tk(2)^3)/tk(1)^3;
gjk(4,1)= -(tk(1)^3*tk(5)-2*tk(1)^2*tk(2)*tk(4)-tk(1)^2*tk(3)^2+...
    3*tk(1)*tk(2)^2*tk(3)-tk(2)^4)/tk(1)^4;
gjk(5,1)= -(tk(1)^4*tk(6)-2*tk(1)^3*tk(2)*tk(5)-2*tk(1)^3*tk(3)*tk(4)+...
    3*tk(1)^2*tk(2)^2*tk(4)+3*tk(1)^2*tk(2)*tk(3)^2-4*tk(1)*tk(2)^3*tk(3)+...
    tk(2)^5)/tk(1)^5;
gjk(6,1)= -(tk(1)^5*tk(7)-2*tk(1)^4*tk(2)*tk(6)-2*tk(1)^4*tk(3)*tk(5)-...
    tk(1)^4*tk(4)^2+3*tk(1)^3*tk(2)^2*tk(5)+6*tk(1)^3*tk(2)*tk(3)*tk(4)+...
    tk(1)^3*tk(3)^3-4*tk(1)^2*tk(2)^3*tk(4)-6*tk(1)^2*tk(2)^2*tk(3)^2+...
    5*tk(1)*tk(2)^4*tk(3)-tk(2)^6)/tk(1)^6;
gjk(7,1)= -(tk(1)^6*tk(8)-2*tk(1)^5*tk(2)*tk(7)-2*tk(1)^5*tk(3)*tk(6)-...
    2*tk(1)^5*tk(4)*tk(5)+3*tk(1)^4*tk(2)^2*tk(6)+6*tk(1)^4*tk(2)*tk(3)*tk(5)+...
    3*tk(1)^4*tk(2)*tk(4)^2+3*tk(1)^4*tk(3)^2*tk(4)-4*tk(1)^3*tk(2)^3*tk(5)-...
    12*tk(1)^3*tk(2)^2*tk(3)*tk(4)-4*tk(1)^3*tk(2)*tk(3)^3+5*tk(1)^2*tk(2)^4*tk(4)+...
    10*tk(1)^2*tk(2)^3*tk(3)^2-6*tk(1)*tk(2)^5*tk(3)+tk(2)^7)/tk(1)^7;
gjk(1,3)= -1/tk(1)^3/t03*(3*t03*tk(1)^2*tk(4)-3*t03*tk(1)*tk(2)*tk(3)+...
    t03*tk(2)^3-3*t02*tk(1)^3*tk(3)+3*t0*tk(1)^4*tk(2)-tk(1)^6);
gjk(2,3)= -3*(tk(1)^3*tk(5)-2*tk(1)^2*tk(2)*tk(4)-tk(1)^2*tk(3)^2+...
    3*tk(1)*tk(2)^2*tk(3)-tk(2)^4)/tk(1)^4;
gjk(3,3)= -3/tk(1)^5*(tk(1)^4*tk(6)-2*tk(1)^3*tk(2)*tk(5)-3*tk(1)^3*tk(3)*tk(4)+...
    4*tk(1)^2*tk(2)^2*tk(4)+5*tk(1)^2*tk(2)*tk(3)^2-7*tk(1)*tk(2)^3*tk(3)+...
    2*tk(2)^5);
gjk(4,3)= -1/tk(1)^6*(3*tk(1)^5*tk(7)-6*tk(1)^4*tk(2)*tk(6)-9*tk(1)^4*tk(3)*tk(5)-...
    6*tk(1)^4*tk(4)^2+12*tk(1)^3*tk(2)^2*tk(5)+36*tk(1)^3*tk(2)*tk(3)*tk(4)+...
    7*tk(1)^3*tk(3)^3-24*tk(1)^2*tk(2)^3*tk(4)-45*tk(1)^2*tk(2)^2*tk(3)^2+...
    42*tk(1)*tk(2)^4*tk(3)-10*tk(2)^6);
gjk(5,3)= -3/tk(1)^7*(tk(1)^6*tk(8)-2*tk(1)^5*tk(2)*tk(7)-3*tk(1)^5*tk(3)*tk(6)-...
    4*tk(1)^5*tk(4)*tk(5)+4*tk(1)^4*tk(2)^2*tk(6)+12*tk(1)^4*tk(2)*tk(3)*tk(5)+...
    7*tk(1)^4*tk(2)*tk(4)^2+8*tk(1)^4*tk(3)^2*tk(4)-8*tk(1)^3*tk(2)^3*tk(5)-...
    33*tk(1)^3*tk(2)^2*tk(3)*tk(4)-13*tk(1)^3*tk(2)*tk(3)^3+...
    15*tk(1)^2*tk(2)^4*tk(4)+36*tk(1)^2*tk(2)^3*tk(3)^2-25*tk(1)*tk(2)^5*tk(3)+5*tk(2)^7);
gjk(1,5)= -1/t05/tk(1)^5*(5*t05*tk(1)^4*tk(6)-5*t05*tk(1)^3*tk(2)*tk(5)-...
    5*t05*tk(1)^3*tk(3)*tk(4)+5*t05*tk(1)^2*tk(2)^2*tk(4)+5*t05*tk(1)^2*tk(2)*tk(3)^2-...
    5*t05*tk(1)*tk(2)^3*tk(3)+t05*tk(2)^5-5*t04*tk(1)^5*tk(5)+5*t03*tk(1)^6*tk(4)+...
    5*t03*tk(1)^5*tk(2)*tk(3)-5*t02*tk(1)^7*tk(3)-5*t02*tk(1)^6*tk(2)^2+...
    5*t0*tk(1)^8*tk(2)-tk(1)^10);
gjk(2,5)= -5/tk(1)^6*(tk(1)^5*tk(7)-2*tk(1)^4*tk(2)*tk(6)-2*tk(1)^4*tk(3)*tk(5)-...
    tk(1)^4*tk(4)^2+3*tk(1)^3*tk(2)^2*tk(5)+6*tk(1)^3*tk(2)*tk(3)*tk(4)+...
    tk(1)^3*tk(3)^3-4*tk(1)^2*tk(2)^3*tk(4)-6*tk(1)^2*tk(2)^2*tk(3)^2+...
    5*tk(1)*tk(2)^4*tk(3)-tk(2)^6);
gjk(3,5)= -5/tk(1)^7*(tk(1)^6*tk(8)-2*tk(1)^5*tk(2)*tk(7)-3*tk(1)^5*tk(3)*tk(6)-...
    3*tk(1)^5*tk(4)*tk(5)+4*tk(1)^4*tk(2)^2*tk(6)+10*tk(1)^4*tk(2)*tk(3)*tk(5)+...
    5*tk(1)^4*tk(2)*tk(4)^2+6*tk(1)^4*tk(3)^2*tk(4)-7*tk(1)^3*tk(2)^3*tk(5)-...
    24*tk(1)^3*tk(2)^2*tk(3)*tk(4)-9*tk(1)^3*tk(2)*tk(3)^3+11*tk(1)^2*tk(2)^4*tk(4)+...
    24*tk(1)^2*tk(2)^3*tk(3)^2-16*tk(1)*tk(2)^5*tk(3)+3*tk(2)^7);
gjk(1,7)= -1/tk(1)^7/t07*(7*t07*tk(1)^6*tk(8)-7*t07*tk(1)^5*tk(2)*tk(7)-...
    7*t07*tk(1)^5*tk(3)*tk(6)-7*t07*tk(1)^5*tk(4)*tk(5)+...
    7*t07*tk(1)^4*tk(2)^2*tk(6)+14*t07*tk(1)^4*tk(2)*tk(3)*tk(5)+...
    7*t07*tk(1)^4*tk(2)*tk(4)^2+7*t07*tk(1)^4*tk(3)^2*tk(4)-...
    7*t07*tk(1)^3*tk(2)^3*tk(5)-21*t07*tk(1)^3*tk(2)^2*tk(3)*tk(4)-...
    7*t07*tk(1)^3*tk(2)*tk(3)^3+7*t07*tk(1)^2*tk(2)^4*tk(4)+...
    14*t07*tk(1)^2*tk(2)^3*tk(3)^2-7*t07*tk(1)*tk(2)^5*tk(3)+...
    t07*tk(2)^7-7*t06*tk(1)^7*tk(7)+7*t05*tk(1)^8*tk(6)+...
    7*t05*tk(1)^7*tk(2)*tk(5)+7*t05*tk(1)^7*tk(3)*tk(4)-...
    7*t04*tk(1)^9*tk(5)-14*t04*tk(1)^8*tk(2)*tk(4)-...
    7*t04*tk(1)^8*tk(3)^2-7*t04*tk(1)^7*tk(2)^2*tk(3)+...
    7*t03*tk(1)^10*tk(4)+21*t03*tk(1)^9*tk(2)*tk(3)+7*t03*tk(1)^8*tk(2)^3-...
    7*t02*tk(1)^11*tk(3)-14*t02*tk(1)^10*tk(2)^2+7*t0*tk(1)^12*tk(2)-tk(1)^14);
end
