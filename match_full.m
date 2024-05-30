%addpath('/home/nico/src/casadi') %add your path to casadi directory
clear all
close all

Zc1 = 16.4-139.5i %Monza R6-P
Rmod = 50; %consider a 50 ohm as modulation resistance (approximation)
Zc2 =  (Zc1*Rmod)/(Zc1+Rmod)
Sc = 10^(-20/10)/1000; %Chip sensitiviy

Pt = 10^(30/10)/1000; %Reader transmitted power
G=10^(6/10); %Reader gain
f0 = 915e6; %Operating frequency
lambda = 3e8/f0;
Pr = 10^(-75/10)/1000; %Reader sensitivity

%Define some constant to shorten the code
R1 = real(Zc1);
R2 = real(Zc2);
X1 = imag(Zc1);
X2 = imag(Zc2);

%Conjugate Matching
Zaf = real(Zc1) - 1j*imag(Zc1)

%Differential Matching
%This method is the one described in the IEEE RFID conference paper (solving a system of 2 equations)
%This method is an approximation
Xar0 = -(X1 + X2)/2; %initial condition (see (8) in conf paper)
Rar0 = ((R1^2+(X1+Xar0)^2)*(R2^2+(X2+Xar0)^2))^(1/4); %initial condition (see (9) in conf paper)
Zarp = complex(Rar0, Xar0);
Ro = [Rar0 Xar0];
fun = @(x)paramfun(x, Zc1, Zc2);
x = fsolve(fun, Ro);
Zar = complex(x(1),x(2));

%Differential Matching using Mayer's equations
%https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiKsNa4n7b_AhUg7rsIHYoBBR8QFnoECA8QAQ&url=https%3A%2F%2Fpublik.tuwien.ac.at%2Ffiles%2FPubDat_176539.pdf&usg=AOvVaw1gb7M6fcai0rP4xV8O64Wy
%This method is exact
Rar = sqrt(R1*R2*((R1+R2)^2+(X1-X2)^2)/(R1+R2)^2);
Xar = -(R2*X1+R1*X2)/(R1+R2);
Zar = complex(Rar, Xar)

%Computation of various constants
K = 16*Pt*Pr*real(Zc1)^2/(Sc^2*abs(Zc2-Zc1)^2)
Ro = -(R1-K*R2)/(1-K);
Xo = -(X1-K*X2)/(1-K);
centero = complex(Ro,Xo); %coord of the center
ro = sqrt(abs(((R1^2+X1^2-K*(R2^2+X2^2))/(1-K))-((R1-K*R2)/(1-K))^2-((X1-K*X2)/(1-K))^2)); %radius of the center

if ((abs(Zaf-centero) < ro) && (abs(Zar-centero) < ro)) %both inside the circle
  if (K>1)
    disp('Friis limited K>1')
    Zao = Zaf;
  else
    disp('Radar limited K<1')
    Zao = Zar;
  end
elseif ((abs(Zaf-centero) > ro) && (abs(Zar-centero) > ro)) %both outside the circle
  if (K>1)
    disp('Radar limited K>1')
    Zao = Zar;
  else
    disp('Friis limited K<1')
    Zao = Zaf;
  end
else
  %display('Optimization') %These lines require the installation of Casadi
  %opti = casadi.Opti();
  %Rao = opti.variable();
  %Xao = opti.variable();

  %opti.minimize(  (-4*R1*Rao/((R1+Rao)**2+(X1+Xao)**2))   );
  %opti.subject_to( ((R1+Rao)**2+(X1+Xao)**2)/((R2+Rao)**2+(X2+Xao)**2)-K==0 );
  %opti.subject_to( Rao>=0 );

  %p_opts = struct('print_time', 0); %quiet plugin
  %s_opts = struct('print_level', 0); %quiet solver
  %opti.solver('ipopt', p_opts, s_opts);

  %sol = opti.solve();

  %Rao = sol.value(Rao);
  %Xao = sol.value(Xao);
  %Zao = complex(Rao, Xao)

  display('Analytic'); %Parametrization of the constraint
  a = -2*ro*Ro*(Xo+X1);
  b = -ro*((Ro+R1)^2+ro^2+(Xo+X1)^2-2*Ro*(R1+Ro));
  c = 2*ro^2*(Xo+X1);
  if (c^2 > a^2+b^2)
    disp('No solutions') %should never happen
    return
  end
  theta1 = +acos(c/(sqrt(a^2+b^2)))+atan2(b,a);
  theta2 = -acos(c/(sqrt(a^2+b^2)))+atan2(b,a);
  Rao1 = Ro + ro*cos(theta1);
  Xao1 = Xo + ro*sin(theta1);
  Zao1 = complex(Rao1, Xao1);
  Rao2 = Ro + ro*cos(theta2);
  Xao2 = Xo + ro*sin(theta2);
  Zao2 = complex(Rao2, Xao2);
  %Normally we should check Zao1 and Zao2 but optimal solution is always Zao1
  Zao = Zao1; 
end
Zao

phi = linspace(0,2*pi, 100);
plot(centero+ro*exp(1j*phi),'m')
hold on
plot(real(Zc1), imag(Zc1),'go')
plot(real(Zc2), imag(Zc2),'co')
plot(real(Zaf), imag(Zaf),'bo')
plot(real(Zar), imag(Zar),'ro')
plot(real(Zao), imag(Zao),'ko')
hold on
xlim([-200 200])
ylim([-200 200])
grid on
