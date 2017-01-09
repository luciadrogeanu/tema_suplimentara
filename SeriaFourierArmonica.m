%semnal dreptunghiular x de perioada 40s, durata 9s, factor de umplere
%9/T*100
T=40
w0=(2*pi)/T
t=0:0.1:100
x=square(w0*t, 9/T*100)
N=50
%C=matricea coeficientilor pe care o egalam initial cu zero
C=zeros(size(N+1))
%x2 va fi semnalul reconstruit pentru SFE
x2=0
%x3 va fi semnalul reconstruit pentru SFA
x3=0
%matricea amplitudinilor
A=zeros(size(N+1))
%matricea fazelor
fi=zeros(size(N+1))
%for pt calculul coeficientilor de la 0 la 50
for k=0:N
    %integrala de la 0 la 9 unde semnalul este 1
    fun2=@(t) 1.*exp(-1j*k*w0*t)
    %integrala de la 9 la 40 unde semnalul este -1
    fun3=@(t) (-1).*exp(-1j*k*w0*t)
    %formula de calcul a coeficientilor (integrala pe o perioada
    C(k+1)=1/T*(integral(fun2, 0, 9)+integral(fun3, 9, 40))
    %semnalul va fi suma modulelor coeficientilor
    x2=x2+real(C(k+1)*exp(1j*k*w0*t))
    %amplitudinea k este dublul modulului coeficientului de ordin k
    A(k+1)=2*abs(C(k+1))
    %calculul fazei dupa formula
    if(2*real(C(k+1))<0)
        fi(k+1)=-atan((-2*imag(C(k+1)))/(2*real(C(k+1))))+pi
    else
        fi(k+1)=-atan((-2*imag(C(k+1)))/(2*real(C(k+1))))
    end
    %semnalul reconstruit
    x3=x3+A(k+1)*cos(k*w0*t+fi(k+1))
end
%componenta continua
fun1=@(t) 1*exp(-1j*0*w0*t)
C0=1/T*integral(fun1, 0, 9)
A0=2*C0
%la semnal adunam componenta continua (nu se poate calcula in for)
x3=x3+A0

subplot(3, 1, 1)
stem((0:N), A); title('Spectru de amplitudini');

subplot(3, 1, 2)
stem((0:N),fi); title('Spectru de faza');

subplot(3,1,3)
%reprezentarea semnalului initial
plot(t, x, '-')
hold on
%reprezentarea semnalului reconstruit
plot(t, x3, '.'); title('Semnalul reconstruit');
hold off

