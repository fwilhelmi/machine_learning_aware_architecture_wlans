% BB - 2019


function [S,D,Utilization,AirTime]=WIFIModel(lambda,L,DBPS,Rates,f,WIFI_std)

% f = 2.4E9 (802.11n / 2.4 GHz), 
% f = 5E9 (802.11ac / 5 GHz)

more off

MaxIter=1000;
ConvergenceInIter=0;

N=length(lambda); % Number of Active Nodes

CWmin=16;
mBO=5;
T0 = 9E-6;
K=100;

%lambda=50.*ones(1,N);
%L=10000;
pe=0.1.*ones(1,N);

% Variables
tau=zeros(N,MaxIter);
pc=zeros(N,MaxIter);
pf=zeros(N,MaxIter);
EBO=zeros(N,MaxIter);
p0=zeros(N,MaxIter);
p1=zeros(N,MaxIter);
p1m=zeros(N,MaxIter);
delta=zeros(N,MaxIter);
ENT=zeros(N,MaxIter);
EDs=zeros(N,MaxIter);
u=zeros(N,MaxIter);
a=zeros(N,MaxIter);
Pb=zeros(N,MaxIter);
ED2s=zeros(N,MaxIter);
ENq=zeros(N,MaxIter);
EDq=zeros(N,MaxIter);
ED=zeros(N,MaxIter);
ES=zeros(N,MaxIter);
AT=zeros(N,MaxIter);

% Initizalization

for n=1:N
   %lambda(n)=20*n;
%------------------------------------%    
%%OPCI� A
%    T(n)=L/Rates(n);
%    Tc(n)=0.1E-3;

%%OPCI� B
   if (WIFI_std == 0)   %IEEE 802.11n (2.4 GHz) + IEEE 802.11ac (5 GHz)
       if(f == 2.4E9)
           [T(n),Tc(n)] = TransmissionTimeRTSCTS11n(L,DBPS(n),min(DBPS));
       else
           [T(n),Tc(n)] = TransmissionTimeRTSCTS11ac(L,DBPS(n),min(DBPS));
       end
   else                 %IEEE 802.11ax (2.4 GHz) + IEEE 802.11ax (5 GHz)
       [T(n),Tc(n)] = TransmissionTimeRTSCTS11ax(1,L,20,1,1,1/2);
   end
%------------------------------------%   
   EBO(n,1)=(CWmin-1)/2;
   a(n,1)= lambda(n) * T(n);
   Pb(n,1)=0.1;
end


for i=1:MaxIter

  %disp('------------------------');
  %disp(i);
  %disp('------------------------');

  % Transmission Probability

  for n=1:N
    tau(n,i)=min(1,a(n,i)*(1-Pb(n,i)))/(EBO(n,i)+1);
  
    % For improving convergence

    if(i>8)
      tau(n,i)=0.25*tau(n,i)+0.25*tau(n,i-1)+0.25*tau(n,i-2)+0.25*tau(n,i-3);
    end
    

  end

  %disp(tau(:,i)');

  % Collision Probability

  for n=1:N
    pc(n,i)=1;
    for m=1:N
      if(m~=n)
        pc(n,i)=pc(n,i)*(1-tau(m,i));
      end
    end
    pc(n,i)=1-pc(n,i);

    if(i>8)
      pc(n,i)=0.25*pc(n,i)+0.25*pc(n,i-1)+0.25*pc(n,i-2)+0.25*pc(n,i-3);
    end
    

  end

  %disp(pc(:,i)');

  % Failure Probability

  for n=1:N
    pf(n,i)=pc(n,i)+(1-pc(n,i))*pe(n);
  end
  
  %disp(pf(:,i)');


  % Average number of transmissions

  for n=1:N
    ENT(n,i) = 1/(1-pf(n,i));
  end
  
  %disp(ENT(:,i));

  % Average Backoff (Max ret = infinite)
  
  for n=1:N
    EBO(n,i+1)=((1-pf(n,i)-pf(n,i)*(2*pf(n,i))^mBO)/(1-2*pf(n,i)))*((CWmin)/2)-(1/2);    
    %for b=0:100
    %  EBO(n,i+1)=EBO(n,i+1) + pf(n,i)^(b)*(1-pf(n,i))*(CWmin*2^min(mBO,b)-1)/2;
    %  disp([(CWmin.*2^min(mBO,b-1)-1)/2 EBO(n,i+1) auxEBO(n,i+1)]);  
    %  pause
    %end
  end

  % Empty slot probability

  for n=1:N
    p0(n,i)=1;
    for m=1:N
      if(m~=n)
        p0(n,i)=p0(n,i)*(1-tau(m,i));
      end
    end
  end

  % Successful slot probability

  for n=1:N
    p1(n,i)=0;
    for m=1:N
      if(m~=n)
        aux=1;
        for k=1:N
          if(k~=m && k~=n)
            aux=aux*(1-tau(k,i));  
          end
        end
        p1(n,i)=p1(n,i)+tau(m,i)*aux;
      end
    end
  end

  % Collision slot probability

  for n=1:N
    p1m(n,i)=1-p0(n,i)-p1(n,i);
  end

  %disp('BO probabilities');
  %disp(p0(:,i)');
  %disp(p1(:,i)');
  %disp(p1m(:,i)');

  % Average duration of 1 transmission

  ET1=zeros(1,N);

  for n=1:N
    for m=1:N
      if(m~=n)
        aux=1;
        for k=1:N
          if(k~=m && k~=n)
            aux=aux*(1-tau(k,i));  
          end
        end
        ET1(n)=ET1(n)+tau(m,i)*aux*T(m);
      end
    end
  end

  % Average slot duration

  for n=1:N
    delta(n,i) = p0(n,i)*T0 + ET1(n) + p1m(n,i)*Tc(n);
  end


  % Average Service Time

  for n=1:N
    EDs(n,i)=ENT(n,i)*(EBO(n,i)*delta(n,i)+T(n));
  end

  % Traffic load

  for n=1:N
   %a(n,i+1)=min(0.99,lambda(n)*EDs(n,i));
    a(n,i+1)=lambda(n)*EDs(n,i);

    if(i>8)
      a(n,i)=0.25*a(n,i)+0.25*a(n,i-1)+0.25*a(n,i-2)+0.25*a(n,i-3);
    end

  end


  % Blocking Probability

  for n=1:N
    if (a(n,i+1)== 1)
      Pb(n,i+1) = 1 /(K+1);
    else  
      Pb(n,i+1)=a(n,i+1)^K*(1-a(n,i+1))/(1-a(n,i+1)^(K+1));
    end
    u(n,i+1)=lambda(n)*EDs(n,i)*(1-Pb(n,i+1));
    %AT(n,i)=lambda(n)*(1-Pb(n,i+1))*ENT(n,i)*T(n);
  end

  % Queueing and Response Delay

  for n=1:N
    if (a(n,i+1)== 1)
      ENq(n,i)=(K)/2;
    else  
      ENq(n,i)=(a(n,i+1)/(1-a(n,i+1)))-((K+1)*a(n,i+1)^(K+1))/(1-a(n,i+1)^(K+1));
    end
    EDq(n,i)=ENq(n,i)/(lambda(n)*(1-Pb(n,i+1)));
    ED(n,i)=EDs(n,i)+EDq(n,i);
  end

  %disp('Queueing Delay and Blocking');
  %disp(Pb(:,i)');
  %disp(ENq(:,i)');
  %disp(EDq(:,i)');

  %disp('EDs');
  %disp(EBO(:,i)');
  %disp(ET1);
  %disp(delta(:,i)');
  %disp(a(:,i)');
  %disp(EDs(:,i)');

  % AirTime

%  for n=1:N
%      AT(n,i)=lambda(n)*(1-Pb(n,i+1))*ENT(n,i)*T(n);
%  end


  % Throughput

  for n=1:N
      ES(n,i)=lambda(n)*L*(1-Pb(n,i+1));
      % Without removing shared airtime 
      % AT(n,i)=lambda(n)*(1-Pb(n,i+1))*((ENT(N,i)-1)*T(n)+T(n));
      % Removing shared airtime (half for each, assuming only 2 nodes
      % collide)
      AT(n,i)=lambda(n)*(1-Pb(n,i+1))*((ENT(N,i)-1)*T(n)/2+T(n));
      %       if isnan(AT(n,i)), AT(n,i) = 0; end
  end
  
  % Convergence condition

  if(i>20)
    if(ceil((pc(1,i))*1E9)==ceil((pc(1,i-1))*1E9))
      %disp('convergence');
      %disp(i);
      ConvergenceInIter=i;
      break;
    %else
      %disp([ceil(pc(1,i)*1E12) ceil(pc(1,i-1)*1E12)]);
      %pause
    end
  end

  %disp('#############');
  %pause

end

if(i==MaxIter)
    ConvergenceInIter = MaxIter;
end

S=ES(:,ConvergenceInIter);
D=ED(:,ConvergenceInIter);
Utilization=u(:,ConvergenceInIter);
AirTime = AT(:,ConvergenceInIter);
sum(AirTime);

if isnan(S), S = 0; end

% for i = 1 : size(AirTime,2)
%     if isnan(AirTime(i)), AirTime(i) = 0; end
% end

%figure
%plot(tau(1,:));
%
%figure
%bar([(lambda.*L) ; ES(:,i)']');
%xlabel('Node');
%ylabel('Throughput');
%
%figure
%bar(EDs(:,i));
%
%figure
%bar(tau(:,i));
%
%
%figure
%bar(ENq(:,i));
%
%
%figure
%bar(EDq(:,i));

end



