function  Yh_phase  = phase_test()
%PHASE_TEST Summary of this function goes here
%   Detailed explanation goes here
% This file tests the phase change when there is a shift of impulse
% in the 1-D input

Yh_phase=[];
for i=1:8
    X=zeros(32,1); 
    X(7+i)=1;
    [Yl,Yh] = dtwavexfm(X,5,'near_sym_b','qshift_b');
    Yh_phase(:,i)=phase(Yh{4});
%    Yh_magnitude=abs(Yh{3});
%    Yh_energy=sum(Yh_magnitude.^2)
end
Yh_increment=[];
% for i=2:32;
% Yh_increment(:,i-1)=Yh_phase(:,i)-Yh_phase(:,i-1)
% end
return;