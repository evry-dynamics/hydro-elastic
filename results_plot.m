% Matlab results of the plate (Solid part)
W_s=load('bg_s_40_40_W.mat').W;
figure
for i = 1:n
    p1=plot([0:60],W_s(i,:)*ScaleFacor,'k-','LineWidth',1.5);%
    hold on
end

% Matlab results of the plate with one side fluid (FSI)
W(4:n,1) = W(4:n,2);
W(4:n,61) = W(4:n,60);
for i = 1:n
    p2=plot([0:60],W(i,:)*ScaleFacor,'LineWidth',1.5); %,'b-'
    hold on
end


% Comsol results - solid element
load("comsol_solid.txt")
Comsol = reshape(comsol_solid(:,2),31,n);
for i = 1:n
    p2=plot([0:2:60],Comsol(:,i),'r:','LineWidth',1.5);
    hold on
end

xlim([0 60])
ylim([0 60])
xticks([0 60/3 60/3*2 60])
xticklabels({'\Gamma ','X','M','\Gamma'})
ylabel('Frequency (kHz)','FontName','Times New Roman','FontSize',15)
% legend([p1,p2],"Solid","FSI",'FontName','Times New Roman','FontSize',15)
legend([p1,p2],"Matlab-Plate","Comsol-Solid",'FontName','Times New Roman','FontSize',15)
axis on
box on
grid on