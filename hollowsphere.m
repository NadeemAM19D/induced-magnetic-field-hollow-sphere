clear all
clc
model = createpde;
r1=0.2;r2=0.5;r3=0.1;
R1 = [3,4,0,r2,r2,0,r2,r2,-r2,-r2]';
C1 = [1,0,0,r1]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
C2 = [1,0,0,r2]';
C2 = [C2;zeros(length(R1) - length(C2),1)];
C3 = [1,0,0,r3]';
C3 = [C3;zeros(length(R1) - length(C3),1)];
gm = [R1,C1,C2,C3];
sf = '(C1+C2+C3)*R1';
ns = char('R1','C1','C2','C3');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
pdegplot(model,"EdgeLabels","on","FaceLabels","on")
axis equal
%figure
mesh = generateMesh(model,GeometricOrder="linear",Hedge={[6,7,10,11],0.001});
pdeplot(mesh)
xc=mesh.Nodes(1,:);
yc=mesh.Nodes(2,:);
connectivity=mesh.Elements';
Nf1=findNodes(mesh,"region","Edge",[6,7]);
Nf2=findNodes(mesh,"region","Edge",[10,11]);
hold on
plot(mesh.Nodes(1,Nf1),mesh.Nodes(2,Nf1),"ok","MarkerFaceColor","g")
plot(mesh.Nodes(1,Nf2),mesh.Nodes(2,Nf2),"ok","MarkerFaceColor","g")
Ne1=findNodes(mesh,"region","Edge",[8,9]);
hold on
plot(mesh.Nodes(1,Ne1),mesh.Nodes(2,Ne1),"ok","MarkerFaceColor","g")
ER=findElements(mesh,"region","Face",2);
pdemesh(mesh.Nodes,mesh.Elements(:,ER),"EdgeColor","green")
muo=1;mur=10;kglobal=zeros(size(xc,2));fglobal=zeros(size(xc,2),1);
for i=1:size(connectivity,1)
ele=connectivity(i,:);
xele=[xc(ele(1,1)),xc(ele(1,2)),xc(ele(1,3))];
yele=[yc(ele(1,1)),yc(ele(1,2)),yc(ele(1,3))];
J=[xele(1)-xele(3),yele(1)-yele(3);xele(2)-xele(3),yele(2)-yele(3)];
Bp=inv(J)*[1,0,-1;0,1,-1];
rc=(xele(1)+xele(2)+xele(3))/3;
if ismember(i,ER)==1
kele=mur*rc*(Bp'*Bp)*0.5*det(J);
else
kele=muo*rc*(Bp'*Bp)*0.5*det(J);
end
% for i1=1:3
% for j1=1:3
% kglobal(ele(i1),ele(j1))=kglobal(ele(i1),ele(j1))+kele(i1,j1);
% end
% end
kglobal(ele,ele)=kglobal(ele,ele)+kele;
end
xb=mesh.Nodes(1,Nf2);yb=mesh.Nodes(2,Nf2);ha=100;
xb(size(xb,2)+1)=xb(1);yb(size(xb,2)+1)=yb(1);
for i=1:size(Nf2,2)-1
xele=[xc(Nf2(i)),xc(Nf2(i+1))];
yele=[yc(Nf2(i)),yc(Nf2(i+1))];
xm=(xele(1)+xele(2))/2;
ym=(yele(1)+yele(2))/2;
le=sqrt((xele(2)-xele(1))^2+(yele(2)-yele(1))^2);
theta=atan(abs(ym)/abs(xm));
aa=(xele(1)+xele(2))/2;bb=(xele(2)-xele(1))/2;
if xm>0 && ym<0
fglobal(Nf2(i),1)=fglobal(Nf2(i),1)-sin(theta)*(le/2)*(1-mur)*ha*(aa-bb/3);
fglobal(Nf2(i+1),1)=fglobal(Nf2(i+1),1)-sin(theta)*(le/2)*(1-mur)*ha*(aa+bb/3);
else 
fglobal(Nf2(i),1)=fglobal(Nf2(i),1)+sin(theta)*(le/2)*(1-mur)*ha*(aa-bb/3);
fglobal(Nf2(i+1),1)=fglobal(Nf2(i+1),1)+sin(theta)*(le/2)*(1-mur)*ha*(aa+bb/3);
end
end
xb=mesh.Nodes(1,Nf1);yb=mesh.Nodes(2,Nf1);ha=100;
xb(size(xb,2)+1)=xb(1);yb(size(xb,2)+1)=yb(1);
for i=1:size(Nf1,2)-1
xele=[xc(Nf1(i)),xc(Nf1(i+1))];
yele=[yc(Nf1(i)),yc(Nf1(i+1))];
xm=(xele(1)+xele(2))/2;
ym=(yele(1)+yele(2))/2;
le=sqrt((xele(2)-xele(1))^2+(yele(2)-yele(1))^2);
theta=atan(abs(ym)/abs(xm));
aa=(xele(1)+xele(2))/2;bb=(xele(2)-xele(1))/2;
if xm>0 && ym<0
fglobal(Nf1(i),1)=fglobal(Nf1(i),1)-sin(theta)*(le/2)*(mur-1)*ha*(aa-bb/3);
fglobal(Nf1(i+1),1)=fglobal(Nf1(i+1),1)-sin(theta)*(le/2)*(mur-1)*ha*(aa+bb/3);
else 
fglobal(Nf1(i),1)=fglobal(Nf1(i),1)+sin(theta)*(le/2)*(mur-1)*ha*(aa-bb/3);
fglobal(Nf1(i+1),1)=fglobal(Nf1(i+1),1)+sin(theta)*(le/2)*(mur-1)*ha*(aa+bb/3);
end
end
kmod=kglobal;fmod=fglobal;
Ne1=sort(Ne1,'descend');
for i=1:size(Ne1,2)
kmod(Ne1(i),:)=[];kmod(:,Ne1(i))=[];
fmod(Ne1(i),:)=[];
end
phim=kmod\fmod;flag=1;
for i=1:size(xc,2)
if ismember(i,Ne1)==1
phi(i)=0;
else
phi(i)=phim(flag);
flag=flag+1;
end
end
for i=1:size(connectivity,1)
ele=connectivity(i,:);
xele=[xc(ele(1,1)),xc(ele(1,2)),xc(ele(1,3))];
yele=[yc(ele(1,1)),yc(ele(1,2)),yc(ele(1,3))];
J=[xele(1)-xele(3),yele(1)-yele(3);xele(2)-xele(3),yele(2)-yele(3)];
Bp=inv(J)*[1,0,-1;0,1,-1];
xg(i)=(xele(1,1)+xele(1,2)+xele(1,3))/3;
yg(i)=(yele(1,1)+yele(1,2)+yele(1,3))/3;
phie=[phi(ele(1,1)),phi(ele(1,2)),phi(ele(1,3))]';
H=Bp*phie;
hx(i)=-H(1,1);
hy(i)=ha-H(2,1);
end
for i=1:size(xc,2)
[val, idx] = find(connectivity==i);
if ismember(i,Ne1)==1
Hy(i)=ha;
Hx(i)=0;
else
Hy(i)=sum(hy(val))/(size(val,1));
Hx(i)=sum(hx(val))/(size(val,1));
end
end
figure
pdeplot(mesh,XYData=Hy)
hold on
axis equal
colormap jet
%pdeplot(mesh)
figure
pdeplot(mesh,XYData=Hx)
hold on
axis equal
colormap jet
%pdeplot(mesh)























































































































