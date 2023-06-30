%dol=0.4
N=size(Particles,2);
for i=1:N % make all the particles around (0,0,0)
    Particles{1,i}.points(:,1:3)= Particles{1,i}.points(:,1:3)-mean(Particles{1,i}.points(:,1:3));
    max_nm(i)=max(max(Particles{1,i}.points(:,1:3)));
end
Cscale=round(max(max_nm));
for i=1:N % to make all the particles inside the unit box
    Particles{1,i}.points(:,1:3)= Particles{1,i}.points(:,1:3)./Cscale-mean(Particles{1,i}.points(:,1:3)./Cscale);
    Particles{1,i}.sigma(:,1:2)=   Particles{1,i}.sigma(:,1:2)./Cscale;
end
minLo=10;
maxLo=10000;%maximumlocalization in one particle
%Prepare the input for the JRMPC
[~,V1,meanlocs,~,~] = ProcessedExperimentalDataPraticle(Particles,1,N,minLo,maxLo);
[Particles1,V2,~,~,~] = ProcessedExperimentalDataPraticle(Particles,Cscale,N,minLo,maxLo);

M=size(V1,1);
raster_spacing=40;       % expected resolution [nm]