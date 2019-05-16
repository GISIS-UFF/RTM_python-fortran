close all
clear all
close all

Ndamper= 100;
center= 350;

wx=50;
Nx = 383;
Nz = 141;


i = [1:Ndamper];
x = (i-1)*pi/Ndamper;
damper = (cos(x)/2+0.5);

mute = zeros(1,Nx);
taper_adjust = zeros(1,2*Nx);

plato_ini=round(Nx/2)+center-round(wx/2);
plato_end=round(Nx/2)+center+round(wx/2);

taper_adjust(plato_ini:plato_end) = 1;
taper_adjust(plato_end+1:plato_end+Ndamper) = damper;
taper_adjust(plato_ini:-1:plato_ini-Ndamper+1) = damper;

figure(1)
subplot(2,1,2)
plot(taper_adjust);

mute = taper_adjust(round(Nx/2)+1:3*round(Nx/2)-1);
subplot(2,1,1)
plot(mute)

%%
migrated = rand(Nz,Nx);
migrated_taper = rand(Nz,Nx);

for i=1:Nz
    migrated_taper(i,:) = migrated(i,:) .* mute;
end
figure(2)
subplot(2,1,1)
imagesc(migrated_taper)
subplot(2,1,2)
imagesc(migrated)





