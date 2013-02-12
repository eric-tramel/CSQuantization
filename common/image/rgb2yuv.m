function [Y U V] = rgb2yuv(RGBFile)
WR = 0.299;
WB = 0.144;
WG = 0.587;
UMAX = 0.436;
VMAX = 0.615;

R = RGBFile(:,:,1);
G = RGBFile(:,:,2);
B = RGBFile(:,:,3);


Y = WR.*R + WG.*G + WB.*B;

U = UMAX .* ((B - Y)./(1 - WB));

V = VMAX .* ((R - Y)./(1 - WR));

