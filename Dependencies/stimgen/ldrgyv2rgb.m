
function [rgb] = ldrgyv2rgb(ld,rg,yv)
% Qasim/Romain color solver
% Converts cartesian coordinates to calibrated colors for display
% import conversion matrix during the paradigm init call.
% See solver documentation for information on generating solver matrix

% Note: this function can handle multiple coordinates simultaneously to improve speed
% when generating large numbers of colors

global g_strctParadigm
ldrgyv = [ld; rg; yv];

%matrix=[	1 1 dRyv;
%1 dGrg dGyv;
%1 dBrg 1];


% Values for MIT, july 2015
matrix = g_strctParadigm.m_strctConversionMatrices.ldgyb;

%{
matrix=  [1	1	0.304455308556149
			1	-0.253409957412193	-0.257060304918859
			1	0.00799723127465672	1];
%}
% These are the correct values given out
% by the calib_correct function.


rgb = matrix*ldrgyv/2.0 + 0.5;
return;

