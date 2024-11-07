function [I,D,P,R,V] = GetTrlInfo(name)

% detection orientation
num_idx = [strfind(name,'Det_') strfind(name,'_Ima')];
D = str2double(name(num_idx(1)+4:num_idx(2)-1));

% imagery orientation
num_idx = [strfind(name,'Ima_') strfind(name,'_P')];
I = str2double(name(num_idx(1)+4:num_idx(2)-1));

% stimulus present or absent
num_idx = strfind(name,'P_');
P = str2double(name(num_idx+2));

% detection response
num_idx = strfind(name,'R_');
R = str2double(name(num_idx+2));

% imagery vividness rating
num_idx = strfind(name,'V');
V = str2double(name(num_idx+1));