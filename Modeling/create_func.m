clc; clear; close all;

syms x y
A = [x^2 + sin(y);cos(x)]; % Example symbolic expression

% Convert symbolic expression to a string
A_str = char(A); 

% Create function file
fid = fopen('getA.m', 'w');
fprintf(fid, 'function A = getA()\n');
fprintf(fid, '\tA = %s;\n', A_str);
fprintf(fid, 'end\n');
fclose(fid);

disp('Function getA.m created successfully.');

% or 

matlabFunction(A, 'File', 'getA_alt');

