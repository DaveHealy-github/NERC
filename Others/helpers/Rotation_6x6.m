clear all; close all; clc;
syms Zn Zt Zn1 Zn3 Z12 Z13 Zt1 Zt2

C=[Zn1 Z12 Z13 0 0 0;
   Z12 Zn1 Z13 0 0 0;
   Z13 Z13 Zn3 0 0 0;
   0 0 0 Zt1 0 0;
   0 0 0 0 Zt1 0;
   0 0 0 0 0 Zt2; ];

A11=0; A12=0; A13=1; A21=0; A22=1; A23=0; A31=-1; A32=0; A33=0;  % rotated with respect to x2 axis (horizontal cracks become haveing x1 normal)
% A11=1; A12=0; A13=0; A21=0; A22=0; A23=1; A31=0; A32=-1; A33=0;  % rotated with respect to x1 axis (horizontal cracks become haveing x2 normal)
% A11=1; A12=0; A13=0; A21=0; A22=1; A23=0; A31=0; A32=0; A33=1; % no rotation

A=[A11 A12 A13;
    A21 A22 A23;
    A31 A32 A33;];

MA= [A11*A11 A12*A12 A13*A13 A12*A13 A11*A13 A11*A12;
     A21*A21 A22*A22 A23*A23 A22*A23 A21*A23 A21*A22;
     A31*A31 A32*A32 A33*A33 A32*A33 A31*A33 A31*A32;
     2*A21*A31 2*A22*A32 2*A23*A33 A22*A33+A23*A32 A21*A33+A23*A31 A21*A32+A22*A31;
     2*A11*A31 2*A12*A32 2*A13*A33 A12*A33+A13*A32 A11*A33+A13*A31 A11*A32+A12*A31;
     2*A11*A21 2*A12*A22 2*A13*A23 A12*A23+A13*A22 A11*A23+A13*A21 A11*A22+A12*A21;];

Crot=MA'*C*MA
