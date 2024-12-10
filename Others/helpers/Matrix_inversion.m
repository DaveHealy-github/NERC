clear all; clc; close all;
syms C11 C12 C23 C22 C44 C55

C_TI_x1 ...
= ...
[C11 C12 C12 0 0 0;
 C12 C22 C23 0 0 0;
 C12 C23 C22 0 0 0;
 0 0 0 C44 0 0;
 0 0 0 0 C55 0;
 0 0 0 0 0 C55;];

inv(C_TI_x1)