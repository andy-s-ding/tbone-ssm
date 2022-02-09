function [Ox,Oy,Oz,sx,sy,sz] = extractMetaFromNrrd(meta)
temp = meta.spaceorigin;
temp(temp=='(')=[];
temp(temp==')')=[];
temp(temp==',')=' ';
[Ox,Oy,Oz] = strread(temp, '%f %f %f');
temp = meta.spacedirections;
temp(temp=='(')=[];
temp(temp==')')=[];
temp(temp==',')=' ';
% [a, b, c, d, e, f, g, h, i] = zeros(9);
if contains(temp, 'none') 
    [a, b, c, d, e, f, g, h, i] = strread(temp, 'none %f %f %f %f %f %f %f %f %f');
else
    [a, b, c, d, e, f, g, h, i] = strread(temp, '%f %f %f %f %f %f %f %f %f');
end
sx = a;
sy = e;
sz = i; 
% sx = temp(1); sy = temp(5); sz = temp(9);