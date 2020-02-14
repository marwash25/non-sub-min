function [W,tab ] = read_genrmf(file);
% read output from genrmf graph data files
x = textread(file,'%s');
a = str2double(x(13));
b = str2double(x(15));
c1 = str2double(x(17));
c2 = str2double(x(19));
seed = str2double(x(21));
p = str2double(x(24));
s = str2double(x(27));
t = str2double(x(30));
tab = [];
for i=33:4:length(x);
    tab = [ tab; [ str2double(x(i)), str2double(x(i+1)), str2double(x(i+2))] ];
end
W =sparse(tab(:,1),tab(:,2),tab(:,3),p,p);


