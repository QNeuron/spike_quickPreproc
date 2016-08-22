function [pos]=subplot_pos(n,p,edgel,edger,edgeh,edgeb,space_h,space_v)
%provide positions for n,p compact subplot 

if nargin<3,
    edgel=0.07;
elseif isempty(edgel),
    edgel=0.07;
end;
if nargin<4,
    edger=0.05;
elseif isempty(edger),
    edger=0.05;
end;
if nargin<5,
    edgeh=0.07;
elseif isempty(edgeh),
    edgeh=0.07;
end;
if nargin<6,
    edgeb=0.07;
elseif isempty(edgeb),
    edgeb=0.07;
end;
if nargin<7,
    space_h=0.07;
elseif isempty(space_h),
    space_h=0.07;
end;
if nargin<8,
    space_v=0.07;
elseif isempty(space_v),
    space_v=0.07;
end;

if edgel>0.5|edgeb>0.5|edger>0.5|edgeh>0.5|edgel<0|edgeb<0|edger<0|edgeh<0|space_h<0|space_h>0.5|space_v<0|space_v>0.5,
    error('values for edges and spaces <0 or >0.5. Are you crazy ?');
end;

width=(1-edger-edgel-(p-1)*space_h)/p;
height=(1-edgeh-edgeb-(n-1)*space_v)/n;
if width<0|height<0,
    error('spaces betwen subplot too large');
end;

ik=0;
for i=1:n,
        for j=1:p,
        ik=ik+1;
        pos{ik}=[edgel+(j-1)*(width+space_h) 1-edgeh+space_v-i*(height+space_v) width height];
    end;
end;

% width=(1-edger-edgel)/p;
% height=(1-edgeh-edgeb)/n;
% 
% ik=0;
% for i=1:n,
%         for j=1:p,
%         ik=ik+1;
%         pos{ik}=[edgel+(j-1)*width 1-edgeh-i*height width height];
%     end;
% end;