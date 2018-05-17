function show_batch(path)
%path=['/Users/yvesrobert/Desktop/Stage/Control_rods/try/twor3/A_two_rows3'];
file=fopen(path,'r');
lines= textscan(file,'%s','Delimiter','\n');
fclose(file);
lines=strtrim(lines{1});

i=1;
 while isempty(strfind(lines{i},'lat 1 '))
     i=i+1;
 end
 while isempty(strfind(lines{i},'0166'))
     i=i+1;
 end
 m=i;

while isempty(strfind(lines{m},'%%'))
    for j=1:8:length(lines{m})
            batch=lines{m}(j:j+3);
            lattice((m-i)/2+1,ceil(j/8))=str2double(batch(end-1:end));            
    end
    m=m+2;
end
while sum(~isnan(lattice(end,:)))==0
    lattice(end,:)=[];
end
while sum(~isnan(latticepath=['/Users/yvesrobert/Desktop/Stage/Control_rods/try/twor3/A_two_rows3'];
file=fopen(path,'r');
lines= textscan(file,'%s','Delimiter','\n');
fclose(file);
lines=strtrim(lines{1});

i=1;
 while isempty(strfind(lines{i},'lat 1 '))
     i=i+1;
 end
 while isempty(strfind(lines{i},'0166'))
     i=i+1;
 end
 m=i;

while isempty(strfind(lines{m},'%%'))
    for j=1:8:length(lines{m})
            batch=lines{m}(j:j+3);
            lattice((m-i)/2+1,ceil(j/8))=str2double(batch(end-1:end));            
    end
    m=m+2;
end
while sum(~isnan(lattice(end,:)))==0
    lattice(end,:)=[];
end
while sum(~isnan(lattice(:,end)))==0
    lattice(:,end)=[];
end
lattice(lattice>50)=0;
plot_core('',ceil(length(lattice)/2),lattice);
(:,end)))==0
    lattice(:,end)=[];
end
lattice(lattice>50)=0;
plot_core('',ceil(length(lattice)/2),lattice);
