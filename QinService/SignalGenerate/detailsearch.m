function [doa_endfind] = detailsearch(doa_detailsearch2,step,Rdoa_Un,location,fre)
%DETAILSEARCH 此处显示有关此函数的摘要
%   此处显示详细说明
count=0;
doa_peak_2 = zeros(1,size(doa_detailsearch2,2));                    %% End find peak
doa_endazmiuth_2 =zeros(1,size(doa_detailsearch2,2));
doa_pitch =zeros(1,size(doa_detailsearch2,2));
lambda = fre\3e8;                                                   %% 波长
am=zeros(8,1);
A_c=3e8;
for i=1:size(doa_detailsearch2,2)
    doa_peak_temp = zeros(1,1600);
    doa_detail2 = cell(1,1600); 
   for j=doa_detailsearch2{i}(1):step:doa_detailsearch2{i}(2)
       for k=doa_detailsearch2{i}(3):step:doa_detailsearch2{i}(4)
           doa_detail2{count+1} =[j,k];
           count2=0;
           for e = 1:2:16
               tao_dir_temple = (location(e,1)*cosd(j)*sind(k)+location(e,2)*sind(j)*sind(k)+location(e,3)*cosd(k))/A_c;
               am(count2+1)=exp(+1i*2*pi*A_c*tao_dir_temple/lambda);
               count2 = count2+1;
           end
           doa_peak_temp(count+1)=1/(am'* (Rdoa_Un* Rdoa_Un')*am);
           count =count+1;
       end
   end
   [doa_peak_value2,peak_index2]=max(doa_peak_temp);
    doa_peak_2(i) = doa_peak_value2;
    doa_endazmiuth_2(i)=doa_detail2{peak_index2}(1);         %% 存取方向角
    doa_pitch(i)=doa_detail2{peak_index2}(2);
    count =0;
end
[~,index]=sort(doa_peak_2,'descend');
doa_endazmiuth_2=doa_endazmiuth_2(1,index);
doa_pitch=doa_pitch(1,index);
doa_endfind = cell(1,size(doa_peak_2,2));

for i=1:size(doa_peak_2,2)
    doa_endfind{i}=[doa_endazmiuth_2(index(i)),doa_pitch(index(i))];
end

end

