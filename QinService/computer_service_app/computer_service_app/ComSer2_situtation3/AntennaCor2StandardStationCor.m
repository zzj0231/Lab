function [sta_fw,sta_fy,alph,beta] = AntennaCor2StandardStationCor(H,fw,fy,ant_fw,ant_fy)

% 根据当前天线阵所在的xyz磁场分量、磁偏角、磁倾角信息，将天线阵doa估计得到的角度
% 信息转换为标准站心坐标系下的方位角和俯仰角。后面的定位都是在标准站心坐标系下进行的

% H：三磁场分量，(Hx,Hy,Hz)1*3的分量
% fw：磁偏角，角度，与y轴夹角（地磁N级）
% fy：磁倾角，角度，与xoy面的夹角
% ant_fw：天线阵坐标系中的方位角，角度
% ant_fy：天线阵坐标系中的俯仰角，角度
% 
% sta_fw：站心坐标系的方位角，角度
% sta_fy：站心坐标系的俯仰角，角度

H_value = sqrt(H*H');
x1 = [H_value*cosd(fw)*cosd(fy);H_value*sind(fw)*cosd(fy);H_value*sind(fy)];
x2 = H';
A = (x1*x2')*(x2*x2')^(-1);
beta = acosd(H_value*sind(fy)/sqrt(H(1:2)*H(1:2)')) + atand(H(1)/H(3));
alph = acosd(H_value*cosd(fy)*sind(fw)/sqrt((H(1,1)*cosd(beta) - H(1,3)*sind(beta))^2 + H(1,2)^2)) - atand(H(1,3)/(H(1,1)*cosd(beta) - H(1,3)*sind(beta)));
ant = [cosd(ant_fy)*cosd(ant_fw);cosd(ant_fy)*sind(ant_fw);sind(ant_fy)];
sta = A*ant;
sta_fw = atand(sta(2,1)/sta(1,1));
sta_fy = asind(sta(3,1));


end

