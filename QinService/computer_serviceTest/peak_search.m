function [music,fangwei,fuyang] = peak_search(pmusic,fangwei_low,fuyang_low,fw_precision,fy_precision)
row = size(pmusic,1);
col = size(pmusic,2);
num = 0;
row_index = zeros(1,1);
col_index = zeros(1,1);
music = zeros(1,1);
for ii = 1:row
    for jj = 1:col
        if ii == 1
            up = 0;
            down = pmusic(ii+1,jj);
            up_left = 0;
            up_right = 0;
            if jj == 1
                left = 0;
                right = pmusic(ii,jj+1);
                down_left = 0;
                down_right = pmusic(ii+1,jj+1);
            elseif jj == col
                left = pmusic(ii,jj-1);
                right = 0;
                down_left = pmusic(ii+1,jj-1);
                down_right = 0;
            else
                left = pmusic(ii,jj-1);
                right = pmusic(ii,jj+1);
                down_left = pmusic(ii+1,jj-1);
                down_right = pmusic(ii+1,jj+1);
            end
        elseif ii == row
            up = pmusic(ii-1,jj);
            down = 0;
            down_left = 0;
            down_right = 0;
            if jj == 1
                left = 0;
                right = pmusic(ii,jj+1);
                up_left = 0;
                up_right = pmusic(ii-1,jj+1);
            elseif jj == col
                left = pmusic(ii,jj-1);
                right = 0;
                up_left = pmusic(ii-1,jj-1);
                up_right = 0;
            else
                left = pmusic(ii,jj-1);
                right = pmusic(ii,jj+1);
                up_left = pmusic(ii-1,jj-1);
                up_right = pmusic(ii-1,jj+1);
            end
        elseif jj == 1
            up = pmusic(ii-1,jj);
            down = pmusic(ii+1,jj);
            left = 0;
            right = pmusic(ii,jj+1);
            up_left = 0;
            up_right = pmusic(ii-1,jj+1);
            down_left = 0;
            down_right = pmusic(ii+1,jj+1);
        elseif jj == col
            up = pmusic(ii-1,jj);
            down = pmusic(ii+1,jj);
            left = pmusic(ii,jj-1);
            right = 0;
            up_left = pmusic(ii-1,jj-1);
            up_right = 0;
            down_left = pmusic(ii+1,jj-1);
            down_right = 0;
        else
            up = pmusic(ii-1,jj);
            down = pmusic(ii+1,jj);
            left = pmusic(ii,jj-1);
            right = pmusic(ii,jj+1);
            up_left = pmusic(ii-1,jj-1);
            up_right = pmusic(ii-1,jj+1);
            down_left = pmusic(ii+1,jj-1);
            down_right = pmusic(ii+1,jj+1);
        end
        tar = pmusic(ii,jj);
        if tar > up && tar > down && tar > left && tar > right
            if tar > up_left && tar > up_right && tar > down_left && tar > down_right
                num = num + 1;
                row_index(1,num) = ii;
                col_index(1,num) = jj;
                music(1,num) = tar;
            end
        end
    end
end
fangwei = (row_index - 1)*fw_precision + fangwei_low;
fuyang = (col_index - 1)*fy_precision + fuyang_low;
end
