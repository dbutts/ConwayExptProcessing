function [xc_pix_abs, yc_pix_abs] = find_sta_center(STA, stim_location)

xleft = stim_location(:,1); % from data.stim_location
ytop = stim_location(:,2);

stim_x = max(xleft) - ScreenSizeX_pix/2;
stim_y = min(ytop)+60 - ScreenSizeY_pix/2;

E = sum(sum(STA.^2,5),4);

% for each STA, shift and find brightest pixels

xc = nan(size(E,1),1);
yc = nan(size(E,1),1);
for i = 1:size(E,1)
    temp = circshift(reshape(E(i,:,:), size(E,2), size(E,3)), [fix(size(E,2))/2, fix(size(E,3))/2]);
    [y,x] = find(temp > prctile(temp(:), 99));

    if std(x) < 5 && std(y) < 5
        xc(i) = mod(mean(x)+fix(size(E,3))/2, size(E,3));
        yc(i) = mod(mean(y)+fix(size(E,2))/2, size(E,2));
    end

end

% corner to which center is closest is the intersection of cloud tiles that
% was placed over handmapped RFs

C = [0 0
    size(E,3) 0
    size(E,3) size(E,2)
    0 size(E,2)
    ];

xc_pix_abs = nan(size(E,1),1);
yc_pix_abs = nan(size(E,1),1);

for i = 1:size(E,1)
    if isfinite(xc(i)) && isfinite(yc(i))
        distFromCorners = sqrt(sum(([xc(i) yc(i)] - C).^2,2));
        closestCorner = C(distFromCorners == min(distFromCorners),:);
        closestCorner = closestCorner(1,:);

        diffFromClosestCorner = [xc(i) yc(i)] - closestCorner;

        xc_pix_abs(i) = stim_x + diffFromClosestCorner(1);
        yc_pix_abs(i) = stim_y + diffFromClosestCorner(2);
    end
end

figure, scatter(xc_pix_abs, yc_pix_abs);

end