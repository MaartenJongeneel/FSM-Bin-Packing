function [H_filt,V,dV] = ReconstructAndFilter(H,time,w, fwp, fwo, p)
%UNTITLED2 Summary of this function goes here



    if nargin == 1
        w = 9;
        fwp = 100;
        fwo = 100;
        p = 2;
    end

    % remove outliers
    [~,TP] = rmoutliers(squeeze(H(1:3,4,:))',"movmean",w,"SamplePoints",time);
    for jj=1:size(H,3)   
        RotVec(jj,:) = vee( logm( H(1:3,1:3,jj)/H(1:3,1:3,1) ) )';
    end
    [~,TPhi] = rmoutliers(RotVec,"movmean",w,"SamplePoints",time);

    % interpolate missing data
    o(1,:) = interp1(time(~TP),squeeze(H(1,4,~TP)),time,'pchip');
    o(2,:) = interp1(time(~TP),squeeze(H(2,4,~TP)),time,'pchip');
    o(3,:) = interp1(time(~TP),squeeze(H(3,4,~TP)),time,'pchip');
    phi(1,:) = interp1(time(~TPhi),RotVec(~TPhi,1),time,'pchip');
    phi(2,:) = interp1(time(~TPhi),RotVec(~TPhi,2),time,'pchip');
    phi(3,:) = interp1(time(~TPhi),RotVec(~TPhi,3),time,'pchip');

    % build back the pose time series
    H_ = zeros(size(H));
    for jj=1:size(H,3)   
        H_(:,:,jj) = [expm(hat(phi(:,jj)))*H(1:3,1:3,1) , o(:,jj); zeros(1,3) 1];
    end

    % filter the data
    [H_filt,V,dV] = FilteringSE3(H_,time, fwp, fwo, p);
                
end

