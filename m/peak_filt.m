function [MaxTab2] = peak_filt(x)

audio_colors;

[MaxTab, MinTab] = peakdet(x, 1e-4);

% figure;
%     plot(x,'color',green2)
%     hold on
%     % plot(locs,pks,'.','MarkerSize',20,'color',orange1);
%     plot(MaxTab(:,1),MaxTab(:,2),'.','MarkerSize',20,'color',blue1);
%     plot(MinTab(:,1),MinTab(:,2),'.','MarkerSize',15,'color',red2);
%     title('\fontsize{16}Autocorrelacion')
%     legend('Autocorr','Maximos locales','Minimos locales')

dists = zeros(size(MaxTab,1),1);
for i=1:size(MaxTab,1)
    min_l = find(MinTab(:,1)<MaxTab(i,1),1,'last');
    min_r = find(MinTab(:,1)>MaxTab(i,1),1,'first');
    if ~isempty(min_l) & ~isempty(min_r)
        dist = (2*MaxTab(i,2)-MinTab(min_l,2)-MinTab(min_r,2))/2;
    elseif isempty(min_l) & ~isempty(min_r)
        dist = MaxTab(i,2)-MinTab(min_r,2);
    elseif ~isempty(min_l) & isempty(min_r)
        dist = MaxTab(i,2)-MinTab(min_l,2);
    else
        disp('Error al eliminar los maximos espureos de la autocorrelacion');
        break
    end
    dists(i) = dist;
end

th = 0.10;
thr = min(dists)+(max(dists)-min(dists))*th;

MaxTab2 = MaxTab(dists>thr,:);