%%
clear all
basedir='Z:\StimulusSet\NTlab_cis\QZ\p\';
cd(basedir)

%%
SFvec=[2,4,8,16,32,64,128]; nSF=length(SFvec);
cirvec=[1,2,3]; nCir=length(cirvec);
dirvec=[-1, 1]; nDir=length(dirvec);
movec=[-1, 1]; nMov=length(movec);
%% test reading a file
sf=1; cir=2; dir=1; mov=1;
cur_params=sprintf('sf_%03d_cir_%d_dire_%d_mo_%d', SFvec(sf), cirvec(cir), dirvec(dir), movec(mov));
% %
cur_img=1;
cur_im=sprintf('im%03d.bmp', cur_img);
[image, colormap]=imread([basedir cur_params '\' cur_im]);
% %
test1=colormap(image+1,:);
test2=reshape(test1,64,64,3);
% %
figure;
imagesc(test2)

%%
ii=1;
for sf=1:nSF
    for cir=1:nCir
        for dir=1:nDir
            for mov=1:nMov
                cur_params=sprintf('sf_%03d_cir_%d_dire_%d_mo_%d', SFvec(sf), cirvec(cir), dirvec(dir), movec(mov));
                all_param_inds(sf,cir,dir,mov)=ii;
                for cur_img=1:128
                    cur_im=sprintf('im%03d.bmp', cur_img);
                    [image, colormap]=imread([basedir cur_params '\' cur_im]);
                    temp1=colormap(image+1,:);
                    all_img{ii,:,:,:}=reshape(temp1,64,64,3);
                    all_img_params(ii,:)=[SFvec(sf), cirvec(cir), dirvec(dir), movec(mov), ii];
                    ii=ii+1;
                end
            end
        end
    end
end
% now get uniform cycles
for cir=1:nCir
    for dir=1:nDir
        cur_params=sprintf('uniform_cir_%d_dire_%d', cirvec(cir), dirvec(dir));
        all_param_inds(nSF+1,cir,dir,1)=ii;
        for cur_img=1:128
            cur_im=sprintf('im%03d.bmp', cur_img);
            [image, colormap]=imread([basedir cur_params '\' cur_im]);
            temp1=colormap(image+1,:);
            all_img{ii,:,:,:}=reshape(temp1,64,64,3);
            all_img_params(ii,:)=[0, cirvec(cir), dirvec(dir), 1, ii];
            ii=ii+1;
        end
    end
end
SFvec=[2,4,8,16,32,64,128,0];
%%
save([basedir 'all_rainbows.mat'], 'all_img', 'all_img_params', 'all_param_inds', 'SFvec', 'cirvec', 'dirvec', 'movec');
%%
disp('Done')
%%
%figure; imagesc(all_img{(128*7)+10})