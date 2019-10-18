%% Demo of PCMI problems

%% Demo of clustering
clear all;

pExampleNames   = {
    'I 1','II 1(a)','II 1(b)'};

fprintf('\n\n Select example to run:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

if (~exist('pExampleIdx') || isempty(pExampleIdx) || pExampleIdx==0)
    pExampleIdx          = input('Pick an example to run:           ');
end


scrsz = get(groot,'ScreenSize');


switch pExampleNames{pExampleIdx}
    case 'I 1'
        gcf1=figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
        gcf2=figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
        for logD = 1:10
            figure(gcf1);
            D        = 2^logD; ID = eye(D);
            x        = (mvnrnd(zeros(D,1),ID,10000))';                                 % guassian
            normxsq  = sum(x.^2,1);                                                    % norm
            subplot(2,5,logD);
            hist(normxsq)                                                      % histogram
            xlabel('||\cdot||^2','interpreter','TeX');title(sprintf('D=%d',D));
            figure(gcf2);
            subplot(2,5,logD);
            histc_res=histc(normxsq/mean(normxsq),linspace(0,2,100));
            bar(linspace(0,2,100),histc_res,'histc');
            xlabel('||\cdot||^2/E[||\cdot||^2]');title(sprintf('D=%d',D));
        end
        figure(gcf1);
    case 'II 1(a)'
        %% disk
        N        = 64;
        r        = 10;     tstep = 2;
        [xx,yy]  = meshgrid(1:N,1:N);
        zz       = sqrt((xx-(r+1)).^2+(yy-(r+1)).^2)<=r; % disk
        trans    = 0:tstep:(N-2*r);
        [tt1,tt2]= meshgrid(trans,trans);
        tt1      = tt1(:);
        tt2      = tt2(:);
        trans    = [tt1 tt2];
        tLen     = size(trans,1);
        X        = zeros(N^2,tLen);
        % generate translation
        for i = 1:tLen
            img    = imtranslate(zz,[trans(i,1), trans(i,2)]);
            X(:,i) = img(:);
        end
        xmean       = sum(X,2)/tLen;
        [U , S, ~]  = svd(bsxfun(@minus,X,xmean));
        S           = diag(S);
        % display image
        figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]); subplot(1,3,1); imagesc(zz); colormap(gray); title('image')
        % Display data matrix and singular values
        subplot(1,3,2); imagesc(X); colormap(gray); title('data matrix')
        subplot(1,3,3); plot(S); title('singular value')
        % display eigenfunctions
        figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4])
        for i = 1:16
            subplot(4,4,i);  imagesc( reshape(U(:,i),[N,N])); colormap(gray); title(num2str(i))
        end
        %% square
        zz       = (xx-(r+1)<=r).*(yy-(r+1)<=r); % square
        trans    = 0:tstep:(N-2*r);
        [tt1,tt2]= meshgrid(trans,trans);
        tt1      = tt1(:);
        tt2      = tt2(:);
        trans    = [tt1 tt2];
        tLen     = size(trans,1);
        X        = zeros(N^2,tLen);
        % generate translation
        for i = 1:tLen
            img    = imtranslate(zz,[trans(i,1), trans(i,2)]);
            X(:,i) = img(:);
        end
        xmean       = sum(X,2)/tLen;
        [U , S, ~]  = svd(bsxfun(@minus,X,xmean));
        S           = diag(S);
        % display image
        figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]); subplot(1,3,1); imagesc(zz); colormap(gray); title('image')
        % Display data matrix and singular values
        subplot(1,3,2); imagesc(X); colormap(gray); title('data matrix')
        subplot(1,3,3); plot(S); title('singular value')
        % display eigenfunctions
        figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4])
        for i = 1:16
            subplot(4,4,i);  imagesc( reshape(U(:,i),[N,N])); colormap(gray); title(num2str(i))
        end
    case 'II 1(b)'
        D      = 225;  ID = rand(D,D); ID=ID*ID';  [V , ~] = eig(ID);
        Rank   = (8:floor(sqrt(D))).^2;  LenR = length(Rank);   % rank
        k      = 100;
        Sample = round(10.^(3:0.2:5));    LenS   = length(Sample);
        Error  = zeros(LenR,LenS);
        for r = 1:LenR
            k = Rank(r);  ID = zeros(D,D);   for i=1:k; ID(i,i) = 1; end
            ID     = V'*ID*V;
            for i = 1:LenS
                n      = Sample(i);
                X      = (mvnrnd(zeros(D,1),ID,n))';
                xmean  = sum(X,2)/n;
                X      = bsxfun(@minus,X,xmean);
                eSigma    = X*X'/n;
                Error(r,i)  = norm(eSigma-ID);
            end
        end
        allcolor = linspace(0,1,LenR);  figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]); leng =cell(LenR,1);
        for r=1:LenR
            plot(log10(Sample),log10(Error(r,:)),'color',[allcolor(r) 0 0]); hold on;
            rate    = polyfit(log10(Sample),log10(Error(r,:)),1);
            leng{r} = ['rank = ' num2str(Rank(r)) ' slope = '  num2str(rate(1))];
        end
        legend(leng)
        xlabel('log(n)'); ylabel('log(error)'); title('Error of sampled covariance matrix versus n')
end

