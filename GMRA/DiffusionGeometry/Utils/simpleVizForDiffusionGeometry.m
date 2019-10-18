%% Simple visualization for Diffusion Geometry
% Original data and diffusion embeddings
scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
subplot(1,3,1);plot3(X(1,:),X(2,:),X(3,:),'.');
xlabel('x_1','Interpreter','tex');ylabel('x_2','Interpreter','tex');zlabel('x_3','Interpreter','tex');
title('Original data, first 3 coords');
subplot(1,3,2);plot3(Data.G.EigenVecs(:,2),Data.G.EigenVecs(:,3),Data.G.EigenVecs(:,4),'.');
xlabel('\phi_2','Interpreter','tex');ylabel('\phi_3','Interpreter','tex');zlabel('\phi_4','Interpreter','tex');
title('Diffusion embedding (2,3,4)');
subplot(1,3,3);plot3(Data.G.EigenVecs(:,5),Data.G.EigenVecs(:,6),Data.G.EigenVecs(:,7),'.');
xlabel('\phi_5','Interpreter','tex');ylabel('\phi_6','Interpreter','tex');zlabel('\phi_7','Interpreter','tex');
title('Diffusion embedding (5,6,7)');


% Original data colored by values of eigenfunctions
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
for k = 1:15
    subplot(3,5,k)
    scatter3(X(1,:),X(2,:),X(3,:),20,Data.G.EigenVecs(:,k),'filled');
    title(sprintf('Eigenfunction %d on the data',k));
end