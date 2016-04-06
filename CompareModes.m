close;
[tf ind] = ismember(mkn,T.mkn,'rows');
subplot(2,1,1);
if tf
    %semilogy(T.t,abs(T.coefs{ind}),'b','linewidth',2);
    plot(T.t,angle(T.coefs{ind}),'bx-');
else
    error('mode missing from Teukolsky');
end
hold on;
[tf ind] = ismember(mkn,K.mkn,'rows');
t = linspace(T.t(1),T.t(end),500);
if tf
    %semilogy(t,abs(ppval(K.coefs_pp{ind},t)),'r','linewidth',2);
    %plot(t,angle(ppval(K.coefs_pp{ind},t)),'r','linewidth',2);
    plot(K.t,angle(K.coefs{ind}),'ro-');
else
    error('mode missing from Kludge');
end
ylabel('phase in radians');
xlabel('t/M');
title(['\omega = ' num2str(mkn(1)) '\omega_\phi +' ...
    num2str(mkn(2)) '\omega_\theta +' ...
    num2str(mkn(3)) '\omega_r']);
legend('Teukolsky','Kludge',3);
[tf ind] = ismember(mkn,T.mkn,'rows');
subplot(2,1,2);
if tf
    semilogy(T.t,abs(T.coefs{ind}),'bx-');
else
    error('mode missing from Teukolsky');
end
hold on;
[tf ind] = ismember(mkn,K.mkn,'rows');
t = linspace(T.t(1),T.t(end),500);
if tf
    %semilogy(t,abs(ppval(K.coefs_pp{ind},t)),'r','linewidth',2);
    semilogy(K.t,abs(K.coefs{ind}),'ro-');
else
    error('mode missing from Kludge');
end
ylabel('(D/M) |h_{mkn}|');
xlabel('t/M');
legend('Teukolsky','Kludge',3);