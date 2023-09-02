figure(2)
subplot(2,2,1)
histogram(tril(rawncallon{1}(:)),distax)
hold on
histogram(tril(rawncalloff{1}(:)),distax)
xlabel('noise corrs')
ylabel('no. of pairs')
legend('Stim ON','Stim OFF')
title('Makkay 1')
grid on;
box off


figure(2)
subplot(2,2,2)
histogram(tril(rawncallon{2}(:)),distax)
hold on
histogram(tril(rawncalloff{2}(:)),distax)
xlabel('noise corrs')
ylabel('no. of pairs')
legend('Stim ON','Stim OFF')
title('Makkay 2')
grid on;
box off
subplot(2,2,3)
histogram(tril(rawncallon{3}(:)),distax)
hold on
histogram(tril(rawncalloff{3}(:)),distax)
xlabel('noise corrs')
ylabel('no. of pairs')
legend('Stim ON','Stim OFF')
title('Dino 1')
grid on;
box off


subplot(2,2,4)
histogram(tril(rawncallon{4}(:)),distax)
hold on
histogram(tril(rawncalloff{4}(:)),distax)
xlabel('noise corrs')
ylabel('no. of pairs')
legend('Stim ON','Stim OFF')
title('Dino 2')
grid on;
box off
suptitle('Distribution of Raw noise corrs')