function spec_new = my_shift(spec_old, frequency_shift)
% lr: shifted amount (+ left) (- right).

k = abs(frequency_shift);

if frequency_shift <= 0
    spec_new = spec_old([end - k + 1 : end, 1 : end - k]); % shift right/down k elements
else
    spec_new = spec_old([k + 1 : end, 1 : k]); % shift left/up k elements
end
end