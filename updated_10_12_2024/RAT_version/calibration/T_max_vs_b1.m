% maximum test sensitivity b1

b1 = -1:0.001:3

s = 1 ./ (1 + exp(-b1))

plot(b1, s)