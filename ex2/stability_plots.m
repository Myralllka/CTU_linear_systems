clear all; close all; clc

A = [1/2, sqrt(3)/2; 
      -sqrt(3)/2, 1/2]

% A = [1/sqrt(2), 1/sqrt(2); 
    % -1/sqrt(2), 1/sqrt(2)]

x = [0.01; 0.04]
n_iter = 100;
X1s = zeros(1,n_iter); % Memory preallocation
X2s = zeros(1,n_iter); % Memory preallocation
X1s(1) = x(1)
X2s(1) = x(2)

for k=2:n_iter     
  x = A*x
  X1s(k) = x(1)
  X2s(k) = x(2)
end

figure(1)

hold on;              
plot(X1s, X2s)
hold off;             
xlabel('X1');         
ylabel('X2');
title('trajectory');
grid on;
daspect([1 1 1])