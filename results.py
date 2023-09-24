dx = [2,5,10,50,100]
expected_d = [32**(1-1/d) for d in dx]
expected_d2 = [2*32**(1-1/d) for d in dx]
min_rate_d= [9.52,24.22,31.12,31.69,31.60]
min_rate_max_d = [15,30,32,32,32]
min_d = [15,30,32,32,32]
min_max_d = [9,30,32,32,32]
rate_d = [8.28,24.03,31.11,31.69,31.60]
rate_max_d = [13,30,32,32,32]

tx = [16,32,64,128]
expected_t = [t**(1/2) for t in tx]
expected_t2 = [2*t**(1/2) for t in tx]
min_rate_t = [5.76,9.51,12.89,18.28]
min_rate_max_t = [8,15,19,26]
min_t = [5.95,8.25,11.71,15.55]
min_max_t=[9,12,17,22]
rate_t = [5.85,8.28,11.76,16.19]
rate_max_t = [8,13,17,22]

nx = [2048,4096,8192,16384]
expected_n = [32**(1/2) for _ in nx]
expected_n2 = [2*32**(1/2) for _ in nx]
min_rate_n = [8.97,9.15,9.52,9.55]
min_rate_max_n = [13,14,15,14]
min_n = [7.60,8.03,8.25,8.59]
min_max_n = [12,12,12,13]
rate_n = [7.90,7.91,8.28,8.6]
rate_max_n = [12,12,13,13]

import matplotlib.pyplot as plt

plt.plot(nx,expected_n,color='black',label='32^(1/2)',linestyle='--')
plt.plot(nx,expected_n2,color='black',label='2*32^(1/2)',linestyle='-.')
plt.plot(nx,min_rate_max_n,color='blue',label='arbitrary min kappa',marker='x')
plt.plot(nx,min_rate_n,color='blue',label='arbitrary min avg',marker='x',linestyle='dotted')
plt.plot(nx,min_max_n,color='green',label='min kappa',marker=5)
plt.plot(nx,min_n,color='green',label='min avg',marker=5,linestyle='dotted')
plt.plot(nx,rate_max_n,color='red',label='uniform kappa',marker=10)
plt.plot(nx,rate_n,color='red',label='uniform avg',marker=10,linestyle='dotted')
plt.legend()
plt.show()
