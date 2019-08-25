# VALSE-EP
Implement the line spectral estimation from quantized measurements.
This code (VALSE-EP) is built based on the VALSE algorithm, and is written by Jiang Zhu, Qi Zhang and Xiangming Meng. 
The VALSE-EP is designed to perform line spectral estimation from quantized samples. 
If you have any problems about the code, please feel free to contact jiangzhu16@zju.edu.cn


Please directly run NMSEvsIter.m to see the performance of the VALSE-EP. 
In addition, the VALSE code (written by Mihai-Alin Badiu) is also provided for performance comparison.


Main function：
out = VALSE_EP( y_q, m, ha, x, Iter_max, B, yy_min, alpha, method_EP )

Input parameters：

y_q：For quantization, it belongs to $0,1,2,\cdots,2^B-1$. For unquantized setting, it is the unquantized measurmenets.

m: The index correspond to incomplete measurements 

ha: Set ha=2, which corresponds to Heuristic 2

x: True spectral signal

Iter_max：The maximum number of iterations

B: The bit-depth. For unquantized setting, set B=inf.

yy_min: The left endpoints of the quantizer

alpha: The stepsize of the quantizer

method_EP：'diag_EP' or scalar_EP’. Please set 'diag_EP'.

Output parameters:

out=struct('freqs',th,'amps',w(s),'x_estimate',xr,'noise_var',nu,'iterations',t,'MSE',mse,'K',Kt);

freqs: The point estimates of the frequency

amps: The complex weight amplitude

xr：Reconstructed line spectral

noise_var: noise variance estimate

t: The iterations exited 

MSE：Normalized MSE

K: Tracking the number of spectral during iteration


