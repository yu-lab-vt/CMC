Input of "CMC_V0.exe": margin totals of a 2D/3D tensor
Output of "CMC_V0.exe": exp(r_i), exp(w_j), and exp(u_k).

We offer a demo located in the "Demo" folder. Please be aware that currently, the demo can only be run on the Windows operating system. Within the "margins" subfolder of the demo, you will find files named "m0.bin," "m1.bin," and "m2.bin," which store the margin totals of each dimension.

To run the demo, simply double-click on "CMC_V0.exe." Then exp(r_i), exp(w_j), and exp(u_k) will be written to the files "exp_r.txt," "exp_w.txt," and "exp_u.txt," respectively. Please note that the files store the values of exp(r_i), exp(w_j), and exp(u_k), not their original values (r_i, w_j, and u_k).

If you wish to try other tensors, we provide a MATLAB code, "Demo/WriteMargins2bin.m," which generates a random tensor and writes its margin into .bin files.

