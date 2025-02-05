import matplotlib.pyplot as plt

h_patch = [0.1, 0.075, 0.05, 0.025, 0.015, 0.01]
error_patch = [7.97E-5, 6.95E-7, 1.95E-7, 3.27E-8, 9.64E-9, 3.26E-9]
plt.plot(h_patch, error_patch, label='$||u^p_h - g_{0.01}||_{L^2(\Omega \cup \mathcal{I})}$', ls='-', color='r',
         marker='o')

h_dir = [0.1, 0.075, 0.05, 0.025, 0.015]
error_dir = [2.97E-5, 1.70E-5, 9.36E-6, 2.32E-6, 3.51E-7]

plt.plot(h_dir, error_dir, label='$||u^D_h - u^D_{0.01}||_{L^2(\Omega \cup \mathcal{I})}$', ls='-', color='b',
         marker='o')

h_dir_2 = [0.1, 0.075, 0.05, 0.025, 0.015]
error_dir_2 = [1.31E-4, 1.01E-4, 6.62E-5, 2.19E-5, 4.38E-6]

plt.plot(h_dir_2, error_dir_2, label='$||u^{D2}_h - u^{D2}_{0.01}||_{L^2(\Omega \cup \mathcal{I})}$', ls='-',
         color='c', marker='o')

h_neu = [0.1, 0.075, 0.05, 0.025, 0.015]
error_neu = [3.76E-5, 1.82E-5, 1.04E-5, 4.23E-6, 1.46E-6]

plt.plot(h_neu, error_neu, label='$||u^N_h - u^N_{0.01}||_{L^2(\Omega \cup \mathcal{I})}$', ls='-', color='g',
         marker='o')

h_ref = [0.1, 0.075, 0.05, 0.025, 0.015]
error_ref = [entry**2 for entry in h_ref]

plt.plot(h_ref, error_ref, label='$h^2$', ls='-', color='gray', marker='o')

plt.xlabel('h')
plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.2), ncol=2, framealpha=1)
plt.tight_layout()
plt.savefig('h_error.pdf')
plt.show()
