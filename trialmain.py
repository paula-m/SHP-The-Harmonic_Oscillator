from metropolis_class import metropolis
import numpy as np
import matplotlib.pyplot as plt
from thermalisation import thermalisation
from error_analysis import error_calculations
#from correlation_function_class import correlation_function
from autocorrelation import autocorrelation
from error_calc import error_calculation
from corr_fun import correlation_function
from puwr import tauint, correlated_data

def main():
    m = 0.1
    lam = 0
    kt = 0
    N = 1200
    h = 3.5
    them = 1000
    sweeps = 1000

    initialise_thermalisation = thermalisation(N, h, m, m, them, lam, kt)
    initial_path = initialise_thermalisation.thermalised_start()

    initialise_metropolis_class = metropolis(initial_path, h, m, m, N, lam, kt)
    output_arrays_x = initialise_metropolis_class.running_the_sweep(sweeps)
    output_arrays_x = np.array(output_arrays_x)

    sweep_len = len(output_arrays_x)
    #print("this is output array", output_arrays_x)
    initialise_error_calculations = error_calculation()

    squared_array_full = initialise_error_calculations.calculate_power(output_arrays_x, 2)

    mean_of_each_x_element = initialise_error_calculations.calculate_the_mean_element(output_arrays_x)
    mean_of_markov_chain_x = initialise_error_calculations.calculate_mean_markov(output_arrays_x)

    print("this is mean of each x", mean_of_each_x_element)
    #print(mean_of_each_x_element.shape)

    mean_of_x2_arrays = np.mean(squared_array_full, axis = 1)
    mean_of_markov_chain_x2 = initialise_error_calculations.calculate_mean_markov(squared_array_full)

    R = 1+ (m**2)*1/2 - m*np.sqrt(1+(m**2)*1/4)
    x2_mean_theoretical = (1/(2*m*m*np.sqrt(1+(m**2)*1/4)))*((1+R**N)/(1-R**N))

    proper_m = np.mean(output_arrays_x, axis = 0)
    markov_mean = (1/(sweep_len*N))*np.sum(proper_m*N)

    print("this is markov mean", markov_mean)

    naive_error_x_array = initialise_error_calculations.calculate_naive_error(output_arrays_x, sweep_len)
    #print("this is the naive error of x", naive_error_x_array)

    markov_arrays = np.mean(output_arrays_x, axis = 1)

    variance_of_each_x = np.var(output_arrays_x, axis=1, ddof=1)
    #print("var shape", variance_of_each_x)
    #naive_error_in_x = np.sqrt(variance_of_each_x/sweep_len)
    print("variance total", np.mean(variance_of_each_x))
    naive_err_try = np.sqrt(np.mean(variance_of_each_x)/N)
    #naive_error_in_MCMC_x_chain = np.mean(naive_error_x_array)
    #print("naive 1", naive_err_try, "mcmc", naive_error_in_MCMC_x_chain)
    print("mean of MCMC" + str(f'{markov_mean: .5}'), '+-', naive_err_try)

    #initialise_autocorrelation = autocorrelation()
    #autocorr = initialise_autocorrelation.calculate_interal_autocor(output_arrays_x[:,0], sweep_len, sweep_len)
    #print("autocorr", autocorr)
    #autocorr_x2 = initialise_autocorrelation.calculate_interal_autocor(mean_of_x2_arrays, sweep_len, 40)

    #try_output = [[output_arrays_x[:,0]]]
    #print(try_output)
    mean_MC, delta, tint, d_tint = tauint([[mean_of_each_x_element]], 0)
    mean_MCx2, varx2, tintx2, dintx2 = tauint([[mean_of_x2_arrays]], 0)
    print("mean, delta, tint, d_tint", mean_MCx2)

    jackknife_error_each_x = []
    jackknife_variance = []
    for i in range(N):
        jaccknife_err, jackknife_var = initialise_error_calculations.jacknife(output_arrays_x[:,i], sweep_len, 78)
        jackknife_error_each_x.append(jaccknife_err)
        jackknife_variance.append(jackknife_var)
    #print("jackknife", jackknife_error_each_x)



    #print("this is the jacknife error for each of the x", jackknife_error_each_x)
    total_error, jackvar = initialise_error_calculations.jacknife(mean_of_x2_arrays, sweep_len, int(tintx2))
    err_jack = np.sqrt(np.mean(jackknife_variance)/N)
    print("err jcakkkf", err_jack)
    print("this is total error", total_error)

    jaccknife_error_written, jaccknife_variance_1 = initialise_error_calculations.jacknife(markov_arrays, sweep_len, int(tint))
    print("this is the original", jaccknife_error_written)
    
    naive_error_x2 = initialise_error_calculations.calculate_naive_error(squared_array_full, sweep_len)

    jackknife_error_each_x2 = []
    jackknife_variance_x2 = []
    for j in range(sweep_len):
        jaccknife_err, jackknife_var = initialise_error_calculations.jacknife(squared_array_full[j], N, 5)
        jackknife_error_each_x2.append(jaccknife_err)
        jackknife_variance_x2.append(jackknife_var)

    
    total_error_x2 = np.sqrt(np.sum(jackknife_variance_x2)/sweep_len)
    #jackknife_error_total_x2, jaccknife_var_2 = initialise_error_calculations.jacknife(mean_of_x2_arrays, sweep_len, int(autocorr_x2))
    print("this is the jacknife error x2", total_error_x2)
    

    #corre_fun = correlation_function(output_arrays_x, N)
    #display = corre_fun.plotting_correlation_function()

    """Plotting of the arrays"""
    
    array_of_theoretical_vals = np.ones(sweep_len)*x2_mean_theoretical
    array_of_mean_markov_value = np.ones(sweep_len)*mean_of_markov_chain_x2
    y = np.arange(sweep_len)
    plt.scatter(y, mean_of_x2_arrays, label=r'$\langle x^2 \rangle$ calculated for each path')
    plt.plot(y, array_of_theoretical_vals, label = r'theoretical value of $\langle x^2 \rangle$' + str(f'{x2_mean_theoretical: .5f}'))
    plt.plot(y, array_of_mean_markov_value, label = r'calculated $\langle x^2 \rangle$' + str(f'{mean_of_markov_chain_x2: .5f}'))
    #plt.errorbar(y, mean_of_x2_arrays, yerr=jackknife_error_each_x2, ecolor='purple')
    #plt.errorb#ar(y, mean_of_x2_arrays, yerr=naive_error_x2, ecolor='pink')
    plt.legend()
    plt.xlabel('Metropolis sweeps')
    plt.ylabel(r"$ \langle x^2 \rangle$")
    plt.title(r'Graph showing the $\langle x^2 \rangle$ vs the number of Metropolis sweeps')
    plt.show()
    

    y = np.arange(N)
    plt.plot(y, mean_of_each_x_element, label='The path of x')
    plt.errorbar(y, mean_of_each_x_element, yerr=jackknife_error_each_x, ecolor='purple', fmt='none', label='Jackknife errors')
    plt.errorbar(y, mean_of_each_x_element,yerr=naive_error_x_array, ecolor='pink' , fmt='none', label='Naive errors')
    plt.xlabel(r'$\tau$')
    plt.ylabel('x value')
    plt.title('The path of x along one element for the full Markov Chain')
    plt.legend()
    plt.show()

    #error_x = np.reshape(jackknife_error_each_x, (-1))
    x_vals_n = np.reshape(output_arrays_x, (-1))
    num_bins = np.arange(-30, 30, 0.5)
    n, bins, patches = plt.hist(x_vals_n, num_bins, facecolor='purple', alpha=0.7, density=True)
    theor = (((m*m/np.pi)**0.25)*np.exp(-m*m*(num_bins**2)/2))**2
    #thr2 = 0.59*np.exp(-1.1*num_bins**2)
    #plt.errorbar(num_bins, x_vals_n, error_x)
    #plt.plot(num_bins, thr2, label='Finite lattice spacing theory curve', color='blue')
    plt.plot(num_bins, theor, label='Infinite lattice spacing theory curve')
    plt.xlabel("x")
    plt.ylabel(r'$|{\psi_0(x)}^2|$')
    plt.title(r'Graph showing ground function $|{\psi}^2|$ vs each value of x')
    plt.legend()
    plt.show()
    print(np.sum(np.diff(bins)*n))

    result = open('resultt.txt', 'a+')
    arr1 = [[m, N, kt, h, them, lam, sweeps, markov_mean, naive_err_try, mean_MCx2, varx2, tintx2, dintx2, mean_MC, delta, tint, d_tint, sweep_len]]
    print(arr1, file= result)
    result.close()
main()