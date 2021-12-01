from math import log
import numpy as np
import matplotlib.pyplot as plt
from error_analysis import error_calculations
from autocorrelation import autocorrelation
from error_calc import error_calculation

from puwr import tauint, correlated_data

class correlation_function():
    def __init__(self, array, nt):
        self.array = array
        self.N_t = nt

    def correlation_function_harmonic(self, tau):
        
        g_vals = []
        for i in range(len(self.array)):
            g_fun = 0
            g_tau = 0
            g_tau_c = 0
            for j in range(len(self.array[i])):
                #for k in range(len(self.array)):
                    #if (k-j)%self.N_t == tau:
                g_val =  self.array[i][j]*self.array[i][(j+tau)%self.N_t]
                g_tau += self.array[i][j]
                g_tau_c += self.array[i][(j+tau)%self.N_t]
                g_fun += g_val
            g_vals.append(g_fun*(1/self.N_t) -(g_tau/self.N_t)*(g_tau_c/self.N_t))
        return g_vals


    def plotting_correlation_function(self):
        
        corr_fun_vals = []
        corr_fun_dist = np.arange(0, 31)
        for i in range(0, 31):
            correlation_fun = correlation_function.correlation_function_harmonic(self, i)
            corr_fun_vals.append(correlation_fun)
            #error_cal.append(np.array(values))
        #print("this is correlation function", corr_fun_vals)
        #print("this is length of correlation function", len(error_cal[0]))

        
        autocorr = autocorrelation()

        err = []
        means = []
        for l in range(len(corr_fun_vals)):
            mean, var, autoc, erautoc = tauint([[corr_fun_vals[l]]], 0)
            print("mean, var, autoc, erauto", mean, var, autoc, erautoc)
            err.append(np.sqrt(np.var(corr_fun_vals[l], ddof= 1)*autoc/len(corr_fun_vals[l])))
            #err.append(np.sqrt(np.var(corr_fun_vals[l], ddof= 1)*autoc/40))
            #err.append(np.sqrt(var*autoc/len(corr_fun_vals[l])))
            means.append(mean)
        #err_val = np.ones(len(corr_fun_vals))
        #print("autocor", autocorr.calculate_interal_autocor(corr_fun_vals[-1], len(corr_fun_vals[-1]), len(corr_fun_vals[-1])))
        #print("np.mean", np.mean(corr_fun_vals, axis =1))
        #print("error", err)
        
        #variance_vals = np.ones(len(corr_fun_vals))
        #err = np.ones(len(corr_fun_vals))
        
        #auto_corr = []
        #for j in range(len(corr_fun_vals)):
            #autocor_tau = autocorr.calculate_interal_autocor(error_cal[j], len(error_cal[j]), (j+1))
            #err_val[j], variance_vals[j] = (analys_err.jacknife(error_cal[j], len(error_cal[j]), int(autocor_tau)))
            #variance_vals[j] = analys_err.calculate_variance(error_cal[j])
            #auto_corr.append(autocor_tau)
            #err[j] = np.sqrt(variance_vals[j]*autocor_tau/20)

        #variance, error_an = analys_err.calculate_variance_of_array(corr_fun_vals)
        #print("this is error anaalysis", err)
        
        #plt.plot( corr_fun_dist, auto_corr)
        #plt.show()
        
        effective_mass_lattice = []
        error = []
        for l in range(1, len(corr_fun_vals)-1):
            val_inside = np.mean(corr_fun_vals[l-1])/np.mean(corr_fun_vals[l+1])
            error.append(0.5*np.log(err[l-1]/err[l+1]))
            #print(val_inside)
            effective_mass = 0.5*np.log(val_inside)
            effective_mass_lattice.append(effective_mass)
    
        array_for_change = np.arange(len(effective_mass_lattice))
        plt.plot(array_for_change, effective_mass_lattice)
        plt.errorbar(array_for_change, effective_mass_lattice, error, fmt='none', ecolor='black')
        plt.xlabel(r'$\Delta \tau$')
        plt.ylabel(r'$m_{eff}(\Delta \tau)$')
        plt.show()
        
        """
        bins = np.linspace(1, 20)
        corr = two_point(self.array, bins)
        plt.plot(bins, corr)
        plt.show()
        """
        plt.plot(corr_fun_dist, means, label='correlation fucntion')
        #plt.errorbar(corr_fun_dist, correlation_of_numbers_sum, yerr=error_an, fmt='none', ecolor='purple')
        plt.errorbar(corr_fun_dist, means, yerr=err, fmt='none', ecolor='black', label = 'error bars')
        #plt.errorbar(corr_fun_dist, correlation_of_numbers_sum, yerr=err_arr_cut1, fmt='none', ecolor='yellow')
        plt.ylabel(r"G($\Delta \tau$)")
        plt.xlabel(r'$\Delta \tau$')
        plt.legend()
        plt.title("Two-point correlation vs lattice spacing change")
        plt.show()
        logar = np.log(means)
        theta = np.polyfit(corr_fun_dist[0:20], logar[0:20], 1)
        poly_fn = np.poly1d(theta)  
        print("theta", theta[0])
        #log__2 = np.log(err_arr_cut1)
        lengs = np.log(err)
        plt.scatter(corr_fun_dist, logar)
        plt.plot(corr_fun_dist[0:20], logar[0:20], corr_fun_dist[0:20], poly_fn(corr_fun_dist[0:20]))
        #plt.errorbar(corr_fun_dist, logar, log_err, fmt = 'none', ecolor='pink')
        #plt.errorbar(corr_fun_dist, logar, log_er_2, fmt = 'none', ecolor='black')
        plt.errorbar(corr_fun_dist, logar, lengs, fmt='none', ecolor='black')
        plt.title("The log plot of the two point correlation vs lattice spacing")
        plt.xlabel(r'$\Delta \tau$')
        plt.ylabel(r'log(G($\Delta \tau$))')
        plt.show()
        

        result = open('resultcorr.txt', 'a+')
        arr1 = [[corr_fun_vals, corr_fun_dist, effective_mass_lattice, err, error, means, theta[0]]]
        print(arr1, file= result)
        result.close()