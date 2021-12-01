import numpy as np

class error_calculation():
    def __init__(self) -> None:
        pass

    def calculate_the_mean_element(self, array):
        mean = np.mean(array, axis = 0)
        return mean
    
    def calculate_mean_markov(self, array):
        mean = np.mean(array)
        return mean
    
    def calculate_variance(self, array):
        return np.var(array, axis=0, ddof=1)

    def calculate_error_markov_chain_naive(self, array):
        variance = np.var(array, ddof=1)
        print("variance", variance)
        print(len(array))
        return np.sqrt((variance/(len(array))))

    def calculate_error(self,array,  no_samples, autocorr):
        return np.sqrt((self.calculate_variance(array)/no_samples)*autocorr)

    def calculate_naive_error(self,array, no_samples):
        return np.sqrt(self.calculate_variance(array)/no_samples)

    def calculate_power(self, array, power):
        squared_array = np.copy(array)**power
        return squared_array

    def jacknife(self, array, N, B):
        new_arr = np.copy(array)
        #print("this is N/B", N%B)
        if (N%B) != 0:
            a = N%B
            for i in range(a):
                new_arr = np.delete(new_arr, 0, 0)
        estimator = np.ones(len(new_arr))
        N = len(new_arr)
        #print(len(new_arr))
        #print(N)
        #return estimator

        sums = np.sum(new_arr)
        mean = sums/len(new_arr)
        #mean_ar = np.ones(len(estimator))*mean
        #print(mean_ar)
        for i in range(int((N/B))):
            #j = i*B
            elim = 0
            for j in range(i*B, (i*B +B)):
                elim += new_arr[j]
            estimator[i] = (sums-elim)/(N-B)

        variance_j = 0
        for k in range(int(N/B)):
            variance_j += (((N/B) -1)/(N/B))*(estimator[k] - mean)**2
        error_j = np.sqrt(variance_j)
        return error_j, variance_j

    def jackknife_analysis(self, array, N, B):
        new_arr = np.copy(array)

        estimate_o = np.mean(new_arr)

        if (N%B) != 0:
            a = N%B
            for p in range(a):
                new_arr = np.delete(new_arr, 0, 0)
        N = len(new_arr)
        N_b = int(N/B)
        block_est = np.ones((N_b))
        for k in range(N_b):
            ok = 0
            for i in range(B):
                ok += new_arr[(k-1)*(B) + i]
            block_est[k] = (1/B)*ok
        
        j_block_est = np.ones(N_b)
        for j in range(N_b):
            O_i = np.sum(new_arr)
            j_block_est[j] = (1/(N-B))*(O_i - B*block_est[j])
        
        sum_j = 0
        for l in range(N_b):
            sum_j += (j_block_est[l] - estimate_o)**2

        stand_dev = ((N_b -1)/(N_b))*sum_j
        return np.sqrt(stand_dev)
