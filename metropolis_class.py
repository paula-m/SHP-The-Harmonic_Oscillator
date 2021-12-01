import numpy as np
import math
import copy

class metropolis:
    def __init__(self, path, h, m, w, path_len_N_t, lamb, z):
        """this initialises the class and the values of m, w, the initial path, the length and h value"""
        self.path = path
        self.m = m
        self.w = w
        self.nt = path_len_N_t
        self.h = h
        self.lamb = lamb
        self.Z = z
        print("this is self m", self.m)

    def define_action(self, x_next, x_i):
        """This is the function which defines the action used for the path integrals"""
        return 0.5*self.m*(x_next - x_i)**2 + 0.5*self.m*(self.w**2)*(x_i**2)

    def anharmonic_oscillator_action(self, x_next, x_i, x_prev):
        return 0.5*self.m*(x_next - x_i)**2 + 0.5*self.m*(self.w**2)*(x_i**2) + self.lamb*(x_i**2 -1)**2 + self.Z*(x_i - x_prev)**2

    def calculate_the_action(self, new_path, chng_i):
        """the method function which calculates the action of the path based on the defined action before"""
        S = 0
        for i in range(chng_i-2, chng_i + 1):   #looks only at indices of the previous, changed and next x values
        #for i in range(self.nt):
            x_val_next = new_path[(i+1) % self.nt]  #this modifies the path given and calculates the next value of x
            x_val_prev = new_path[(i-1)%self.nt]
            #S += self.anharmonic_oscillator_action(x_val_next, new_path[i], x_val_prev)    #finds the contribution to the action S and adds it to the change in action
            S += self.define_action(x_val_next, new_path[i])
        return S

    def metropolis_sweep(self, path):
        """the method function which creates the metropolis sweep every single time. """
        new_path = copy.deepcopy(path)
        counter = 0     #for calculating the acceptance rate
        for i in range(self.nt):
            trial_path = new_path
            calc_s_initial = self.calculate_the_action(new_path, i) #calculates initial action
            small_x_change = np.random.uniform(low=-self.h, high=self.h)    #random uniform number to find the change in x
            trial_path[i] = trial_path[i] + small_x_change  #creates trial x
            calc_s_change = self.calculate_the_action(trial_path, i)    #calculates the new path
            delta_s = calc_s_change - calc_s_initial #change in path
            if math.exp(-delta_s) > np.random.uniform(low=0, high=1):   #acceptance conditions
                new_path[i] = trial_path[i]     #the accepted value
                counter += 1/self.nt    #acceptance rate
            else:
                new_path[i] = path[i]
        #self.h = self.h*(counter/0.6)   #adjusts the h after the sweep to keep acceptance rate around 70%
        print("this is the acceptance rate", counter)
        return new_path


    def running_the_sweep(self, num_runs):
        """the function for rerunning the sweep multiple times and returning the array of paths."""
        paths = []
        for i in range(num_runs):
            self.path = metropolis.metropolis_sweep(self, self.path)
            if i%5 == 0:
                paths.append(self.path)
        return paths
