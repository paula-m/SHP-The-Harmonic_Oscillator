from metropolis_class import metropolis
import numpy as np
import matplotlib.pyplot as plt


class thermalisation:
    def __init__(self, nt, h, m, w, run_no, lam, z):
        self.nt = nt
        self.h = h
        self.m = m
        self.w = w
        self.r = run_no
        self.lam = lam
        self.z = z

    def create_empty_array(self):
        return np.zeros(self.nt)

    def create_hot_start(self):
        return np.random.randn(self.nt)

    def hot_start(self):
        starting_array = self.create_empty_array()
        for i in range(self.nt):
            starting_array[i] = np.random.uniform(low=-self.h, high=self.h)
        return starting_array

    def thermalised_start(self):
        path = self.hot_start()
        metro = metropolis(path, self.h, self.m, self.w, self.nt, self.lam, self.z)
        array_path = metro.running_the_sweep(self.r)
        return array_path[-1]