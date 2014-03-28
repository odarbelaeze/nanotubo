# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import numpy as np

class Sample(object):
    """docstring for Sample"""
    def __init__(self, width, height, mu = 10.0, sigma = 1.0):
        super(Sample, self).__init__()
        self.width = width
        self.height = height
        self.mu = mu
        self.sigma = sigma
        self.count = int(np.floor(0.9 * float(width ** 2) * height / ((4.0 / 3) * np.pi * (mu ** 3))))

        self.rads = np.sort(np.random.normal(self.mu, self.sigma, self.count))[::-1]
        self.centers = []

        self.find_centers()

    def _pass_constraints(self, center):
        for prev in self.centers:
            if np.sqrt(np.sum((prev[1] - center[1]) ** 2)) < prev[0] + center[0]:
                return False
        return True

    def find_centers(self):
        for rad in self.rads:
            while True:

                x = np.random.uniform(self.mu, self.width - self.mu)
                y = np.random.uniform(self.mu, self.width - self.mu)
                z = np.random.uniform(self.mu, self.height - self.mu)

                center = [rad, np.array([x, y, z])]
                if self._pass_constraints(center):
                    self.centers.append(center)
                    print len(self.centers), 'of', self.count, center
                    break
        

def main():
    sample = Sample(100.0, 10.0)
    print sample.rads

if __name__ == '__main__':
    main()
