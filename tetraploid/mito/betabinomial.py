#! /usr/bin/env python
import numpy as np
import scipy
from scipy.optimize import minimize
import json

class BetaBinomialDistribution():
    def __init__(self, alpha, beta):
        self.alpha = alpha
        self.beta = beta
        self.parameters = (self.alpha, self.beta)
        self.summaries = np.zeros(3)
        self.d = 2

    def probability(self, X):
        if self.alpha <= 0 or self.beta <= 0:
            # Return NA value
            return np.full(X.shape, 0)
            #return 0
        return np.exp(self.log_probability(X))

    def log_probability(self, X):
        if self.alpha <= 0 or self.beta <= 0:
            pass
        return scipy.stats.betabinom.logpmf(X[:,0], X[:,1], self.alpha, self.beta)
    
    def neg_ll(self, params, data, w):
        alpha, beta = params
        lls = scipy.stats.betabinom.logpmf(data[:,0], data[:,1], alpha, beta)
        return -np.sum(w*lls)

    def summarize(self, X, w=None):
        if self.alpha <= 0 or self.beta <= 0:
            return

        if w is None:
            w = np.ones(X.shape[0])
        bounds = [(0, None), (0, None)]
        
        result = minimize(self.neg_ll, [self.alpha, self.beta], args=(X,w), method='L-BFGS-B', bounds=bounds)
        alpha_mle, beta_mle = result.x
        self.summaries[0] = alpha_mle
        self.summaries[1] = beta_mle

    def from_summaries(self, inertia=0.0):
        if self.alpha <= 0 or self.beta <= 0:
            return

        if self.summaries[0] == 0:
            self.alpha = -1
            self.beta = -1
            return
        
        self.alpha = self.summaries[0]
        self.beta = self.summaries[1]
        

    def clear_summaries(self, inertia=0.0):
        self.summaries = np.zeros(3)

    def to_dict(self):
        self_dict = {}
        self_dict['type'] = 'BetaBinomialDistribution'
        self_dict['parameters'] = (self.alpha, self.beta)
        return self_dict

    def __repr__(self):
        return json.dumps(self.to_dict(), indent=4)

    @classmethod
    def from_samples(cls, X, weights=None):
        d = BetaBinomialDistribution(0, 0)
        d.summarize(X, weights)
        d.from_summaries()
        return d

    @classmethod
    def blank(cls):
        return BetaBinomialDistribution(0, 0)


