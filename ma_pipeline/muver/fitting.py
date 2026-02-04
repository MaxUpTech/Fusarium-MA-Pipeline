#!/usr/bin/env python3
"""
Fitting Functions for Muver Statistical Models

Provides Gaussian and Logistic fitting functions used throughout
the muver statistical models.
"""

import math
import numpy as np


def gaussian(x, mu, sigma):
    """
    Gaussian (normal) probability density function.
    
    Parameters
    ----------
    x : float or array-like
        Input value(s)
    mu : float
        Mean of the distribution
    sigma : float
        Standard deviation of the distribution
    
    Returns
    -------
    float or array-like
        PDF value(s) at x
    """
    return (1 / math.sqrt(2 * math.pi * sigma ** 2)) * \
        np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def logistic(x, x0, L, M, k):
    """
    Logistic function for fitting indel rates in repeat regions.
    
    The logistic function models the relationship between repeat tract
    length and indel error rates.
    
    Parameters
    ----------
    x : float or array-like
        Input value(s) - typically repeat tract length
    x0 : float
        Midpoint of the logistic curve (inflection point)
    L : float
        Maximum value of the curve (amplitude)
    M : float
        Minimum value of the curve (baseline)
    k : float
        Steepness of the curve (growth rate)
    
    Returns
    -------
    float or array-like
        Logistic function value(s) at x
    
    Notes
    -----
    Formula: M + L / (1 + exp(-k * (x - x0)))
    
    For indel rate correction, rates are log-transformed before fitting,
    so the output represents log10(rate).
    """
    return M + (L / (1 + np.exp(-k * (x - x0))))


def inverse_logistic(y, x0, L, M, k):
    """
    Inverse of the logistic function.
    
    Parameters
    ----------
    y : float
        Output value to invert
    x0, L, M, k : float
        Logistic function parameters
    
    Returns
    -------
    float
        x value corresponding to y
    """
    if y <= M or y >= M + L:
        return None
    return x0 - (1/k) * np.log(L / (y - M) - 1)


def fit_gaussian(data, initial_guess=None):
    """
    Fit a Gaussian distribution to data.
    
    Parameters
    ----------
    data : array-like
        Data to fit
    initial_guess : tuple, optional
        Initial guess for (mu, sigma)
    
    Returns
    -------
    tuple
        (mu, sigma) fitted parameters
    """
    from scipy.stats import norm
    from scipy.optimize import curve_fit
    
    # Initial estimate using scipy.stats
    p0_mu, p0_sigma = norm.fit(data)
    if p0_sigma == 0:
        p0_sigma = 0.01
    
    if initial_guess is not None:
        p0_mu, p0_sigma = initial_guess
    
    # Create histogram for fitting
    hist, bin_edges = np.histogram(data, bins='auto', density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    try:
        popt, pcov = curve_fit(gaussian, bin_centers, hist, 
                               p0=[p0_mu, p0_sigma], maxfev=100000)
        mu, sigma = popt
        sigma = abs(sigma)
        return mu, sigma
    except:
        return p0_mu, abs(p0_sigma)


def fit_logistic(x_data, y_data, initial_guess=None):
    """
    Fit a logistic function to data.
    
    Parameters
    ----------
    x_data : array-like
        Independent variable data
    y_data : array-like
        Dependent variable data
    initial_guess : tuple, optional
        Initial guess for (x0, L, M, k)
    
    Returns
    -------
    dict
        Dictionary with fitted parameters: x0, L, M, k
    """
    from scipy.optimize import curve_fit
    
    if len(x_data) < 4:
        return None
    
    # Generate initial guesses if not provided
    if initial_guess is None:
        p0_k = 0.1
        p0_L = max(y_data) - min(y_data)
        p0_M = min(y_data)
        
        # Find midpoint
        mid_value = (max(y_data) + min(y_data)) / 2
        mid_diff = float('inf')
        peak_idx = np.argmax(y_data)
        p0_x0 = x_data[0]
        
        for i, val in enumerate(y_data):
            diff = abs(val - mid_value)
            if diff < mid_diff and i <= peak_idx:
                mid_diff = diff
                p0_x0 = x_data[i]
        
        initial_guess = (p0_x0, p0_L, p0_M, p0_k)
    
    try:
        popt, pcov = curve_fit(
            logistic,
            x_data,
            y_data,
            p0=initial_guess,
            maxfev=1000000,
            method='trf',
        )
        x0, L, M, k = popt
        
        return {
            'x0': x0,
            'L': L,
            'M': M,
            'k': k,
        }
    except:
        return None
