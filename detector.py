#!/usr/bin/env python
#-*- coding: utf-8 -*-

#import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import os

import timeseries


# global variable
margin = 0  # additional term for both ends of detected outbursts
win_width = 300  # window size
win_slide = 100  # window-migration-length in 1 step
minimum_length = 10  # minimum term that we regard the detected term as outburst
maximum_length = 1000  # maximum term that we regard the detected term as outburst
sigma_level = 2  # stddev ratio for threshold
smoothing_level = 2  # repeats of smoothing
wma_width = 7  # [wma_width] points weighted moving average


class Detector:

    def __init__(self, src):

        # read a datafile
        ts = timeseries.fromfile(src)

        # normalize, interpolate, smoothing
        ts.normalize()
        ts.interpolate()
        ts.smooth(smoothing_level, wma_width)

        # target to detect outlier
        self.ts = ts
        
        self.detect(ts.rate, margin, win_width, win_slide, minimum_length, maximum_length, sigma_level)


    def detect(self, seq, width, slide, ratio, min_len, max_len, margin):
        self.windows = self.get_windows(seq, width, slide, ratio)
        burst_h, burst_l = self.merge_window(windows)
        burst_h = filter_burst(burst_h, min_len, max_len, margin)
        burst_l = filter_burst(burst_l, min_len, max_len, margin)
        
        
        
    def get_windows(self, seq, width, slide, ratio):
        """Detect outbursts in a TimeSeries.

        Args:
            seq: sequence to detect bursts
            width: window size
            slide: window-migration-length in 1 step
            ratio: std ratio for threshold

        Returns:
            windows: list of Window objects
        """
        windows = []

        head = 0  # head index of window in seq
        tail = head + width
        while head < len(seq):
            median = np.median(seq[head : tail])
            threshold = np.std(seq[head : tail]) * ratio
            burst_h = [x if x > threshold_h else np.nan for x in seq[head : tail]]
            burst_l = [x if x > threshold_l else np.nan for x in seq[head : tail]]

            windows.append(Window(head, tail, median, threshold, burst_h, burst_l))
            head += slide
            tail += slide
            
        return windows
            
            
    def merge_window(self, windows):
        """Extract and merge bursts from Window objects.

        Args:
            windows: list of Window objects
        
        Returns:
            burst_h: list of rate, only with higher bursts
            burst_l: list of rate, only with lower bursts
        """
        burst_h = []
        burst_l = []
        index = 0
        for w in windows:
            diff = index - w.head
            burst_h += w.burst_h[index : w.tail]
            burst_l += w.burst_l[index : w.tail]
            index = w.tail
        return (burst_h, burst_l)


    def get_burst_indexes(self, burst):
        """Get head/tail indexes of the continuous bursts."""
        indexes = []
        i = 0
        while i < len(burst):
            if burst[i] != np.nan:
                j = i
                while burst[j] != np.nan:
                    j += 1
                    if j == len(burst):
                        break
                indexes.append((i, j))
                i = j
            else:
                i += 1
        return indexes
        

    def filter_burst(self, seq, min_len, max_len):
        """Remove bursts with too long/short length.

        Args:
            seq: seq of bursts
            min_len: minimum length to regard outliers term as an outburst
            max_len: maximum ...
        
        Returns:
            filtered: filtered indexes of bursts
        """
        indexes = get_burst_indexes(burst)
        filtered = []

        last = 0
        for i in indexes:
            filtered += [np.nan] * (i[0] - last)
            if min_len <= i[1] - i[0] <= max_len:
                filtered += seq[i[0] : i[1]]
            else:
                filtered += [np.nan] * (i[1] - i[0])
            last = i[1]
        filtered += [np.nan] * (len(seq) - last)
        return filtered

        

    def plot(self, title, dst, width, slide, sigma):
        """plot bursts with matplotlib."""

        path = dst + "/" + title
        if os.path.exists(path) == False:
            os.mkdir(path)

        for w in self.windows:
            tmp = copy.copy(ts)
            tmp.rate = ([np.nan] * w.index +
                        tmp.rate[w.head : w.tail] +
                        [np.nan] * (len(tmp.rate) - w.tail))
            self.plot_window
            

        self.plot_burst(ts, title, path + "/total.png")
            

    def plot_burst(self, ts, title, name):
        """plot 1 burst with matplotlib.

        Args:
            ts: original TimeSeries
            title: figure title
            name: name for PNG file
        """
        
        plt.clf()
        #plt.xticks([self.time.index(x) for x in xtick], tuple([str(x) for x in xtick]), rotation = 30)
        plt.xlim([0, len(ts.time) - 1])
        plt.ylim([min(ts.rate), max(ts.rate)])
        plt.plot(ts.rate, 'c-')
        plt.plot(ts.smoothed, 'g-')
        plt.plot(ts.burst, 'r-', lw = 2)

        if dst == None:
            plt.show()
        else:
            plt.savefig(path + "/step" + str(i) + ".png")

        
    def plot_window(self, window):
        start_edge = slide * i
        if slide * i + width < len(self.time):
            end_edge = start_edge + width
        else:
            end_edge =len(self.time)
            
        plt.plot([start_edge, start_edge], [-10, 10], 'k-', lw = 2)
        plt.plot([end_edge, end_edge], [-10, 10], 'k-', lw = 2)
        plt.plot([self.median[i] + self.std[i] * sigma] * len(self.rate), 'r--')
        plt.plot([self.median[i]] * len(self.rate), 'k--')        
    


