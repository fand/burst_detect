#!/usr/bin/env python
#-*- coding: utf-8 -*-

#import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import os

from timeseries import TimeSeries


# global variable
win_width = 300  # window size
win_slide = 100  # window-migration-length in 1 step
minimum_length = 10  # minimum term that we regard the detected term as outburst
maximum_length = 1000  # maximum term that we regard the detected term as outburst
sigma_level = 2  # set for decide standard line
smoothing_level = 2  # what time we will smooth the time-series-data(TSD)
wma_width = 7  # [wma_width] points weighted moving average


class Detector:

    def __init__(self, src):

        global win_width
        global win_slide
        global minimum_length
        global sigma_level
        global smoothing_level
        global wma_width

        # read a datafile
        self.ts = TimeSeries(src)

        # normalize, interpolate, smoothing
        self.ts.normalize()
        self.ts.interpolate()
        self.ts_smoothed = self.ts.smooth(smoothing_level, wma_width)

        # target to detect outlier
        target = self.ts_smoothed
        
        self.detect(target, 300, win_width, win_slide, minimum_length, maximum_length, sigma_level)  # detect outbursts


    def detect(self, src, margin, width, slide, min_len, max_len, sigma):
        """Detect outbursts in a TimeSeries.

        Args:
            src: smoothed source
            margin: ? (now, no use)
            width: window size
            slide: window-migration-length in 1 step
            min_len: minimum term that we regard the detected term as outburst
            max_len: maximum term that we regard the detected term as outburst
            sigma: set for decide standard line (red broken line)

        Returns:
            a
        """

        std = []  # standard deviation for each step [float list]
        median = []  # median for eacph step [float list]
        burst_day_h = []  # days of high bursts [int list]
        burst_day_l = []  # days of low bursts [int list] (now, no use)
        burst_term_h = []  # terms of high bursts [(float, float) list]
        event_term_burst_l = []  # terms of low bursts [(float, float) list] (now, no use)
        burst_h = []  # extract only outburst from original data [float list]

        # window position
        start_pnt = 0
        end_pnt = width - 1

        # what time the window has slided
        step = 0

        while end_pnt < len(src) + slide - 1:

            if end_pnt >= len(src):
                end_pnt = len(src) - 1

            std.append(np.std(self.error[start_pnt : end_pnt]))
            md.append(np.median(src[start_pnt : end_pnt]))

            burst_day_h_step = []  # detected outburst day for each step [int list]

            for i in range(start_pnt, end_pnt):
                if src[i] > median[step] + std[step] * sigma:
                    burst_day_h_step.append(i)

            burst_day_h.append(burst_day_h_step)

            start_pnt += slide
            end_pnt += slide

            # remove the detect term that length is shorter than [min_len]
            burst_term_h_step = []  # detected outburst term for each step [(float, float) list]
            burst_flag = True
            start_day = 0 if len(burst_day_h[step]) == 0 else burst_day_h[step][0]
            start_day = 0 if len(burst_day_h[step]) == 0 else burst_day_h[step][0]

            for i in range(0, len(burst_day_h[step]) - 1):

                # indexes of detected points are continuous ?
                #     Yes: continue getting next index
                #     No : judge the term is a outburst
                if burst_day_h[step][i + 1] - burst_day_h[step][i] == 1 and i + 1 != len(burst_day_h[step]) - 1:
                    if burst_flag == False:
                        start_day = burst_day_h[step][i]
                        burst_flag = True
                else:
                    # [burst_flag] = True ?
                    #     Yes: judge the term length
                    #     No : delete [i]th point
                    if burst_flag == True:
                        if burst_day_h[step][i + 1] - burst_day_h[step][i] == 1 and i + 1 != len(burst_day_h[step]) - 1:
                            end_day = burst_day_h[step][i + 1]
                        else:
                            end_day = burst_day_h[step][i]

                        # continuous term is longer than [min_len] and shorter than [max_len] ?
                        #     Yes: add the term to <burst_term_h_step>
                        #     No : detele data points in the term (substitute -1)
                        if end_day - start_day + 1 >= min_len and end_day - start_day + 1 <= max_len:
                            burst_term_h_step.append((start_day, end_day))
                        else:
                            idx_s = burst_day_h[step].index(start_day)
                            idx_e = burst_day_h[step].index(end_day) + 1
                            for j in range(idx_s, idx_e + 1):
                                burst_day_h[step][j] = -1

                        burst_flag = False

                    else:
                        burst_day_h[step][i] = -1
                        if i + 1 == len(burst_day_h[step]) - 1:
                            burst_day_h[step][i + 1] = -1

            burst_term_h.append(burst_term_h_step)

            # extract outburst form original data
            burst_h_step = []
            for i in range(0, len(src)):
                if i in burst_day_h[step]:
                    burst_h_step.append(self.rate[i])
                else:
                    burst_h_step.append(float('nan'))
            burst_h.append(burst_h_step)

            step += 1
            # end: [while end_pnt < len(src) + slide - 1:]

        # merge extracted outbursts in all steps
        burst_all = copy.copy(burst_h[0])
        for i in range(1, len(burst_h)):
            for j in range(0, len(burst_h[i])):
                if np.isnan(burst_h[i][j]) == False:
                    burst_all[j] = burst_h[i][j]

        # terms of all outbursts
        burst_term_all = []
        burst_flag = False
        start_day = 0
        end_day = 0
        for i in range(0, len(burst_all)):
            if np.isnan(burst_all[i]) == False and i != len(burst_all) - 1:
                if burst_flag == False:
                    start_day = i
                    burst_flag = True
            else:
                if burst_flag == True or i == len(burst_all) - 1:
                    end_day = i
                    if burst_flag == True and i != len(burst_all) - 1:
                        if end_day - start_day + 1 >= min_len and end_day - start_day + 1 <= max_len:
                            burst_term_all.append((start_day, end_day))
                    burst_flag = False

        # evaluate each outburst
        # <valuation basis (temporary)>
        #   * a longer outburst has higher value
        #   * higher average that higher than it in neighber terms makes value higher
        #     + weighted average using error value (weight: 1/error)

        print("\nvalue = Length * Diff Average / 10.0")

        burst_value = []
        for i in range(0, len(burst_term_all)):
            burst = copy.copy(burst_all)[burst_term_all[i][0] : burst_term_all[i][1] + 1]
            lgth = float(len(burst))  # length
            weight = 0.0
            average = 0.0
            r_average = 0.0
            l_average = 0.0

            for j in range(0, len(burst)):
                if np.isnan(burst[j]):
                    burst[j] = 0
                average += burst[j] / self.error[burst_term_all[i][0] + j]
                weight += 1.0 / self.error[burst_term_all[i][0] + j]
            average /= weight

            # average in right term
            if len(self.time[: burst_term_all[i][0]]) < len(burst):
                r_term = len(self.time[: burst_term_all[i][0]])
            else:
                r_term = len(burst)
            for j in range(0, r_term):
                r_average += self.rate[burst_term_all[i][0] - r_term + j] / self.error[burst_term_all[i][0] + j]
                weight += 1.0 / self.error[burst_term_all[i][0] -r_term + j]
            r_average /= weight

            # average in left term
            if len(self.time[burst_term_all[i][1] :]) < len(burst):
                l_term = len(self.time[burst_term_all[i][1] :])
            else:
                l_term = len(burst)
            for j in range(0, l_term):
                l_average += self.rate[burst_term_all[i][0] - l_term + j] / self.error[burst_term_all[i][0] + j]
                weight += 1.0 / self.error[burst_term_all[i][0] -l_term + j]
            l_average /= weight

            diff_average = ((average - r_average) + (average - l_average)) / 2.0
            value = diff_average * lgth / 10.0
            print("        burst %2d: %f * %f = %f" % (i + 1, lgth, diff_average, value))
            burst_value.append(value)

        self.std = td
        self.median = median
        self.bst_h = burst_h
        self.bst_term = burst_term_h
        self.bst_all = burst_all
        self.bst_term_all = burst_term_all
        self.bst_val = burst_value


    def plot(self, title, dst, width, slide, sigma):
        """plot bursts with matplotlib."""
        
        path = dst + "/" + title
        if os.path.exists(path) == False:
            os.mkdir(path)

        # plot bursts for every sliding window
        for i in range(0, len(self.std)):
            plt.clf()

            plt.xlim([0, len(self.time) - 1])
            plt.ylim([np.min(self.rate), np.max(self.rate)])

            plt.plot(self.rate, 'c-')
            plt.plot(self.smoothed, 'g-')
            plt.plot(self.bst_h[i], 'r-', lw = 2)

            start_edge = slide * i
            end_edge = slide * i + width if slide * i + width < len(self.time) else len(self.time) - 1
            plt.plot([start_edge, start_edge], [-10, 10], 'k-', lw = 2)
            plt.plot([end_edge, end_edge], [-10, 10], 'k-', lw = 2)

            plt.plot([self.median[i] + self.std[i] * sigma] * len(self.rate), 'r--')
            plt.plot([self.median[i]] * len(self.rate), 'k--')        

            plt.suptitle(title)

            if dst == None:
                plt.show()
            else:
                plt.savefig(path + "/step" + str(i) + ".png")

        # plot all outbursts in 1 figure
        plt.clf()
        plt.xlim([0, len(self.time) - 1])
        plt.ylim([np.min(self.rate), np.max(self.rate)])
        plt.plot(self.rate, 'c-')
        plt.plot(self.smoothed, 'g-')
        plt.plot(self.bst_all, 'r-', lw = 2)
        plt.suptitle(title)
        if dst == None:
            plt.show()
        else:
            plt.savefig(path + "/overall.png")
            

    def plot_burst(self, ts, ts_smoothed, ts_burst, title, xlim, ylim):
        """plot 1 burst with matplotlib."""
        
        plt.clf()
        plt.xlim([0, len(self.time) - 1])
        plt.ylim([np.min(self.rate), np.max(self.rate)])
        plt.plot(self.rate, 'c-')
        plt.plot(self.smoothed, 'g-')
        plt.plot(self.bst_h[i], 'r-', lw = 2)

        if dst == None:
            plt.show()
        else:
            plt.savefig(path + "/step" + str(i) + ".png")

        
    def plot_window(self, start, end):
        start_edge = slide * i
        if slide * i + width < len(self.time):
            end_edge = start_edge + width
        else:
            end_edge =len(self.time)
            
        plt.plot([start_edge, start_edge], [-10, 10], 'k-', lw = 2)
        plt.plot([end_edge, end_edge], [-10, 10], 'k-', lw = 2)
        plt.plot([self.median[i] + self.std[i] * sigma] * len(self.rate), 'r--')
        plt.plot([self.median[i]] * len(self.rate), 'k--')        
        plt.suptitle(title)
    
    
        # plot all outbursts in 1 figure
        plt.clf()
        plt.xlim([0, len(self.time) - 1])
        plt.ylim([np.min(self.rate), np.max(self.rate)])
        plt.plot(self.rate, 'c-')
        plt.plot(self.smoothed, 'g-')
        plt.plot(self.bst_all, 'r-', lw = 2)
        plt.suptitle(title)
        if dst == None:
            plt.show()
        else:
            plt.savefig(path + "/overall.png")
    
            



