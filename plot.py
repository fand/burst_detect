#!/usr/bin/env python
#-*- coding: utf-8 -*-

# ver.20130726(hayashi)
#   remove small bugs
#   remove never-used elements
#   modify the fomula for outbursts evaluation
#     (length * average * standard diviation -> length * diff average)
#   put out result as CSV file
#   use only WMA in smoothing
#
# ver.20130725(hayashi)
#   available for "BAT" data
#   interpolate lost data points
#   make plot figure for overall outbursts
#   evaluate for detected outbursts


#import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import os

# global variable
win_width = 300  # window size
win_slide = 100  # window-migration-length in 1 step
minimum_length = 10  # minimum term that we regard the detected term as outburst
maximum_length = 1000  # maximum term that we regard the detected term as outburst
sigma_level = 2  # set for decide standard line
smoothing_level = 2  # what time we will smooth the time-series-data(TSD)
wma_width = 7  # [wma_width] points weighted moving average


class Detector:

    # make structure: "self"
    def __init__(self, src):

        global win_width
        global win_slide
        global minimum_length
        global sigma_level
        global smoothing_level
        global wma_width

        # read a datafile
        self.readTxt(src)

        # normalization
        self.normalize(self.origin_rate, self.origin_error)

        # interpolate falled out data point
#        self.interpolate(self.origin_time, self.origin_rate, self.origin_error)
        self.interpolate(self.origin_time, self.normal_rate, self.normal_error)
        
        # smoothing
        self.smoothed = self.smooth(self.rate, self.error, smoothing_level)

        # target to detect outlier
        target = self.smoothed
        
        self.detect(target, 300, win_width, win_slide, minimum_length, maximum_length, sigma_level)  # detect outbursts


    # read a datafile        
    def readTxt(self, src):
        with open(src) as f:
            lines = filter(lambda l: l[0] != '#', f.readlines()[:-1])

        zipped = map(lambda l: map(float, l.strip().split()[:3]), lines)
        self.origin_time, self.origin_rate, self.origin_error = map(list, zip(*zipped))

    # normalization (ignore minus-value- point)
    def normalize(self, orate, oerror):

        nrate = []
        nerror = []
        max_rate = np.max(orate)

        for i in range(0, len(orate)):
            nrate.append(orate[i] / max_rate)
            nerror.append(oerror[i]/ max_rate)

        self.normal_rate = nrate
        self.normal_error = nerror                        

    # interpolate lost data point (liner interpolation)
    def interpolate(self, otime, orate, oerror):

        itime = []
        irate = []
        ierror = []
        new_idx = 0

        for i in range(0, len(otime) - 1):
            time_diff = otime[i + 1] - otime[i] - 1
            itime.append(otime[i])
            irate.append(orate[i])
            ierror.append(oerror[i])
            if time_diff != 0.0:
                for j in range(1, int(time_diff) + 1):
                    itime.append(otime[i] + float(j))
                    rate_variation = (orate[i + 1] - orate[i]) / (time_diff + 1)
                    error_variation = (oerror[i + 1] - oerror[i]) / (time_diff + 1)
                    irate.append(orate[i] + rate_variation * j)
                    ierror.append(oerror[i] + error_variation * j)
                new_idx = len(itime)
            else:
                new_idx += 1

        self.time = itime
        self.rate = irate
        self.error = ierror

    
    # smoothing(in [smoothness] time)
    def smooth(self, src, err, smoothness):

        global wma_width
        s = copy.copy(src)

        for i in range(smoothness):
            s = self.wma_err(s, err, wma_width)
        
        return s
    

    # Weighted Moving Average([n] points, use error breadth)
    #     [n] is odd number
    def wma_err(self, src, err, n):
        s = copy.copy(src)
        n0 = (int)(n / 2)

        # both ends
        for i in range(1, n0 - 1):
            l_end = (src[ i] / err[ i]) * float(n0 + 1)
            r_end = src[-i] / err[-1] 
            l_weightsum = 1.0 / err[ i]
            r_weightsum = 1.0 / err[-1]
            for j in range(1, i + 1):
                l_end += ((src[i - j] / err[i - j]) + (src[i + j] / err[i + j])) * float(n0 + 1 - j)
                r_end += ((src[-(i - j + 1)] / err[-(i - j + 1)]) + (src[-(i + j + 1)] / err[-(i + j + 1)])) * float(n0 + 1 - j)
                l_weightsum += (1.0 / err[i - j]) + (1.0 / err[i + j]) * float(n0 + 1 - j)
                l_weightsum += (1.0 / err[-(i - j + 1)]) + (1.0 / err[-(i + j + 1)]) * float(n0 + 1 - j)
            s[ i] = l_end / l_weightsum
            s[-i] = r_end / r_weightsum

        # except ends
        for i in range(n0, len(src) - n0):
            t = src[i] / err[i]
            weightsum = (1.0 / err[i])
            for j in range(1, n0):
                t += (src[i - j] / err[i - j]) + (src[i + j] / err[i + j]) 
                weightsum += (1.0 / err[i - j]) + (1.0 / err[i + j])
            s[i] = t / weightsum

        return s
    

    # outburst detection
    #     src: smoothed source
    #     margin: ? (now, no use)
    #     width: window size
    #     slide: window-migration-length in 1 step
    #     min_len: minimum term that we regard the detected term as outburst 
    #     max_len: maximum term that we regard the detected term as outburst 
    #     sigma: set for decide standard line (red broken line)
    def detect(self, src, margin, width, slide, min_len, max_len, sigma):
        se = []  # standard error for each step [float list]
        md = []  # median for each step [float list]
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

            se.append(np.std(self.error[start_pnt : end_pnt]))
            md.append(np.median(src[start_pnt : end_pnt]))

            burst_day_h_step = []  # detected outburst day for each step [int list]

            for i in range(start_pnt, end_pnt):
                if src[i] > md[step] + se[step] * sigma:
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

        self.std_err = se
        self.median = md
        self.bst_h = burst_h
        self.bst_term = burst_term_h
        self.bst_all = burst_all
        self.bst_term_all = burst_term_all
        self.bst_val = burst_value

        step = 0


    # plot
    def plot(self, title, dst, width, slide, sigma):

        path = dst + "/" + title
        if os.path.exists(path) == False:
            os.mkdir(path)

        # plot for each step
        for i in range(0, len(self.std_err)):
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

            plt.plot([self.median[i] + self.std_err[i] * sigma] * len(self.rate), 'r--')
            plt.plot([self.median[i]] * len(self.rate), 'k--')        

            plt.suptitle(title)

            if dst == None:
                plt.show()
            else:
                plt.savefig(path + "/step" + str(i) + ".png")

        # plot all outburst
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

        
if __name__=='__main__':

    global win_width
    global win_slide
    global sigma_level

    if len(sys.argv) != 3:
        print("Usage: # python %s source_file dst_dir" % sys.argv[0])
        print("%s" % sys.argv[2])
        quit()

    d = Detector(sys.argv[1])

    # result (standard) output
    print("\nMedian + Standard error * %.1f in each step:" % sigma_level)
    for i in range(0, len(d.std_err)):
        line = d.median[i] + d.std_err[i]
        print("  step %2d: %f" % (i, line))

    print("\nDetected outburst in each step:")
    for i in range(0, len(d.std_err)):
        if i < len(d.bst_term):
            for j in range(0, len(d.bst_term[i])):
                print("  step %2d: %.1f-%.1f [length: %d]" % (i, d.time[d.bst_term[i][j][0]], d.time[d.bst_term[i][j][1]], int(d.bst_term[i][j][1] - d.bst_term[i][j][0] + 1.0)))

    print("\nDetected outburst (all):")
    for i in range(0, len(d.bst_term_all)):
        print("  %2d: %.1f-%.1f [length: %3d]   value: %f" % (i + 1, d.time[d.bst_term_all[i][0]], d.time[d.bst_term_all[i][1]], int(d.bst_term_all[i][1] - d.bst_term_all[i][0] + 1.0), d.bst_val[i]))
    
    # plot result
    d.plot(sys.argv[1][8:-4], sys.argv[2], win_width, win_slide, sigma_level)

    # write result to CSV file
    f = open(sys.argv[2] + "/" + sys.argv[1][8:-4] + "/" + "result.csv", "w")
    f.write(sys.argv[1][8:-4] + "\n")
    f.write("Observartion Period: %.1f - %.1f\n" % (d.time[0], d.time[-1]))
    f.write("No,Start,End,Length,Value\n")
    for i in range(0, len(d.bst_term_all)):
        f.write("%d,%f,%f,%d,%f\n" % (i + 1, d.time[d.bst_term_all[i][0]], d.time[d.bst_term_all[i][1]], int(d.bst_term_all[i][1] - d.bst_term_all[i][0] + 1.0), d.bst_val[i]))
