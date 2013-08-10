#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
from detector import Detector


# global variable
win_width = 300  # window size
win_slide = 100  # window-migration-length in 1 step
minimum_length = 10  # minimum term that we regard the detected term as outburst
maximum_length = 1000  # maximum term that we regard the detected term as outburst
sigma_level = 2  # set for decide standard line
smoothing_level = 2  # what time we will smooth the time-series-data(TSD)
wma_width = 7  # [wma_width] points weighted moving average



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
    
        
