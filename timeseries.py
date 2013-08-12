import numpy as np
import copy


def _fromFile(self, src):
    """Create a TimeSeries object from file."""
    with open(src) as f:
        lines = filter(lambda l: l[0] != '#', f.readlines()[:-1])
    zipped = map(lambda l: map(float, l.strip().split()[:3]), lines)
    s = map(list, zip(*zipped)))
    return TimeSeries(s[0], s[1], s[2])


class TimeSeries:

    def __init__(self, time, rate, error, smoothed = None, burst = None):
        self.time = time
        self.rate = ratio
        self.error = error
        self.smoothed = smoothed if smoothed else [np.nan] * len(time)
        self.burst = burst if burst else [np.nan] * len(time)        

        
    def normalize(self):
        """normalize self.rate and self.error. (ignore negative-values)"""
        nrate = []
        nerror = []
        max_rate = np.max(self.rate)

        for i in range(0, len(self.rate)):
            nrate.append(self.rate[i] / max_rate)
            nerror.append(self.error[i]/ max_rate)

        self.rate = nrate
        self.error = nerror                        

        
    def interpolate(self):
        """interpolate lost data point (liner interpolation)."""
        itime = []
        irate = []
        ierror = []
        new_idx = 0

        for i in range(0, len(self.time) - 1):
            time_diff = self.time[i + 1] - self.time[i] - 1
            itime.append(self.time[i])
            irate.append(self.rate[i])
            ierror.append(self.error[i])
            if time_diff != 0.0:
                for j in range(1, int(time_diff) + 1):
                    itime.append(self.time[i] + float(j))
                    rate_variation = (self.rate[i + 1] - self.rate[i]) / (time_diff + 1)
                    error_variation = (self.error[i + 1] - self.error[i]) / (time_diff + 1)
                    irate.append(self.rate[i] + rate_variation * j)
                    ierror.append(self.error[i] + error_variation * j)
                new_idx = len(itime)
            else:
                new_idx += 1

        self.time = itime
        self.rate = irate
        self.error = ierror

    
    def smooth(self, smoothness, wma_width):
        """smooth self.rate with wma."""

        s = copy.copy(self.rate)
        for i in range(smoothness):
            s = self.wma_err(wma_width)
        self.smoothed = s

    
    def wma_err(self, n):
        """
        Weighted Moving Average ([n] points, use error breadth).
        [n] must be an odd number.
        """
        s = copy.copy(self.rate)
        n0 = int(n / 2)

        # both ends
        for i in range(1, n0 - 1):
            l_end = (self.rate[i] / self.error[i]) * float(n0 + 1)
            r_end = self.rate[-i] / self.error[-1] 
            l_weightsum = 1.0 / self.error[ i]
            r_weightsum = 1.0 / self.error[-1]
            
            for j in range(1, i + 1):
                l_end += ((self.rate[i - j] / self.error[i - j]) +
                          (self.rate[i + j] / self.error[i + j])) * float(n0 + 1 - j)
                r_end += ((self.rate[-(i - j + 1)] / self.error[-(i - j + 1)]) +
                          (self.rate[-(i + j + 1)] / self.error[-(i + j + 1)])) * float(n0 + 1 - j)
                l_weightsum += (1.0 / self.error[i - j]) + (1.0 / self.error[i + j]) * float(n0 + 1 - j)
                l_weightsum += (1.0 / self.error[-(i - j + 1)]) + (1.0 / self.error[-(i + j + 1)]) * float(n0 + 1 - j)
                
            s[ i] = l_end / l_weightsum
            s[-i] = r_end / r_weightsum

        # except ends
        for i in range(n0, len(self.rate) - n0):
            t = self.rate[i] / self.error[i]
            weightsum = (1.0 / self.error[i])
            for j in range(1, n0):
                t += (self.rate[i - j] / self.error[i - j]) + (self.rate[i + j] / self.error[i + j]) 
                weightsum += (1.0 / self.error[i - j]) + (1.0 / self.error[i + j])
            s[i] = t / weightsum

        return s
        
