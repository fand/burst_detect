

class Window:

    def __init__(self, head, tail, median, threshold, burst_h, burst_l):
        self.head = head  # head index of window in seq (Note that this is index, not date)
        self.tail = tail
        self.median = median
        self.threshold = threshold
        self.burst_h = burst_h
        self.burst_l = burst_l
