# TWED: Time Warp Edit Distance

Te Time Warp Edit Distance (TWED) is a distance measure for discrete time series matching with time 'elasticity'. In comparison to other distance measures, (e.g. DTW (Dynamic Time Warping) or LCS (Longest Common Subsequence Problem)), TWED is a metric. Its computational time complexity is O(n^2), but can be drastically reduced in some specific situation by using a corridor to reduce the search space. Its memory space complexity can be reduced to O(n). It was first proposed in [1].

[1] Marteau, P.; F. (2009). "Time Warp Edit Distance with Stiffness Adjustment for Time Series Matching". IEEE Transactions on Pattern Analysis and Machine Intelligence. 31 (2): 306–318.
