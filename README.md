# SILICON SIGNALS RISE TIME
# Advance Scientific Python Course Project

The project consists of reading a bunch of binary files that contain millions of detectors signals from an experiment and studying their rise time as a function of the energy and mass. The signals are orginized by mass and energy on individual binary files. 

The idea is to make a program that finds the rise time of the signals and then creates a histogram with the founded values. The rise time of a signal is the number of time samples that there are from its baseline until its saturation. After the histograms are produced, a gaussian fit is performed and then the mean of the fits are stored on an object which is used to plot the evolution of the rise time as a function of the energy. 

<!-- After this the program needs to list a bunch of ions with a similar energy (E - 1 < E < E + 1), with this information and the mean value of the rise time for a given energy a plot of the rise time as a function of the mass is going to be performed. -->


```
	                  
	                    SILICON DETECTOR SIGNAL
         ______________________________________________________
         |                                                     | 
         |                                *********************|
         |                               *|                    |
         |                              * |                    |
    V    |                             *  |                    |
    O    |                            *   |                    |
    L    |                           *    |                    |
    T    |                          *     |                    |
    A    |                         *      |                    |
    G    |                        *       |                    |
    E    |                       *        |                    |
         |                      *         |                    |
   [V]   |                     *          |                    |
         |*********************           |                    |
         |                     |----------|                    |
         |                      RAISE TIME                     |
         |_____________________________________________________|
                                                    TIME [ns]           

```

# Ana Gomez






