import numpy as np
from matplotlib import pyplot as plt

def loop_scheduling():
    time1 = [31.793, 17.597, 11.708, 9.080, 7.344, 6.138, 5.296, 4.697, 4.4183, 3.810, 3.432] #static
    time2 = [32.025, 16.664, 11.507, 8.654, 6.972, 5.940, 5.154, 4.636, 4.523, 4.281, 3.406] #dynamic
    time3 = [30.330, 16.933, 11.561, 8.744, 7.124, 5.901, 5.091, 4.577, 4.118, 3.744, 3.470] #guided

    nb_cores = [1,2,3,4,5,6,7,8,9,10, 11]

    plt.plot(nb_cores, time1, 'b-', label='static')
    plt.plot(nb_cores, time1, 'bo')
    plt.plot(nb_cores, time2, 'r-', label='dynamic')
    plt.plot(nb_cores, time2, 'ro')
    plt.plot(nb_cores, time3, 'g-', label='guided')
    plt.plot(nb_cores, time3, 'go')
    plt.legend()
    plt.xlabel('Number of cores')
    plt.ylabel('Real Time (s)')
    plt.title('Execution time for different number of cores with different kind of loop')
    file_name = "dif_loop"
    plt.savefig("%s.pdf" %file_name)

def strong_scaling():
    time1 = [78.382, 41.497, 28.978, 21.502, 16.591, 13.866, 12.185, 10.822, 9.686, 8.661, 7.605] #avec O0
    time2 = [30.839, 16.989, 11.844, 9.083, 7.211, 6.058, 5.524, 5.112, 4.4182, 3.767, 3.426] #avec O2

    nb_cores= [1,2,3,4,5,6,7,8,9,10, 11]

    plt.plot(nb_cores, time1, 'b-', label='O0')
    plt.plot(nb_cores, time1, 'bo')
    plt.plot(nb_cores, time2, 'r-', label='O2')
    plt.plot(nb_cores, time2, 'ro')
    plt.legend()
    plt.xlabel('Number of cores')
    plt.ylabel('Real Time (s)')
    plt.title('Execution time for different number of cores')


    for x,y in zip(nb_cores,time1):

        label = "{:.2f}".format(y)

        plt.annotate(label,  # this is the text
                    (x,y), # this is the point to label
                    textcoords="offset points", # how to position the text
                    xytext=(0,4), # distance from text to points (x,y)
                    ha='left') # horizontal alignment can be left, right or center

    for x,y in zip(nb_cores,time2):

        label = "{:.2f}".format(y)

        plt.annotate(label,  # this is the text
                    (x,y), # this is the point to label
                    textcoords="offset points", # how to position the text
                    xytext=(0,4), # distance from text to points (x,y)
                    ha='left') # horizontal alignment can be left, right or center



    file_name = "strong_scaling"
    plt.savefig("%s.pdf" %file_name)

if __name__ == "__main__":
    strong_scaling()
    #loop_scheduling()