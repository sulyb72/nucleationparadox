import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns
import scipy.optimize

def histoPlot(temp, maxPlot, τ0):
    e=math.e
    
    file = temp+"K-waitingtimes.txt"
    waitingTimes = np.loadtxt(file)
    
    binsNum=100
    rangeMin=0.0
    rangeMax=1.0
    
    #sns.displot(waitingTimes, bins=binsNum, kde=1, color="firebrick").set(title="Histogram of Waiting Times for "+temp+"K Data")
    
    #filteredWTs = waitingTimes[(waitingTimes >= 0) & (waitingTimes < 1)]
    #sns.displot(filteredWTs, bins=binsNum, kde=1, color="royalblue").set(title="Histogram of Waiting Times for "+temp+"K Data")
    
    #plt.show()
    
    plt.hist(waitingTimes, bins=binsNum, range=(rangeMin,rangeMax), color="royalblue", label="Histogram of Waiting Time Data")
    y, x = np.histogram(waitingTimes, bins=binsNum, range=(rangeMin,rangeMax))
    x_adjust=x[:-1]
    #plt.plot(x_adjust, y, color="royalblue"")
    
    def exponential(x, m, t):
        y_vals = m * np.exp(-t * x)
        return y_vals
    
    initialParams = ((1/τ0), (1/τ0))
    params, cv = scipy.optimize.curve_fit(exponential, x_adjust, y, initialParams)
    M, T = params
    y_curve = exponential(x_adjust, M, T)  
    
    initial_y = exponential(x_adjust, (1/τ0), (1/τ0))
    plt.plot(x_adjust, initial_y, color="darkorange", label="Predicted Exponential Waiting Time Distribution")
    
    plt.plot(x_adjust, y_curve, color="firebrick", label="Fitted Exponential Waiting Time Distribution")
    
    params_error = np.sqrt(np.diag(cv))
    M_err, T_err = params_error
    
    #sqr_misfits = np.square(y_curve-initial_y)
    #rms_misfit = (np.average(sqr_misfits))**0.5
    #print(rms_misfit)
    
    #sqr_diffs = np.square(y_curve - y)
    #sqr_diffs_from_mean = np.square(y - np.mean(y))
    #rSquared = 1 - np.sum(sqr_diffs) / np.sum(sqr_diffs_from_mean)
    #print("R² Coefficient of Determination: "+str(rSquared))
    
    print("-------------------------")
    print("Parameters of Exponential Fit to "+temp+"K data [m, t]: "+str(M)+" ± "+str(M_err)+", "+str(T)+" ± "+str(T_err))
    
    approxMean = (np.sum(x_adjust*y_curve))/(np.sum(y_curve))
    print("Approximate Mean of Exponential Fit: "+str(approxMean))
    
    plt.title("Histogram of Waiting Times for "+temp+"K Data")
    plt.xlabel("Waiting Time to Observe a Freezing Event, τv (ns)")
    plt.ylabel("Counts")
    #plt.xlim(rangeMin,rangeMax)
    plt.xlim(0, maxPlot)
    plt.ylim(0, max(y))
    plt.legend()
    
    filename="histogram_"+temp+".png"
    #plt.savefig(filename,dpi=300)
    
    plt.show()
    
histoPlot("3800", 0.2, 0.02253918)
histoPlot("4000", 1, 0.37497022)
histoPlot("4200", 0.2, 0.07094806)