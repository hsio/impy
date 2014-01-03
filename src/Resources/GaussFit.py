import math
import numpy
import scipy.optimize
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import sys


def Gaussian(x, A, mu, sigma):
    """Evaluate a Gaussian distribution.

    :param x: the independent variable
    :param A: amplitude of the distribution
    :param mu: mean of the distribution
    :param sigma: standard deviation of the distribution
    :returns: distribution evaluated at x
    """
    # 2.50663 = sqrt(2*Pi) evaluated for speed in python (interpreted lang)
    # evaluated a sqrt 2 in exp to speed things up
    return A * math.exp(-math.pow( (x - mu)/(1.41421 * sigma), 2)) / (2.50663 * sigma)


class GaussFit(object):
    """Wrapper class for performing Gaussian fits to data.

    :param data: The data to fit. Must be 3 x n array (i.e. 3 columns, n row).
    :param guess: (optional) The initial guess for the fitting routine. Default [A,mu,sigma] = [5e7,10,1]
    :param restrict_chi2: (optional) Whether we should restrict chi2 calculations to +/- 5 sigma of the mean. Default is true.
    :param name: (optional) a text string describing this dataset / fit
    :author: Alex Zylstra
    :date: 2013/07/05
    """

    def __init__(self, data, guess=[5e7, 10, 1], restrict_chi2=True, name=""):
        """Constructor."""
        # super constructor:
        super(GaussFit, self).__init__()

        # initialization:
        self.data = []  # multi-dimensional array containing data
        self.data_x = []  # x values of the data
        self.data_y = []  # y values of the data
        self.data_err = []  # error bars for the data
        self.fit = []  # result of the fit
        self.covariance = []  # covariance results of the fit

        # Generic guesses for fitting:
        self.GUESS = []

        # if the chi2 calculations should restrict limits to 5sigma of the mean
        self.OPT_RESTRICT_CHI2 = True

        # name of the dataset / fit
        name = ""

        # copy data:
        self.data = numpy.copy(data)

        # set guess parameters:
        self.GUESS = guess

        # set name:
        self.name = name

        # set restrict chi2 option:
        self.OPT_RESTRICT_CHI2 = restrict_chi2

        # split data into components:
        self.data_x = numpy.ndarray(len(data))
        self.data_y = numpy.ndarray(len(data))
        self.data_err = numpy.ndarray(len(data))
        for i in range(len(data)):
            self.data_x[i] = data[i][0]
            self.data_y[i] = data[i][1]
            self.data_err[i] = data[i][2]

        # perform the fit:
        self.fit, self.covariance = self.do_fit()

    def do_fit(self, method='chi^2', guess=[], fixed=[]) -> tuple:
        """Perform the fit. Called automatically by the constructor.
        Can be invoked again with different initial guess if desired.

        :param method: (optional) Which method to use. Default is chi^2. Options are:

            'fmin' downhill simplex minimization: minimizes chi^2 for the dataset using `scipy.optimize.fmin`

            'leastsq' least-squares fit (warning: does not use error bars). Uses `scipy.optimize.leastsq`

            'chi^2' chi^2 minimization using the `scipy.optimize.curve_fit` algorithm

        :param guess: (optional) The initial guess for the fitting routine.
        :returns: tuple containing best fit parameters, and covariance matrix
        """
        if len(guess) != 3:
            guess = self.GUESS
        if len(fixed) != len(guess):
            fixed = []
            for i in guess:
                fixed.append(False)

        if method == 'fmin':
            # Do fit using downhill simplex:
            fit = scipy.optimize.fmin(
                func=self.chi2_other, # function to minimize
                x0=guess, # initial guess
                disp=False) # turn off console output

            return fit, []

        # least squares minimization in scipy
        if method == 'leastsq':
            def gaussian(B,x):
                return B[0]/(B[2]*numpy.sqrt(2*numpy.pi))*numpy.exp(-((x-B[1])**2/(2*B[2]**2)))
            def func(p, x, y):
                return y-gaussian(p, x)

            result = scipy.optimize.leastsq(func, guess, args=(self.data_x, self.data_y), full_output=True)

            if result[4] != (1 or 2 or 3 or 4):
                print("Error fitting data\n",result[3])
            return result[0], result[1]

        # curve fit routine in scipy
        else:
            def gaussian2(x, A, mu, sigma):
                return A/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-((x-mu)**2/(2*sigma**2)))

            result = scipy.optimize.curve_fit(gaussian2, self.data_x, self.data_y, p0=guess, sigma=self.data_err)
            return result[0], result[1]

    def get_fit(self) -> list:
        """Get the calculated best fit.

        :returns: a python list containing [a,mu,sigma] from the best fit.
        """
        return self.fit

    def get_covariance(self):
        """Get covariance matrix for the best fit.

        :returns: covariance matrix (see scipy.optimize.curve_fit documentation)
        """
        return self.covariance

    def eval_fit(self, x):
        """Evaluate the best fit.

        :param x: the independent variable
        :returns: the fit evaluated at x
        """
        return Gaussian(x, self.fit[0], self.fit[1], self.fit[2])

    def chi2(self):
        """Calculate total chi2 for best fit.

        :returns: the value of chi2
        """
        return self.chi2_other(self.fit)

    def chi2_other(self, fit):
        """Calculate chi2 for any fit to this data.

        :param fit: a list containing the fit parameters [A,mu,sigma]
        :returns: chi2 for class data and fit
        """
        chi2 = 0

        # iterate over all data:
        for point in self.data:
            # sanity check that error bar is not zero:
            if point[2] != 0:
                # if we want to restrict chi2, only count if we are within 5 sigma of mean:
                if ( self.OPT_RESTRICT_CHI2 and math.fabs(point[0] - fit[1]) / fit[2] <= 5 ):
                    chi2 += ( Gaussian(point[0], fit[0], fit[1], fit[2]) - point[1] ) ** 2 / (point[2]) ** 2
                # otherwise, count all points:
                else:
                    chi2 += ( Gaussian(point[0], fit[0], fit[1], fit[2]) - point[1] ) ** 2 / (point[2]) ** 2

        return chi2

    def red_chi2(self):
        """Calculate reduced chi2 for the best fit.

        :returns: value of reduced chi2 for the best fit
        """
        return self.red_chi2_other(self.fit)

    def red_chi2_other(self, fit):
        """Calculate reduced chi2 for any fit to this data.

        :param fit: an array containing the fit parameters [A,mu,sigma]
        :returns: reduced chi2 for class data and fit
        """
        return self.chi2_other(fit) / ( len(self.data) - 3 )

    def delta_chi2_amp(self, dA, sign=1):
        """Calculate increase in chi2 due to a change in amplitude.

        :param dA: dA change in amplitude (A' = A + dA)
        :param sign: (optional) whether amplitude should be increased or decreased (default=1 -> increase)
        :returns: the change in chi2
        """
        A = self.fit[0] + sign * dA
        # calculate a new best fit:
        def temp_gauss(x, mu, sigma):
            return A/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-((x-mu)**2/(2*sigma**2)))
        guess = [self.fit[1], self.fit[2]]
        result = scipy.optimize.curve_fit(temp_gauss, self.data_x, self.data_y, p0=guess, sigma=self.data_err)
        fit = result[0]

        chi2 = self.chi2_other([A, fit[0], fit[1]])

        return chi2 - self.chi2()

    def delta_chi2_mu(self, dmu, sign=1):
        """Calculate increase in chi2 due to a change in mean.

        :param dmu: change in mean (mu' = mu + dmu)
        :param sign: (optional) whether mean should be increased or decreased (default=1 -> increase)
        :returns: change in chi2
        """
        mu = self.fit[1] + sign * dmu
        # calculate a new best fit:
        def temp_gauss(x, A, sigma):
            return A/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-((x-mu)**2/(2*sigma**2)))
        guess = [self.fit[0], self.fit[2]]
        result = scipy.optimize.curve_fit(temp_gauss, self.data_x, self.data_y, p0=guess, sigma=self.data_err)
        fit = result[0]

        chi2 = self.chi2_other([fit[0], mu, fit[1]])

        return chi2 - self.chi2()

    def delta_chi2_sigma(self, ds, sign=1):
        """Calculate increase in chi2 due to a change in sigma.

        :param ds: change in sigma (sigma' = sigma + ds)
        :param sign: (optional) whether sigma should be increased or decreased (default=1 -> increase)
        :returns: change in chi2
        """
        sigma = self.fit[2] + sign * ds
        # calculate a new best fit:
        def temp_gauss(x, A, mu):
            return A/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-((x-mu)**2/(2*sigma**2)))
        guess = [self.fit[0], self.fit[1]]
        result = scipy.optimize.curve_fit(temp_gauss, self.data_x, self.data_y, p0=guess, sigma=self.data_err)
        fit = result[0]

        chi2 = self.chi2_other([fit[0], fit[1], sigma])

        return chi2 - self.chi2()

    def chi2_fit_unc(self):
        """Calculate uncertainty in the fit parameters.
        Routine: each fit parameter is varied to produce an increase of 1 in total chi2

        :returns: list containing uncertainties [ [-A,+A] , [-mu,+mu] , [-sigma,+sigma] ]
        """
        # return value: delta unc in each parameter
        delta = [[], [], []]

        # Calculate uncertainties
        # need to do for both + and -
        for sign in [-1, 1]:
            # calculate uncertainty in the amplitude:
            dA = scipy.optimize.fminbound(
                func=(lambda x: math.fabs(self.delta_chi2_amp(x, sign) - 1)), # min @ delta chi2 = 1
                x1=0, x2=2*math.fabs(self.fit[0]), # min is 0, max is 2x nominal
                full_output=False) # suppress full output

            # calculate uncertainty in the mean:
            dmu = scipy.optimize.fminbound(
                func=(lambda x: math.fabs(self.delta_chi2_mu(x, sign) - 1)), # min @ delta chi2 = 1
                x1=0, x2=self.fit[2], # min/max change +/- sigma
                full_output=False) # suppress full output

            # calculate uncertainty in sigma:
            ds = scipy.optimize.fminbound(
                func=(lambda x: math.fabs(self.delta_chi2_sigma(x, sign) - 1)), # min @ delta chi2 = 1
                x1=0, x2=self.fit[2] * 0.99, # min/max change +/- sigma
                full_output=False) # suppress full output

            # add info to delta:
            delta[0].append(sign * dA)
            delta[1].append(sign * dmu)
            delta[2].append(sign * ds)

        return delta

    def plot_file(self, fname):
        """Save a plot to file.

        :param fname: the file to save
        """
        # import matplotlib
        import matplotlib
        import matplotlib.pyplot as plt
        if matplotlib.get_backend() != 'agg':
            plt.switch_backend('Agg')
        # get the figure:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.plot(ax)

        # save to file:
        fig.savefig(fname)

    def plot_window(self, interactive=False):
        """Make a plot in a new UI window

        :param interactive: (optional) Whether to show the plot in interactive mode {default = false}
        """
        # import matplotlib
        import matplotlib
        import matplotlib.pyplot as plt

        # os detection
        if sys.platform.startswith('linux'):  # generic *nix
            plt.switch_backend('TkAgg')
        elif sys.platform.startswith('darwin'):  # Mac OS X
            if matplotlib.get_backend() != 'MacOSX':
                plt.switch_backend('MacOSX')
        # use interactive mode if requested:
        plt.interactive(interactive)

        # get the figure:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.plot(ax)

        plt.show()

    def plot(self, ax=None):
        """Make a plot of the data and fit, drawn into given Axes.

        :param ax: Axes instance to plot into
        """
        # sanity check:
        if ax is None:
            ax = plt.gca()

        # plot the data with error bars
        ax.errorbar(
            self.data_x, # x values
            self.data_y, # y values
            yerr=self.data_err, # y error bars
            marker='s', # square markers
            lw=0, # no lines connecting points
            elinewidth=1, # error bar line width
            mfc='black', # marker color
            mec='black', # marker line color
            ecolor='black') # error bar color

        # make an evaluted fit dataset for plotting
        fit_x = []
        fit_y = []
        x = self.data_x[0]
        dx = (self.data_x[1] - x) / 10
        while x < self.data_x[-1]:
            fit_x.append(x)
            fit_y.append(self.eval_fit(x))
            x += dx

        ax.plot(fit_x, fit_y, 'r--')

        # add text with fit parameters
        # location for the text: (based on line position)
        x = self.fit[1] + 4 * self.fit[2]
        y = self.fit[0] * 6 / 8
        # construct a text string to display:
        text = r'$Y_p$ = ' + '{:.2e}'.format(self.fit[0])
        ax.text(x, y, # data
                text, # text to display
                backgroundcolor='white') # fill background
        # new line of text
        y = self.fit[0] * 5 / 8
        text = r'$E_p$ = ' + '{:.2f}'.format(self.fit[1]) + ' MeV'
        ax.text(x, y, # data
                text, # text to display
                backgroundcolor='white') # fill background
        # new line of text
        y = self.fit[0] * 4 / 8
        text = r'$\sigma_p$ = ' + '{:.2f}'.format(self.fit[2]) + ' MeV'
        ax.text(x, y, # data
                text, # text to display
                backgroundcolor='white') # fill background
        # new line of text
        y = self.fit[0] * 3 / 8
        text = r'$\chi^2$ red. = ' + '{:.2f}'.format(self.red_chi2())
        ax.text(x, y, # data
                text, # text to display
                backgroundcolor='white') # fill background

        # write text to the plot:
        # ax.text(x, y, # data
        #         text, # text to display
        #         backgroundcolor='white', # fill background
                #bbox=dict(fc='white', ec='black', alpha=1.0)) # add black boundary

        ax.grid(True)
        ax.set_xlabel('Energy (MeV)')
        ax.set_ylabel('Yield / MeV')
        ax.set_title(self.name + ' Fit')

        #return fig
