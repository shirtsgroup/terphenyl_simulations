import os
import sys
import glob
import pymbar
import alchemlyb
import time as timer
import subprocess
import argparse
import natsort
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm


def initialize(args):
    parser = argparse.ArgumentParser(
        description='This code analyzes the log file generated by replica \
                    exchange molecular dynamics (REMD) simulations and plot \
                    the trajectory of each replica as a function of time.')
    parser.add_argument('-l',
                        '--log',
                        type=str,
                        nargs='+',
                        help='The file name of the log file.')
    parser.add_argument('-p',
                        '--prefix',
                        type=str,
                        help='The common prefix of the simulation files.')
    parser.add_argument('-d',
                        '--diag',
                        default=False,
                        action='store_true',
                        help='Whethter to plot the diagonal of the transition \
                            matrix as a function of time. "-d" specified meands \
                            that the plot will be generated.')

    args_parse = parser.parse_args(args)

    return args_parse


class LogInfo:
    def __init__(self, logfiles):
        """
        Gets the needed parameters and data from the log file and set up
        relevant attributes to run the analysis

        Parameters
        ----------
        logfiles : list
            A list of the filenames of log files
        """
        # the info from all th log files should be the same
        f = open(logfiles[0], 'r')
        lines = f.readlines()
        f.close()

        self.nex = None
        self.replex = None
        self.N_states = None
        self.dt = None

        line_n = 0
        for l in lines:
            line_n += 1

            if 'dt  ' in l and self.dt is None:
                self.dt = float(l.split('=')[1])

            if 'Command line' in l:
                if 'nex' in lines[line_n]:
                    self.nex = True  
                else:
                    self.nex = False            

                if 'replex' in lines[line_n]:
                    self.replex = float(lines[line_n].split('replex')[1].split()[0])

            if 'Replica exchange in ' in l:
                self.N_states = len(lines[line_n].split())

            if 'Started mdrun' in l:
                self.start = line_n
                # the line number that the simulation got started
                break


class REMDAnalysis(LogInfo):
    """
    A class for state-time analysis of replica exchange molecular dynamics simulations. When
    instantiating this class, one parameter (logfiles) is required.
    """

    def __init__(self, logfiles):
        """
        Sets up the properties of the instance of REMDAnalysis
        """
        self.sample_all = None
        self.finish = False   # if the simulation finishes all the steps specified
        self.n_ex = 0        # number of exchanges
        self.final_t = None
        LogInfo.__init__(self, logfiles)

    def get_replica_data(self, logfiles, calc_diag=False):
        """
        Gets the data of of each replica from the log file, including the state-time 
        data and transition matrix. Note that according to the log file in GROMACS, 
        index (i, j) in the transition matrix represents the number of exchanges 
        "BETWEEN SPOT i AND j" divided the total number of exchanges.

        Parameters
        ----------
        logfiles : list
            A list of the names of log files
        calc_diag: bool
            Whethter to calculate the diagonal of the transition matrix as a function 
            of time.

        Returns
        -------
        transition_matrix : np.array
            The transition matrix.
        """

        for k in range(len(logfiles)):
            f = open(logfiles[k], 'r')
            lines = f.readlines()
            f.close()

            # Part 1: Get state-time data
            line_n = 0
            
            if k == 0:
                state_data = []
                if calc_diag is True:
                    diag_data = []

                for i in range(self.N_states):
                    # data_all[i] represents the states that replica i has visited
                    state_data.append([])
                    state_data[i].append(i) # replica i starts from state i  
                    if calc_diag is True:
                        diag_data.append([])

            for l in lines[self.start:]:
                line_n += 1
                if 'Order After Exchange' in l:
                    self.n_ex += 1   # this will accumulate across iterations
                    state_list = [int(i) for i in l.split(':')[1].split()]
                    for i in range(self.N_states):
                        state_data[i].append(state_list.index(i))
        
        # for state-time plotting
        self.final_t = self.n_ex * self.replex * self.dt  # units: ps
        time = np.linspace(0, self.final_t, self.n_ex + 1)

        # Part 2: Get the data of transition matrix
        transition_matrix = np.zeros([self.N_states, self.N_states])

        # if calc_diag is True, then we have to go over the whole file to calculate 
        # transition matrix anyway, so it doesn't matter if self.finish is True or not
        if calc_diag is False:
            try:
                # Use external shell command to search the string would be faster
                subprocess.check_output("grep 'Replica exchange statistics' %s" % logfiles[-1], shell=True)
                self.finish = True
            except:
                # Handle the error which would occur if the string is not found
                self.finish = False

            if self.finish is True:
                # Then read from the bottom to get the transition matrix
                # (since if the simulation is finished, the overlap matrix will be calculated automatically)
                # lines here are from the last logfile
                lines.reverse()    # from this point, lines have been reversed
                line_n = 0
                for l in lines:
                    # print(l)    # this will print from the bottom
                    line_n += 1
                    if 'Replica exchange statistics' in l:
                        # data starts from lines[line_n - 5]
                        for i in range(self.N_states):
                            transition_matrix[i] = [float(v) for v in lines[line_n - 5 -
                                                                            i].split('Repl')[1].split()[:self.N_states]]

        if self.finish is False or calc_diag is True:
            # Then we have to calculate transition matrix on our own as follows.
            count_matrix = np.zeros([self.N_states, self.N_states])

            for k in range(len(logfiles)):
                f = open(logfiles[k], 'r')
                lines = f.readlines()
                f.close()

                line_n = 0
                for l in lines[self.start:]:  # skip the metadata
                    line_n += 1
                    if 'Accepted Exchanges' in l:
                        rep_list = [int(i) for i in l.split(':')[1].split()]
                        for i in rep_list:
                            # Note that for the diagnal: count_matrix[i, i] += 2
                            count_matrix[i, rep_list.index(i)] += 1
                            count_matrix[rep_list.index(i), i] += 1

                        if calc_diag is True:
                            # diag_prob = the diagonal probability of each state at a certain time frame
                            diag_prob = [count_matrix[i, i] / sum(count_matrix[i]) for i in range(self.N_states)]
                            
                            for i in range(len(diag_prob)):
                                diag_data[i].append(diag_prob[i])

            # get transition_matrix from the data of count_matrix
            for i in range(self.N_states):
                for j in range(self.N_states):
                    transition_matrix[i, j] = count_matrix[i, j] / sum(count_matrix[i])
        
        if calc_diag is False:
            return time, state_data, transition_matrix
        elif calc_diag is True:
            return time, state_data, np.array(diag_data), transition_matrix

    def plot_replica_data(self, time, data, png_name, diag=None, plot_type=None, n_subplots=None, start_idx=1):
        if n_subplots is None:
            n_subplots = self.N_states

        if int(np.sqrt(n_subplots) + 0.5) ** 2 == n_subplots:
            # perfect sqaure number
            n_cols = int(np.sqrt(n_subplots))
        else:
            n_cols = int(np.floor(np.sqrt(n_subplots))) + 1

        if n_subplots % n_cols == 0:
            n_rows = int(n_subplots / n_cols)
        else:
            n_rows = int(np.floor(n_subplots / n_cols)) + 1

        _, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(2.5 * n_cols, 2.5 * n_rows))
        if plot_type == 'state':
            plt.suptitle('Exploration of states as a function of time', weight='bold', fontsize=14)
        if plot_type == 'diag':
            plt.suptitle('Probability of staying at the original state as a function of time', weight='bold', fontsize=14)

        for i in range(n_subplots):
            plt.subplot(n_rows, n_cols, i + 1)
            plt.plot(np.array(time) / 1000, data[i])
            
            if diag is not None:
                # r = probability of staying at the original state
                plt.annotate('(Replica %s, r = %s%%)' % (i + start_idx, diag[i]), xy=(0, 0), xytext=(0, self.N_states * 0.05))

            if (i + 1) % n_cols == 1:
                if plot_type == 'state':
                    plt.ylabel('State')
                elif plot_type == 'diag':
                    plt.ylabel('Diagonal probability (%)')
                    
            if (i + 1) > (n_rows - 1) * n_cols:
                    plt.xlabel('Time (ns)')
            
            if plot_type == 'state':
                plt.ylim([0, self.N_states - 1])
            elif plot_type == 'diag':
                plt.ylim([0, 100])
            plt.grid()

        # Remove redundant subplots
        n_rm = n_cols * n_rows - n_subplots
        for i in range(n_rm):
            ax.flat[-1 * (i + 1)].set_visible(False)

        plt.tight_layout(pad=5.0, w_pad=0.5, h_pad=2.0)
        plt.savefig(png_name, dpi=600)
        # plt.show()
        plt.close()
    
    def plot_matrix(self, matrix, png_name, start_idx=0):
        sns.set_context(rc={
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 5
        })

        K = len(matrix)
        plt.figure(figsize=(K / 3, K / 3))
        annot_matrix = np.zeros([K, K])   # matrix for annotating values

        mask = []
        for i in range(K):
            mask.append([])
            for j in range(len(matrix[0])):
                if matrix[i][j] < 0.005:            
                    mask[-1].append(True)
                else:
                    mask[-1].append(False)

        for i in range(K):
            for j in range(K):
                annot_matrix[i, j] = round(matrix[i, j], 2)

        x_tick_labels = y_tick_labels = np.arange(start_idx, start_idx + K)
        ax = sns.heatmap(matrix, cmap="YlGnBu", linecolor='silver', linewidth=0.25,
                        annot=annot_matrix, square=True, mask=matrix < 0.005, fmt='.2f', cbar=False, xticklabels=x_tick_labels, yticklabels=y_tick_labels)
        ax.xaxis.tick_top()
        ax.tick_params(length=0)
        cmap = cm.get_cmap('YlGnBu')   # to get the facecolor
        ax.set_facecolor(cmap(0))      # use the brightest color (value = 0)
        for _, spine in ax.spines.items():
            spine.set_visible(True)    # add frames to the heat map
        plt.annotate('$\lambda$', xy=(0, 0), xytext=(-0.45, -0.20))
        plt.title('Transition matrix', fontsize=14, weight='bold')
        plt.tight_layout(pad=1.0)

        plt.savefig(png_name, dpi=600)
        # plt.show()
        plt.close()

class MBARAnalysis(REMDAnalysis):
    """
    A class using MBAR to perform free energy calculations.
    """

    def __init__(self, logfile):
        """
        Sets up the properties of the instance of MBARAnalysis
        """
        LogInfo.__init__(self, logfile)
        REMDAnalysis.__init__(self, logfile)

    def decorrelate_data(self, start, end):
        K = self.N_states
        u_kln = np.zeros([K, K, max(end-start), np.float64])
        N_k = np.zeros(K, int)  # the number of uncorrelated samples from state k
        g = np.zeros(K, float)  # correlation times for the data
        

    def get_overlap_matrix():
        # this matrix can be obtained only if MBAR is used
        pass

        #f = open(log, 'r')
        #lines = f.readlines()
        # f.close()


def main():
    start = timer.time()

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    args = initialize(sys.argv[1:])

    if args.log is None:
        args.log = []
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log.append(file)
        if not args.log:
            print('No log files provided or found! Please check if the dirctory is correct or specify the name of the log file.')
            sys.exit()
        else:
            args.log = natsort.natsorted(args.log, reverse=False)

    if args.prefix is None:
        args.prefix = args.log[0].split('.')[0]

    log_str = ''
    for i in range(len(args.log)):
        log_str += args.log[i]
        if i == len(args.log) - 2:
            log_str += ', and '
        elif i == len(args.log) - 1:
            pass
        else:
            log_str += ', '

    result_str = '\nData analysis of the file(s): %s' % log_str
    print(result_str)
    print('=' * (len(result_str) - 1))  # len(result_str) includes \n 
    
    print('Analyzing the log file ...')
    RA = REMDAnalysis(args.log)

    if args.diag is False:
        time, state, t_matrix = RA.get_replica_data(args.log, calc_diag=args.diag)
    elif args.diag is True:
        time, state, diag_prob, t_matrix = RA.get_replica_data(args.log, calc_diag=args.diag)

    # the probability of staying at the original state (at the last time frame)
    diag = np.round(100 * np.diagonal(t_matrix), 1) 

    print('Simulation length: %s ns (%s exchanges occured.)\n' % (RA.final_t / 1000, RA.n_ex))
    
    # Some parameters to use for larger numbers of states
    if RA.N_states > 80: # using 49 as a unit, separate the plot into several
        n_figs = int(np.ceil(RA.N_states / 49))
        n_subplots = int(np.floor(RA.N_states / n_figs))  # min number of subplots
        n_list = np.zeros(n_figs)  # a list of number of subplots of each figure

        # Consider the following equation
        # x: number of figures that contain n_subplots subplots
        # y: number of figures that contain (n_subplots + 1) subplots
        # nx + (n+1)y = self.N_states (n: n_subplots)
        # x + y = n_figs
        # For example, if n_figs = 3, self.N_states = 83, then n_list = np.array(27, 28, 28])
        # and bounds = [0, 27, 55, 83]

        y = RA.N_states - n_subplots * n_figs
        x = n_figs - y

        for i in range(x):
            n_list[i] = n_subplots
        for i in range(y):
            n_list[i + x] = n_subplots + 1
        
        bounds = [0]
        for i in range(n_figs):
            bounds.append(int(np.sum(n_list[:i + 1])))

    print('Plotting the exploration of states as a function of time ...')
    if RA.N_states <= 80:
        RA.plot_replica_data(time, state, 'state_time_%s.png' % args.prefix, diag=diag, plot_type='state')
        print('The state time plot, state_time_%s.png, has been generated.\n % args.prefix')
    else:   
        for i in range(n_figs):
            RA.plot_replica_data(time, state[bounds[i]:bounds[i + 1]], 'state_time_%s_part%s.png' % (args.prefix, i + 1), diag=diag, plot_type='state', n_subplots=int(n_list[i]), start_idx=bounds[i])
            print('The state time plot, state_time_%s_part%s.png, has been generated.\n' % (args.prefix, i + 1))
    
    if args.diag is True:
        print('Plotting the diagonal of the transition matrix as a function of time ...')
        if RA.N_states <= 80:
            RA.plot_replica_data(time[1:], diag_prob * 100, 'diagprob_time_%s.png' % args.prefix, diag=diag, plot_type='diag')
            print('The plot of diagonal proability, diag_prob_time_%s.png, has been generated.\n' % args.prefix)
        else:
            for i in range(n_figs):
                RA.plot_replica_data(time[1:], diag_prob[bounds[i]:bounds[i + 1]] * 100, 'diagprob_time_%s_part%s.png' % (args.prefix, i + 1), diag=diag, plot_type='diag', n_subplots=int(n_list[i]), start_idx=bounds[i])
                print('The plot of diagonal probability, diag_prob_time_%s_part%s.png, has been generated.' % (args.prefix, i + 1))
    
    print('Plotting the transition matrix as a heat map ...')
    if RA.N_states <= 60:
        RA.plot_matrix(t_matrix, 'transition_matrix_%s.png' % args.prefix)
        print('The heat map of the transition matrix, transition_matrix_%s.png, has been generated.\n' % args.prefix)
    else:
        # Note that generally a transition matrix for a REMD simulation has a lot of zero elements
        # so for a REMD with a larger number of states, we only plot the elements at the diagonal +- 7
        # Therefore, n_fig * x - (n_fig - 1) * 7 = self.N_states  (x: n_replica, # of replicas to plot)
        n_figs = int(np.ceil(RA.N_states / 36))
        n_list = np.zeros(n_figs)
        n_replicas = int(np.floor((RA.N_states + 7 * (n_figs - 1)) / n_figs))  # analogous to n_subplots when using plot_replica_data
        y = (RA.N_states + 7 * (n_figs - 1)) - n_replicas * n_figs
        x = n_figs - y

        for i in range(x):
            n_list[i] = n_replicas
        for i in range(x, x + y):
            n_list[i] = n_replicas + 1
        
        low_bounds = [0]
        # For example, for 126 states, n_list = np.array([46, 47, 47])
        # nn_list = np.array([39, 40, 40])
        # so low_bounds = [0, 39, 79], high_bounds = [46, 86, 126]
        nn_list = n_list - 7
        for i in range(n_figs - 1):
            low_bounds.append(int(np.sum(nn_list[:i + 1])))
        high_bounds = low_bounds + n_list

        subdata = []
        for i in range(n_figs):
            subdata.append(np.zeros([int(n_list[i]), int(n_list[i])]))
            for j in range(int(n_list[i])):
                subdata[i][j] = t_matrix[j + low_bounds[i]][low_bounds[i]:int(high_bounds[i])]
            RA.plot_matrix(subdata[i], 'transition_matrix_%s_part%s.png' % (args.prefix, i + 1), start_idx=low_bounds[i])
            print('The heat map of the transition matrix, transition_matrix_%s_part%s.png, has been generated.\n' % (args.prefix, i + 1))


    print('Plotting the histogram of the diagonal transition probability ...')
    plt.figure()
    plt.hist(diag, bins=20, edgecolor='black')
    plt.xlabel('Probability of staying at the original state (%)')
    plt.ylabel('Count of states')
    plt.title('Histogram of the diagonal transition probaility')
    plt.grid()
    plt.savefig('hist_diag_%s.png' % args.prefix, dpi=600)
    # plt.show()
    plt.close()
    print('The histogram of the diagonal transition probability, hist_diag_%s.png, has been generated.\n' % args.prefix)

    end = timer.time()
    print('Total time elasped: %s seconds.\n' % (end - start))