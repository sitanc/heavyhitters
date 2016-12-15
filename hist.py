from collections import defaultdict
import operator
from itertools import izip
import os.path
import pickle
import numpy as np
from pqdict import pqdict
import subprocess
import matplotlib.pyplot as plt

# with open('temp.pickle', 'rb') as handle:
#     adj_list = pickle.load(handle)

def data_path(time, fname, trial_number, win_dep=True):
    if win_dep:
        return 'data/'+str(time)+'/trial_'+str(trial_number)+'/'+fname
    else:
        return 'data/'+str(time)+'/'+fname

def dump(time, fname,obj, trial_number):
    path = data_path(time, fname, trial_number)
    with open(path, 'wb') as handle:
        pickle.dump(obj, handle)

def load(time, fname, trial_number):
    path = data_path(time, fname, trial_number)
    with open(path, 'rb') as handle:
        return pickle.load(handle)

# compute minimum perc% vertex cover (greedily)
def percentile(adj_list,perc):
    # pqdict.popitem pops things of the lowest priority, hence the negative sign and addition below
    degrees = pqdict({k: -sum(adj_list[k].values()) for k in adj_list.keys()})
    total_degrees = -sum(degrees.values())/2
    cover = {}
    mass_covered = 0
    while True:
        best, degree = degrees.popitem()
        degree = -degree
        cover[best] = True
        for neighbor, weight in adj_list[best].items():
            try:
                adj_list[neighbor].pop(best)
                degrees[neighbor] += adj_list[best][neighbor]
            except:
                print "asdf"
        adj_list.pop(best)
        try:
            mass_covered += degree
        except:
            print mass_covered, degree
        print "MASS COVERED: %d" % mass_covered
        if mass_covered > perc * total_degrees:
            last_degree = degree
            break
    return cover, total_degrees, last_degree

# regard adjacency list just as a (double cover) of list of pairs
# literally just compute the 95th percentile over union of elements in pairs
def naive_percentile(adj_list,perc):
    # dictionary where keys are nodes and values are their weighted degrees, i.e. frequencies
    hist = {k: sum(adj_list[k].values()) for k in adj_list.keys()}
    mass = sum(hist.values())
    hist = sorted(hist.items(), key=operator.itemgetter(1))
    hist.reverse()
    so_far = 0
    ids = {}
    for i in range(len(hist)):
        so_far += hist[i][1]
        ids[hist[i][0]] = True
        if so_far >= perc*mass:
            last_degree = hist[i][1]
            break
    return ids, mass, last_degree

def add_to(l,x,y):
    # there are weird cases where an IP sends to itself, and we'll ignore these cases
    if x != y:
        if y not in l[x].keys():
            l[x][y] = 1
        else:
            l[x][y] += 1
        if x not in l[y].keys():
            l[y][x] = 1
        else:
            l[y][x] += 1

def local_get_heavys(window_size,perc,time,naive,saving=False):
    if saving:
        heavy_fname = 'heavys_'+str(window_size) + '.pickle'
    heavys = []
    i = 0
    adj_list = defaultdict(dict)
    adj_list_number = 0

    # list where entry i is total edge weight in adj_list i+1 covered by heavys[i]
    prev_covered_list = []
    # list of total degrees of adj_list i
    total_degrees_list = []
    # number of vertices in adj_list i
    graph_sizes = []
    # size of 95% cover for adj_list i
    cover_sizes = []
    # intersection of cover i with cover i+1
    stables = []
    # mass(heavys[0], adj_list[i]) / total_degree(adj_list[i])
    correlations = [0.95]
    prev_heavys = None

    if naive:
        perc_f = naive_percentile
    else:
        perc_f = percentile

    with open(data_path(time, 'output_longer.txt',0,False)) as f:
        for line in f.readlines():
            line = line.strip().split()
            if line == ['snapshot'] or not line:
                i += 1
                continue
            try:
                src = line[0]
                dst = line[2]
                add_to(adj_list,src,dst)
            except:
                print line
            i += 1
            if i >= window_size:
                print "Processing graph %d" % adj_list_number
                ''' take statistics/make updates while adj_list is intact'''
                # optionally store graph of the most recent window
                if saving:
                    dump(time, 'adj_list' + str(adj_list_number) + '.pickle', adj_list, 5)
                # compute how much previous cover covers new graph
                if adj_list_number > 0:
                    prev_covered = mass(adj_list, prev_heavys)
                    prev_covered_list.append(prev_covered)
                    # compute numerator for correlations entry
                    correlation_with_start = mass(adj_list, heavys[0])
                    correlations.append(correlation_with_start)
                # store size of graph
                graph_sizes.append(len(adj_list))
                ''' make updates once percentile destroys adj_list '''
                # compute new cover, store statistics
                cover, total_degrees, last_degree = perc_f(adj_list,perc)
                # once we've found total_degree of adj_list, we divide most recent correlations entry by that
                if adj_list_number > 0:
                    correlations[-1] /= float(total_degrees)
                cover_sizes.append(len(cover))
                total_degrees_list.append(total_degrees)
                # optionally store new cover
                if saving or adj_list_number == 0:
                    heavys.append(cover)
                # compute how much covers have changed
                if adj_list_number > 0:
                    stable = len(set(cover) & set(prev_heavys))
                    stables.append(stable)
                # setup for next run
                adj_list = defaultdict(dict)
                i = 0
                adj_list_number += 1
                prev_heavys = cover

    if saving:
        print "DUMPING HEAVYS"
        dump(time, heavy_fname, heavys, 5)
    return heavys, prev_covered_list, total_degrees_list, graph_sizes, cover_sizes, stables, correlations, last_degree

def global_get_heavys(times,perc,naive,saving=False):
    if saving:
        heavy_fname = 'heavys.pickle'
    heavys = []
    adj_list = defaultdict(dict)
    adj_list_number = 0

    # list where entry i is total edge weight in adj_list i+1 covered by heavys[i]
    prev_covered_list = []
    # list of total degrees of adj_list i
    total_degrees_list = []
    # number of vertices in adj_list i
    graph_sizes = []
    # size of 95% cover for adj_list i
    cover_sizes = []
    # intersection of cover i with cover i+1
    stables = []
    # mass(heavys[0], adj_list[i]) / total_degree(adj_list[i])
    correlations = [0.95]
    prev_heavys = None

    if naive:
        perc_f = naive_percentile
    else:
        perc_f = percentile

    for time in times:
        print "ANALYZING TIME %d" % time
        if os.path.isfile('data/%d/trial_7/adj_list.pickle' % (time)):
            adj_list = load(time, 'adj_list.pickle', 7)
        else:
            with open(data_path(time, 'output_longer.txt',0,False)) as f:
                for line in f.readlines():
                    line = line.strip().split()
                    if line == ['snapshot'] or not line:
                        continue
                    try:
                        src = line[0]
                        dst = line[2]
                        add_to(adj_list,src,dst)
                    except:
                        print line
        print "DONE READING FILE"
        ''' take statistics/make updates while adj_list is intact'''
        # optionally store graph of the most recent window
        # if saving:
        #     dump(time, 'adj_list.pickle', adj_list, 7)
        # compute how much previous cover covers new graph
        if adj_list_number > 0:
            prev_covered = mass(adj_list, prev_heavys)
            prev_covered_list.append(prev_covered)
            # compute numerator for correlations entry
            correlation_with_start = mass(adj_list, heavys[0])
            correlations.append(correlation_with_start)
        # store size of graph
        graph_sizes.append(len(adj_list))
        ''' make updates once percentile destroys adj_list '''
        print "PERC_F'ing"
        # compute new cover, store statistics
        cover, total_degrees, last_degree = perc_f(adj_list,perc)
        # once we've found total_degree of adj_list, we divide most recent correlations entry by that
        if adj_list_number > 0:
            correlations[-1] /= float(total_degrees)
        cover_sizes.append(len(cover))
        total_degrees_list.append(total_degrees)
        # optionally store new cover
        if saving or adj_list_number == 0:
            heavys.append(cover)
        # compute how much covers have changed
        if adj_list_number > 0:
            stable = len(set(cover) & set(heavys[0]))
            stables.append(stable)
        # setup for next run
        adj_list = defaultdict(dict)
        i = 0
        adj_list_number += 1
        prev_heavys = cover
    if saving:
        print "DUMPING HEAVYS"
        dump(times[0], heavy_fname, heavys, 7)
    return heavys, prev_covered_list, total_degrees_list, graph_sizes, cover_sizes, stables, correlations, last_degree
    

# w1 and w2 are integers from 0 to packets/WINDOW_SIZE
# returns size of intersection of cores
def change(i, j):
    w1 = heavys[i]
    w2 = heavys[j]
    stable = len(set(w1) & set(w2))
    lost = (len(w1) - stable)/float(len(w1))
    added = (len(w2) - stable)/float(len(w1))
    sim_diff = (len(w1) - stable) + (len(w2) - stable)
    print "LOST PCT: "+str(lost)
    print "ADDED PCT: "+str(added)
    print "SIM_DIFF: "+ str(sim_diff)
    print "WINDOW_1: "+ str(len(w1))
    print "WINDOW_2: " + str(len(w2))
    print "STABLE PCT: "+ str(stable)
    return {"lost": lost, "added": added, "sim_diff": sim_diff, "lenw1": len(w1), "lenw2": len(w2), "stable": stable}

def test(key):
    extreme_val = -np.inf
    for i in range(15):
        for j in range(i+1,15):
            output = change(i,j)
            if output[key] > extreme_val:
                extreme_ij = (i,j)
                extreme_val = output[key]
    return extreme_val, extreme_ij

def total_degree(adj_list):
    return sum([-sum(adj_list[k].values()) for k in adj_list.keys()])/2
    
def mass(adj_list, cover):
    output = 0
    for el in cover:
        for v in adj_list[el].keys():
            if v in cover:
                output += adj_list[el][v]/2.
            else:
                output += adj_list[el][v]
    return output

PERC = 0.95
naive = False

def process(times, n, saving, naive, trial_number):
    # if experimental results for trial already computed and pickled, just load them, else compute
    if os.path.isfile('data/results/output%d' % trial_number):
        with open('data/results/output%d' % trial_number, 'rb') as f:
            output = pickle.load(f)
            heavys, prev_covered_list, total_degrees_list, graph_sizes, cover_sizes, stables, correlations, last_degree = output
    else:
        if len(times) == 1:
            time = times[0]
            lines = int(subprocess.check_output('wc -l ' + data_path(time, 'output_longer.txt',0,False), shell=True).split()[0])
            WINDOW_SIZE = int(lines/float(n))
            output = local_get_heavys(WINDOW_SIZE, PERC, time, naive, saving)
            heavys, prev_covered_list, total_degrees_list, graph_sizes, cover_sizes, stables, correlations, last_degree = output
            with open('data/results/output%d' % trial_number, 'wb') as f:
                pickle.dump(output, f)
        else:
            n = len(times)
            output = global_get_heavys(times, PERC, naive, saving)
            heavys, prev_covered_list, total_degrees_list, graph_sizes, cover_sizes, stables, correlations, last_degree = output
            with open('data/results/output%d' % trial_number, 'wb') as f:
                pickle.dump(output, f)
    cover_densities = []
    prev_mass_densities = []
    stable_densities = []
    if len(times) == 1:
        # replace singleton list "times" with list of windows within that trace
        times = [int(times[0]) + 49./n*i for i in range(n)]
    for i in range(n):
        print "WINDOW %d STATISTICS" % (i)
        print "TOTAL # IPs: %d" % (graph_sizes[i])
        print "TOTAL MASS: %d" % (total_degrees_list[i])
        print "COVER SIZE: %d" % (cover_sizes[i])
        cover_density = (cover_sizes[i]/float(graph_sizes[i])*100.)
        print "COVER DENSITY: %.2f%%" % cover_density
        cover_densities.append(cover_density)
        print "LAST DEGREE: %d" % last_degree
        if i > 0:
            prev_mass_density = (prev_covered_list[i-1] / float(total_degrees_list[i])*100.)
            print "PREVIOUS COVER MASS DENSITY: %.2f%%" % prev_mass_density
            prev_mass_densities.append(prev_mass_density)
            stable_density = (stables[i-1] / float(cover_sizes[i-1]) * 100.)
            print "STABLE HEAVYS: %.2f%%" % stable_density
            stable_densities.append(stable_density)
            print "CORRELATION WITH EARLIEST CORE: %.2f%%" % (correlations[i]*100.)
        print "\n"
    print "ALL STATISTICS IN ROW FORM: "
    print "PREV_COVERED_LIST: ", prev_covered_list
    print "TOTAL_DEGREES_LIST: ", total_degrees_list
    print "GRAPH_SIZES: ", graph_sizes
    print "COVER_SIZES: ", cover_sizes
    print "STABLES: ", stables
    print "LAST_DEGREE: ", last_degree
    print "CORRELATIONS: ", correlations
    print "\n"
    x_axis = [(time/100 - 1300) for time in times]
    x_axis = [(x/100)*60 + (x % 100) for x in x_axis]
    line_plot(x_axis, cover_densities, "Cover Densities", "time (min)", "density", "plots/cover_densities%d.png" % trial_number)
    line_plot(x_axis[1:], prev_mass_densities, "Persistent Mass Density", "time (min)", "density", "plots/pers_densities%d.png" % trial_number)
    line_plot(x_axis[1:], stable_densities, "Stable Set Density", "time (min)", "density", "plots/stable_densities%d.png" % trial_number)
    if n > 1:
        line_plot(x_axis, correlations, "Correlation Decay", "time (min)", "density", "plots/cor_decay%d.png" % trial_number)

def hist_plot(l, num_bins, title, xaxis, yaxis, fname):
    fig = plt.figure()
    plt.hist(l, num_bins)
    fig.suptitle(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    fig.savefig(fname, dpi=fig.dpi)
    plt.close

def line_plot(x, y, title, xaxis, yaxis, fname):
    fig = plt.figure()
    plt.plot(x,y,'ro',markersize=10,linestyle='-')
    fig.suptitle(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    fig.savefig(fname, dpi=fig.dpi)
    plt.close(fig)

def main():
    # '''125911'''
    # print "ANALYZING TIME 12:59:11"
    # # Run with 6 bins
    # print "--------------------RUNNING WITH 6 BINS, NON-NAIVE--------------------"
    # process(['125911'], 6, False, False, 1)

    # print "--------------------RUNNING WITH 6 BINS, NAIVE--------------------"
    # process(['125911'], 6, False, True, 2)

    # # Run with 2 bins
    # print "--------------------RUNNING WITH 2 BINS, NON-NAIVE--------------------"
    # process(['125911'], 2, False, False, 3)

    # print "--------------------RUNNING WITH 2 BINS, NAIVE--------------------"
    # process(['125911'], 2, False, True, 4)

    # # Run with 1 bin, for non-naive save the adj_lists
    # print "--------------------RUNNING WITH 1 BIN, NON-NAIVE--------------------"
    # process(['125911'], 1, True, False, 5)

    # print "--------------------RUNNING WITH 1 BIN, NAIVE--------------------"
    # process(['125911'], 1, False, True, 6)

    '''130000'''
    print "ANALYZING TIME 13:00:00"
    # Run with 6 bins
    print "--------------------RUNNING WITH 6 BINS, NON-NAIVE--------------------"
    process(['130000'], 6, False, False, 9)

    print "--------------------RUNNING WITH 6 BINS, NAIVE--------------------"
    process(['130000'], 6, False, True, 10)

    # Run with 2 bins
    print "--------------------RUNNING WITH 2 BINS, NON-NAIVE--------------------"
    process(['130000'], 2, False, False, 11)

    print "--------------------RUNNING WITH 2 BINS, NAIVE--------------------"
    process(['130000'], 2, False, True, 12)

    '''130000 to 140000'''
    print "ANALYZING TIMES 13:00:00 to 14:02:00"
    times = [130000,130100,130400,131600,133200,134200,135200,140200]
    print "--------------------RUNNING WITH 1 BIN, NON-NAIVE--------------------"
    process(times, 8, True, False, 7)

    print "--------------------RUNNING WITH 1 BIN, NAIVE--------------------"
    process(times, 8, False, True, 8)


if __name__ == "__main__":
    main()