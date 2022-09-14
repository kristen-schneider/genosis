import argparse
import collections
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    parser.add_argument('--samples_dir')
    parser.add_argument('--popA')
    parser.add_argument('--popB')
    parser.add_argument('--pop_file')
    parser.add_argument('--png')
    return parser.parse_args()

def main():
    args = get_args()
    A_file = args.samples_dir + args.popA + '.samples'
    popA_samples = read_samples(A_file)
    #B_file = args.samples_dir + args.popB + '.samples'
    #popB_samples = read_samples(B_file)
    
    plink_dict = read_plink_dict(args.plink)
    [pop_summary_dict, pop_ID_name_dict, pop_ID_super_name_dict] = read_pop_file(args.pop_file)
    
    A_all_pop_data = defaultdict(dict)
    for X_pop in pop_summary_dict:
        X_file = args.samples_dir + X_pop + '.samples'
        print(X_file)
        X_super_pop = pop_ID_super_name_dict[X_pop]
        popX_samples = read_samples(X_file)
        A_X_data = make_vp_data(popA_samples, popX_samples, plink_dict)
        A_all_pop_data[X_super_pop][X_pop] = A_X_data

    # sort by super population
    all_plot_data = dict()
    for sp in dict(sorted(A_all_pop_data.items())):
        for p in A_all_pop_data[sp]:
            all_plot_data[p] = A_all_pop_data[sp][p]
            print(sp, p, len(A_all_pop_data[sp][p]))
    plot_pop_violinplot(args.popA, all_plot_data, pop_ID_name_dict, pop_ID_super_name_dict, args.png)

def plot_pop_violinplot(popA, plot_data, pop_ID_name_dict, pop_ID_super_name_dict, png):
    # organizing data to plot
    violin_data = [plot_data[popX] for popX in plot_data]
    #violin_data = plot_data

    # plotting data
    png_name = png + popA + '_all.relatedness.png'
    plt.figure(figsize=(30, 20))
    ax1 = plt.subplot(111)
    png_title = 'Plink reported relatedness for\n' + pop_ID_name_dict[popA] + ' samples'
    ax1.set_title(png_title, fontsize=35)
    vp = ax1.violinplot(violin_data)
    # colors
    super_pops = []
    for pop_ID in pop_ID_super_name_dict:
        if pop_ID_super_name_dict[pop_ID] not in super_pops:
            super_pops.append(pop_ID_super_name_dict[pop_ID])
    colors = ['firebrick', 'olivedrab', 'steelblue', 'gold', 'rebeccapurple']
    ordered_colors = []
    for v in plot_data:
        v_sp = pop_ID_super_name_dict[v]
        sp_i = super_pops.index(v_sp)
        ordered_colors.append(colors[sp_i])
    p = 0
    for v in vp['bodies']:
        v.set_facecolor(ordered_colors[p])
        p += 1
    # legend
    legend = []
    for i in range(len(colors)):
        legend.append(mpatches.Patch(color=colors[i], label=super_pops[i]))
    ax1.legend(handles=legend, fontsize=20)
    # formatting
    ax1.set_xlabel('population', fontsize=30)
    ax1.set_ylabel('plink similarity', fontsize=30)
    x_label = []
    for pop_id in plot_data.keys():
        x_label.append(pop_ID_name_dict[pop_id])
    ax1.set_xticks(np.arange(1, len(plot_data)+1), labels=x_label,  rotation=90, fontsize=20)
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    plt.subplots_adjust(bottom=0.19)
    plt.savefig(png_name)

def make_vp_data(popA_samples, popX_samples, plink_dict):
    popA_popX = []
    for sample_A in popA_samples:
        curr_sample = plink_dict[sample_A]
        for sample_X in popX_samples:
            try:
                A_X = curr_sample[sample_X]
                popA_popX.append(A_X)
            except KeyError:
                continue
            
    return popA_popX


def read_pop_file(pop_summary):
    # pop_ID, samples
    summary_dict = dict()
    # pop_ID, pop_name
    pop_ID_name_dict = dict()
    # pop_ID, pop_super_name
    pop_ID_super_name_dict = dict()

    f = open(pop_summary)
    for line in f:
        L = line.strip().split(',')
        sample_ID = L[0]
        pop_ID = L[1]
        pop_name = L[2]
        pop_super_ID = L[3]
        pop_super_ID_name = L[4]

        pop_ID_name_dict[pop_ID] = pop_name
        pop_ID_super_name_dict[pop_ID] = pop_super_ID_name
        #pop_ID_super_ID_dict[pop_ID] = pop_super_ID 
        #pop_super_ID_name_dict[pop_super_ID] = pop_super_ID_name
        try:
            summary_dict[pop_ID].append(sample_ID)
        except KeyError:
            summary_dict[pop_ID] = [sample_ID]
    f.close()
    return summary_dict, pop_ID_name_dict, pop_ID_super_name_dict

def read_plink_dict(plink_file):
    plink_dict = defaultdict(dict)
    f = open(plink_file, 'r')
    for h in range(3):
        f.readline()
    for line in f:
        L = line.strip().split()
        SID1 = L[0]
        SID2 = L[1]
        dist = float(L[2])
        # add self to dictionary
        #plink_dict[SID1][SID1] = 1.0
        plink_dict[SID1][SID2] = dist
        plink_dict[SID2][SID1] = dist

    f.close()
    return plink_dict

def read_samples(samples_file):
    samples_list = []
    f = open(samples_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        SID = L[0]
        samples_list.append(SID)
    f.close()
    return samples_list


if __name__ == '__main__':
    main()
