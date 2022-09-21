import argparse
import distance_calculations
import read_encoding
import matplotlib.colors as c
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    parser.add_argument('--query')
    parser.add_argument('--dad')
    parser.add_argument('--mom')
    parser.add_argument('--better0')
    parser.add_argument('--better1')
    return parser.parse_args()

def main():
    args = get_args()
    encoding_dict = read_encoding.read_encoding_file(args.encoding_file)
    
    #trio_dict_hap0 = {args.query+'_'+str(0): encoding_dict[args.query+'_'+str(0)],
    #                    args.dad+'_'+str(0): encoding_dict[args.dad+'_'+str(0)],
    #                    args.dad+'_'+str(1): encoding_dict[args.dad+'_'+str(1)],
    #                    args.mom+'_'+str(0): encoding_dict[args.mom+'_'+str(0)],
    #                    args.mom+'_'+str(1): encoding_dict[args.mom+'_'+str(1)]}
    #trio_data_hap0 = make_plot_data(trio_dict_hap0)
    #trio_dict_hap1 = {args.query+'_'+str(1): encoding_dict[args.query+'_'+str(1)],
    #                    args.dad+'_'+str(0): encoding_dict[args.dad+'_'+str(0)],
    #                    args.dad+'_'+str(1): encoding_dict[args.dad+'_'+str(1)],
    #                    args.mom+'_'+str(0): encoding_dict[args.mom+'_'+str(0)],
    #                    args.mom+'_'+str(1): encoding_dict[args.mom+'_'+str(1)]}
    #trio_data_hap1 = make_plot_data(trio_dict_hap1)

    #plot_data = make_plot_data(encoding_dict)
    
    self_hap0 = encoding_dict[args.query+'_'+str(0)]
    self_hap1 = encoding_dict[args.query+'_'+str(1)]
    dad_hap0 = encoding_dict[args.dad+'_'+str(0)]
    dad_hap1 = encoding_dict[args.dad+'_'+str(1)]
    mom_hap0 = encoding_dict[args.mom+'_'+str(0)]
    mom_hap1 = encoding_dict[args.mom+'_'+str(1)]
    better_hap0 = encoding_dict[args.better0]
    better_hap1 = encoding_dict[args.better1]

    trio_encodings_0 = [self_hap0, dad_hap0, better_hap0, mom_hap0]
    trio_samples_0 = ['self0', 'dad0', 'better0', 'mom0']
    trio_encodings_1 = [self_hap1, mom_hap1, better_hap1, dad_hap1]
    trio_samples_1 = ['self1', 'mom1', 'better1', 'dad1']
    
    # compute distances
    compute_distances(trio_encodings_0, trio_samples_0, 0)
    compute_distances(trio_encodings_1, trio_samples_1, 1)
    # plot
    print('plotting...')
    plot_encodings(trio_encodings_0, trio_samples_0, 0)
    plot_encodings(trio_encodings_1, trio_samples_1, 1)

def compute_distances(trio_encodings, trio_samples, hap):
    print('computing distances for hap: ', hap)
    print(f"\nid\tid\t?\ted")
    start = 0
    end = len(trio_encodings[0])-1
    for s in range(len(trio_encodings)):
        sv = distance_calculations.shared_variants(trio_encodings[0][start:end], trio_encodings[s][start:end])
        ed = distance_calculations.euclidean_distance(trio_encodings[0][start:end], trio_encodings[s][start:end])
        L = trio_samples[0], trio_samples[s], sv, ed
        print(f"{L[0]}\t{L[1]}\t{L[2]}\t{L[3]}")
    

def plot_encodings(trio_encodings, trio_samples, hap):
    plt.figure(figsize=(60, 10), dpi=200)
    color_scale = 'Pastel1'
    cmap = c.ListedColormap(['mistyrose', 'maroon'])
    #print(trio_encodings[0])
    #ax1 = plt.subplot(1, 1, 1)
    for sp in range(len(trio_encodings)):
        ax_sp = plt.subplot(len(trio_encodings), 1, sp+1)
        if sp == 0:
            ax_sp.set_title('Segment 25, Haplotype '+str(hap), fontsize=40)
        ax_sp.imshow([trio_encodings[sp]], cmap=cmap)#colors[sp])
        if sp == len(trio_encodings)-1:
            ax_sp.set_xticks(range(0, len(trio_encodings[0]), 500), pad=30)
            ax_sp.set_xlabel('basepair', fontsize=30) 
        else: ax_sp.set_xticks([])
        ax_sp.set_yticks([])
        ax_sp.set_ylabel(trio_samples[sp], rotation=90, labelpad=30, fontsize=30) 
        ax_sp.spines.top.set_visible(False)
        ax_sp.spines.bottom.set_visible(False)
        ax_sp.spines.left.set_visible(False)
        ax_sp.spines.right.set_visible(False)
        ax_sp.set_aspect('auto') 

    plt.subplots_adjust(wspace=0.4, hspace=0.6, bottom=0.10, top=0.85, left=0.05, right=0.95)
    plt.savefig('encoding.'+str(hap)+'.png')

def make_plot_data(encoding_dict):
    plot_data = []
    for s in encoding_dict:
        plot_data.append(encoding_dict[s])
    return plot_data

if __name__ == '__main__':
    main()
