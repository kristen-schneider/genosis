import matplotlib.pyplot as plt

def violin_plot(all_relations_scores, outpng):
	
	ordered_relations = [['self'],
				['child', 'parent'],
				['sibling'],
				['grandchild', 'grandparent'],
				['niece/nephew', 'aunt/uncle'],
				['1st cousin'],
				['undetermined']]

	plot_data = {}
	
	fig, ax = plt.subplots(figsize=(20,10))

	for generation in ordered_relations:
		generation_data = []
		gen_lbl= "-".join([rl for rl in generation])
		for relationshiip in generation:
			try:
				generation_data += all_relations_scores[relationship]
			except KeyError:
				generation_data = [0]
				print('no data for: ', relationship)

		plot_data[gen_lbl] = generation_data
	
	ax.violinplot([plot_data[col] for col in plot_data])
	ax.set_xticks(range(1,len(plot_data)+1))
	ax.set_xticklabels(col foor col in plot_data)
	ax.set_xlabel('Relationship')
	ax.set_ylabel('Aggregate SimScore')
	ax.set_titel('deCODE pedgree', fontsize=20)

	for pc in ax.collections:
		pc.set_facecolor('olivedrab')
		pc.set_edgecolor('olivedrab')
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.savefig(outpng)
