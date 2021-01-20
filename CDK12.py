import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors



start = int(37617764)
end = int(37721160)

f_table = pd.read_table("CDK12.txt",index_col = 0)

for i in range(start,end):
    if i not in f_table.index:
        f_table.loc[i] = 0
        
f_table.sort_index(inplace = True)



fig, ax = plt.subplots(figsize=(20, 2.5))
font = {'family' : 'Times New Roman', 'Weight' : 'normal', 'Size': 15,}
sns.heatmap(f_table.T,cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#A6A4A4', '#FF0000']))

ax.set_title("CDK12:Chr17:37617764-37712160",font)
ax.set_yticks([])
fig.savefig("cdk12-editing.pdf")
                                         
