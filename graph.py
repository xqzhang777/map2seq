#import matplotlib.pyplot as plt
import pickle
from bokeh.plotting import ColumnDataSource, figure, output_file, save
from bokeh.models import Label
import os.path
from os import path

matches_found = "Best matches"
no_matches_found = "No matches found"

# def graph_name():
#     num = 0
#     graph_file_path = f"../Map_Files/Graphs/graph{num}.png"
#     while path.exists(graph_file_path):
#         num += 1
#         graph_file_path = f"../Map_Files/Graphs/graph{num}.png"
#     return graph_file_path


def make_graph(ids, e_vals, outputFile):
    
    #https://docs.bokeh.org/en/latest/docs/user_guide/tools.html

    output_file('{}.html'.format(outputFile))

    source = ColumnDataSource(data=dict(x=range(len(ids)),y=e_vals,ID=ids))
    top_source = ColumnDataSource(data=dict(x=[0],y=[e_vals[0]],ID=[ids[0]]))
    label = Label(x=0, y=e_vals[0], text='Best Match', x_offset=10, y_offset=-5, render_mode='canvas')
  
    TOOLTIPS = [('index','$index'),('ID','@ID'),('E-val','@y')]
   
    p = figure(width=400,height=400,tooltips=TOOLTIPS,y_axis_type='log', title='Ranked Sequences')
    p.circle('x','y',source=source)
    p.circle('x','y',source=top_source, size=10,line_color='red',fill_color='red')
    p.yaxis.axis_label = 'E-values'
    p.xaxis.axis_label = 'Rank Order'
    p.y_range.flipped = True
    p.add_layout(label)
    
    save(p)
    
    with open('{}_x.pkl'.format(outputFile),'wb') as o:
        pickle.dump(ids,o,pickle.HIGHEST_PROTOCOL)
    with open('{}_y.pkl'.format(outputFile),'wb') as o:
        pickle.dump(e_vals,o,pickle.HIGHEST_PROTOCOL)
    
    #show(p)

#    plt.clf()
#    plt.scatter(ids, e_vals, c = "blue")
#    plt.scatter(ids[0], e_vals[0], c='red')
#    plt.gca().set_xticks(plt.gca().get_xticks()[::1000])
#    plt.gca().axes.xaxis.set_ticklabels([])
#    plt.xticks(rotation = 90) 
    # plt.gca().axes.get_xaxis().set_ticks([])
#    best_match = ids[0].split('|')
#    plt.annotate(best_match[1], (ids[0], e_vals[0]))
#    plt.ylabel("E-Values")
#    plt.yscale('log')
#    plt.title("Ranked Sequences")
#    plt.gca().invert_yaxis()
#    plt.savefig(outputFile)


def parse_file(outputFile, filepath):
    # print(outputFile)
    # print(filepath)
    ids = []
    e_vals = []
    with open(filepath, 'r') as file:
        firstline = file.readline().rstrip()
        if firstline.find(matches_found) == -1:
            print(no_matches_found)
            return -1
        for line in file:
            #print(line)
            line = line.rstrip()
            #list = line.split('|')
            #print(list)
            #list[0:3] = ["|".join(list[0:3])]
            #print(list)
            #list[0] = list[0].strip()
            #list[1] = list[1].removeprefix('E-value=')
            #list[1] = float(list[1])
            #ids.append(list[0])
            #e_vals.append(list[1])
            list = line.split(' ')
            ids.append(list[0])
            e_vals.append(float(list[1]))
            
        make_graph(ids, e_vals, outputFile)
    return 1

def main():
    parse_file()


if __name__ == "__main__":
   main()