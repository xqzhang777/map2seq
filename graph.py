import matplotlib.pyplot as plt
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
    plt.clf()
    plt.scatter(ids, e_vals, c = "blue")
    plt.scatter(ids[0], e_vals[0], c='red')
    plt.gca().set_xticks(plt.gca().get_xticks()[::1000])
    plt.gca().axes.xaxis.set_ticklabels([])
    plt.xticks(rotation = 90) 
    # plt.gca().axes.get_xaxis().set_ticks([])
    best_match = ids[0].split('|')
    plt.annotate(best_match[1], (ids[0], e_vals[0]))
    plt.ylabel("E-Values")
    plt.yscale('log')
    plt.title("Ranked Sequences")
    plt.gca().invert_yaxis()
    plt.savefig(outputFile)


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