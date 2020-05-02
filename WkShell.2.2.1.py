__author__ = "Nilesh Kumar"
__email__ = "nilesh.iiita@gmail.com"
__copyright__ = "Copyright 2020, MIT License"
__credits__ = ["Nilesh Kumar"]
__license__ = "MIT License"
__version__ = "2.2.1"
__maintainer__ = "Nilesh Kumar"
__status__ = "Production"



"""
Requires-Python: >=3
Requires-Python:Networkx >=2.4
Requires-Python:Networkx >=2.4
Requires-Python:pandas

Takes csv (Comma-separated values) files as input network (edges list).
Updated on : May 1 2020
"""

import sys, os, re, time
import networkx as nx
from tqdm import tqdm
from collections import defaultdict
import pandas as pd


Input_file_name = "PPIN1.csv"
Java_ref_directory = "Wk_shell_files"#"Java_compatibale_files"

if not os.path.exists(Java_ref_directory):
    os.makedirs(Java_ref_directory)
Java_ref_directory += "/"

print(sys.version)


def remove_self_loops(File):

    G = nx.read_edgelist(File, delimiter=",")
    print(nx.info(G))
    print("removing self loops")
    G.remove_edges_from(nx.selfloop_edges(G)) 
    print(nx.info(G))
    
    File = Java_ref_directory + File+".unsloop.csv"
    print("writing network to",File)
    
    #text_box.delete(1.0, "end-1c")
    #text_box.insert("end-1c", "You win!")

    nx.write_edgelist(G, File, delimiter=',', data=False)
    #Ecount = G.number_of_edges()
    return File, G

def numerical_network(Filename):
    Data = open(Filename).read().splitlines()
    List = []
    ref_dic_nodes = {}
    ref_dic_nodes_rev = {}
    for i in Data:
        i = i.rstrip()
        if not len(i):continue
        a,b = i.split(",")
        if a not in List: List.append(a)
        if b not in List: List.append(b)

    fh = open(Java_ref_directory+"Gene_list.txt", "w")
    print(f"Total number of nodes : {len(List)}")
    for n in tqdm(range(len(List))):
        i = List[n]
        print(i, file=fh)
        I = n+1
        ref_dic_nodes[i] = str(I)
        ref_dic_nodes_rev[str(I)] = i
    fh.close()
    #G = nx.Graph(name="Numnet")
    print(f"Indexed gene list is written..")

    fh = open(Java_ref_directory+"edges_num.csv", "w")

    for i in Data:
        i = i.rstrip()
        if not len(i):continue
        a,b = i.split(",")
        a, b = ref_dic_nodes[a], ref_dic_nodes[b]
        print(a, b, sep=",", file=fh)    
    print(f"Numerical network is written..")

    fh.close()


    return ref_dic_nodes, ref_dic_nodes_rev
    

def wkshell(ref_dic_nodes, alpha = 0.5):
    old_min = -2147483648
    shl = 1;
    #alpha = 0.5
    filename = Java_ref_directory+"edges_num.csv"

    ##1-Read edges into the Graph
    tic = time.process_time()


    #print('Start of #1 %f sec' % (time.process_time() - tic))

    G = nx.read_edgelist(filename,delimiter = ",")
    Nodes = G.number_of_nodes()
    shells = []

    Degree = G.degree()

    def weight(G,i,w):
        return Degree[str(i)] + Degree[str(w)]
    wgt = 0



    W_e = {}
    N_w = {}
    for n in G:
        N_w[n] = 0
        for w in G[n]:
            W = weight(G,n,w)
            L = [n,w]
            s = "_".join([str(L[0]),str(L[1])])
            s1 = "_".join([str(L[1]),str(L[0])])
            if not(s in W_e or s1 in W_e):
                W_e[s] = [W,0]
    List = defaultdict(list)
    for i in W_e:
        a,b  = i.split("_")
        List[a].append(i)
        List[b].append(i)

    #print('Start of #3 %f sec' % (time.process_time() - tic))
    while G.number_of_nodes():
        Time = (time.process_time() - tic)
        #print('Start of #3.1 %f sec' % (time.process_time() - tic))
        for i in range(1,len(N_w)+1):
            i = str(i)
            we = 0
            if i in G:
                for w in G.neighbors(str(i)):
                    if w in G:
                        we += weight(G,i,w)
                    if str(i)+"_"+str(w) in W_e:
                        W_e[str(i)+"_"+str(w)][0] = weight(G,i,w)
                        W_e[str(i)+"_"+str(w)][1] = 1

                for j in List[i]:
                    a,b = j.split("_")  
                    if((a ==i or b ==i ) and int(W_e[j][1])==0):
                        we += W_e[j][0];

                N_w[i] = int(alpha*Degree[str(i)] + (1-alpha)*we)

        #print('Start of #4 %f sec' % (time.process_time() - tic))
        #Find lowest weight in N_w
        Min = min(N_w.items(), key=lambda x: x[1])[1]
        #5- Remove it in a shell
        sh = "";
        for i in N_w:
            if(N_w[i] == Min):
                sh += str(i)+","
                G.remove_node(str(i))
                for j in W_e:
                    n,w = j.split("_")
                    if(n == i or w == i):
                        W_e[j][1] = 0;
        ss = "";
        if(Min <= old_min):
            Min = old_min;
            ss = shells[len(shells)-1]
            shells.pop(len(shells)-1);
            ss+=sh;
            sh=ss;

        shells.append(sh+"("+str(Min)+")");

        sh = "";
        for i in N_w:
            N_w[i] = 2147483647
        old_min = Min
        #print((time.process_time() - tic),(time.process_time() - tic) - Time,G.number_of_nodes(),(Nodes-G.number_of_nodes())*100./Nodes,"%")
        msg = (time.process_time() - tic),(time.process_time() - tic) - Time,G.number_of_nodes(),(Nodes-G.number_of_nodes())*100./Nodes,"%"
        print(msg)
        sys.stdout.write("\033[F")

        
        
    with open(Java_ref_directory+ 'test_out.txt', 'w') as f:
        for s in shells:
            f.write("("+str(shl)+")-"+s+'\n')
    #fhn = open(Java_ref_directory+"Wk-Shell_numerical.csv", "w"); print("Node,Shell", file= fhn)
    nstup = []
    tup = []

    for s in shells:
        numbers = list(set(re.findall(r'[0-9]+',s)))
        shell_set = set(re.findall(r'\(\d*\)', s))
        #if len(shell_set) > 1:
        #    print(s)
        shell_set = list(shell_set)[0].replace("(","").replace(")","")
        #print(shell_set)
        numbers.remove(shell_set)
        #print(numbers, shell_set, s)
        for n in numbers:
            nstup.append((n,int(shell_set)))
            tup.append((ref_dic_nodes[n],int(shell_set)))

    dfn = pd.DataFrame(nstup, columns = ['Node_num', 'Shell']) 
    df = pd.DataFrame(tup, columns = ['Node', 'Shell']) 

    return df, dfn
    #fhn.close()



if __name__ == "__main__":
    File = Input_file_name#"PPIN1.csv"#"Network_1000.csv"
    #File = tkinter.filedialog.askopenfilename(initialfile = "ai1.csv", title = "Select network file", filetypes = (("CSV files","*.csv"),("all files","*.*")))
    #File = os.path.basename(File)
    print(File)
    Filename, G = remove_self_loops(File)
    ref_dic_nodes, ref_dic_nodes_rev = numerical_network(Filename)
    df, dfn = wkshell(ref_dic_nodes_rev)
    df = df.sort_values(by=['Shell'], ascending=True)
    dfn = dfn.sort_values(by=['Shell'], ascending=True)
    List = df["Shell"].values.tolist()
    Ulist = list(set(List))
    Ulist.sort()
    Dic_index = {}
    Indexs = []
    for i in List:
        I = Ulist.index(i)+1
        Indexs.append(I)
    #print(Indexs)
    df["Shell_number"] = Indexs
    df = df.sort_values(by=['Shell_number'], ascending=False)
    

    dfn.to_csv(Java_ref_directory+"Wkshell_num.csv", index=False)
    df["Wk_mean-normalized"] = (df["Shell_number"]-df["Shell_number"].mean())/df["Shell_number"].std()
    df["Wk_min-max-normalization"] = (df["Shell_number"]-df["Shell_number"].min())/(df["Shell_number"].max()-df["Shell_number"].min())
    saveFilePath = Java_ref_directory+"Wkshell.csv"

    #saveFilePath  =  tkinter.filedialog.asksaveasfilename(initialfile = "Wkshell.csv", initialdir = Java_ref_directory, title = "Save as",filetypes = (("CSV files","*.csv"),("all files","*.*")))
    #df.to_csv(Java_ref_directory+"Wkshell.csv", index=False)
    print(f"Writing to CSV file...")
    df.to_csv(saveFilePath, index=False)

    print(f"Writing to Microsoft excle (xlsx sheet) file...")
    df.to_excel(saveFilePath.replace(".csv",".xlsx"), sheet_name='Weighted k-shell decomposition', index=False, float_format="%.2f")

    print(f"Writing Html table")
    df.to_html(saveFilePath.replace(".csv",".html"), index=False, float_format="%.2f", justify= "justify-all", notebook = True, show_dimensions = True)
    
    print(f"Data written to disk..")

    print(df.describe())
    ##Plots
    #ax = df["Shell_number"].plot.hist(bins=12, alpha=0.5)
    reseponse = input("Enter Y/yes to plot distribution or N/No to exit\n(y/n):").upper()
    
    if reseponse == "Y" or reseponse == "YES":
        X = [i for i in range(df.shape[0])]
        df["Nodes"] = X
        ##ax = df.plot.hexbin(y = "Nodes",x = "Shell_number", gridsize=20, colormap = 'GnBu')
        ax = df.plot.scatter(x = "Nodes",y = "Shell_number", c = "Shell", colormap='viridis')
        fig = ax.get_figure()
        fig.savefig(saveFilePath.replace(".csv",".png"))
    else:
        print(f"You have saved lots of time..\nPlot skipped")

  
    

    




    
