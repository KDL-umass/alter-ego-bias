import pandas as pd
import matplotlib.pyplot as plt

def combine_files(filenames):
    dfs = pd.DataFrame([])
    for filename in filenames:
        df = pd.read_csv(filename) 
        dfs = pd.concat([dfs, df])
    return dfs

def column_ops(dfs):
    dfs["bias"] = dfs["ATE.true"] - dfs["ATE.adv.gui"]
    dfs["est.diff"] = dfs["nonadv.ATE"] - dfs["ATE.adv.gui"]
    dfs["diff.norm"] = dfs["est.diff"] / dfs["nonadv.ATE"]
    dfs["pt.adversaries"] = (dfs["index"]*2)/dfs["n"]
    return dfs

def plot_dfs(dfs, ax, color, linestyle):
    mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")
    stdev = dfs.groupby("pt.adversaries")["diff.norm"].std()
    ax.plot(num_adversaries, mean, color=color, linestyle=linestyle)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(0,1)
    ax.set_xlim(0,0.5)
    ax.fill_between(num_adversaries, mean-stdev, mean+stdev, facecolor=color, alpha=0.2)

def total_plot(model, color):
    dom_filenames = []
    rand_filenames = []
    _, ax = plt.subplots(nrows=4, ncols=4)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}
    lambda1 = {"0.25":0, "0.5":1, "0.75":2, "1":3}

    for l2 in lambda2.keys():
        for l1 in lambda1.keys():
            for i in range(10, 51, 10):
                dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+model+"-"+l1+"-"+l2+"-"+str(i)+".csv")
                dom_dfs = combine_files(dom_filenames)
                dom_dfs = column_ops(dom_dfs)

                rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-"+model+"-"+l1+"-"+l2+"-"+str(i)+".csv")
                rand_dfs = combine_files(rand_filenames)
                rand_dfs = column_ops(rand_dfs)

            print(l1)
            print(l2)
            print(lambda1[l1])
            print(lambda2[l2])
            print()

            plot_dfs(dom_dfs, ax[lambda1[l1],lambda2[l2]], color, "solid")
            plot_dfs(rand_dfs, ax[lambda1[l1],lambda2[l2]], color, "dashed")

            ax[lambda1[l1],lambda2[l2]].set_xticklabels([])
            ax[lambda1[l1],lambda2[l2]].set_yticklabels([])
            if lambda1[l1]==3:
                ax[lambda1[l1],lambda2[l2]].set_xticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"])
            if lambda2[l2]==0:
                ax[lambda1[l1],lambda2[l2]].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])

            if lambda2[l2]==3:
                ax[lambda1[l1],lambda2[l2]].set_ylabel("\u03BB1="+l1,rotation=-90,labelpad=15) # \u03BB is escaped unicode for lowercase lambda
                ax[lambda1[l1],lambda2[l2]].yaxis.set_label_position("right")
            if lambda1[l1]==0:
                ax[lambda1[l1],lambda2[l2]].set_xlabel("\u03BB2="+l2)
                ax[lambda1[l1],lambda2[l2]].xaxis.set_label_position("top")

    plt.show()

if __name__=='__main__':
    total_plot("forest-fire", "green")