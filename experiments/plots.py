import pandas as pd
import numpy as np
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

def synthetic_plot_dfs(dfs, ax, color, linestyle):
    mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")
    ci = dfs.groupby("pt.adversaries")["diff.norm"].apply(lambda group: 1.96*(group.std()/np.sqrt(len(group))) )
    handle = ax.plot(num_adversaries, mean, color=color, linestyle=linestyle)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(0,1)
    ax.set_xlim(0,0.5)
    ax.fill_between(num_adversaries, mean-ci, mean+ci, facecolor=color, alpha=0.2)
    return handle

def synthetic_total_plot():
    dom_filenames = []
    rand_filenames = []
    fig, ax = plt.subplots(nrows=4, ncols=3)
    fig.set_size_inches(6,6.5)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}
    # lambda1 = {"0.25":0, "0.5":1, "0.75":2, "1":3}
    model_color = {"small-world": "blue", "forest-fire": "green", "sbm": "orange"}
    model_name = {"small-world": "small world", "forest-fire": "forest fire", "sbm": "SBM"}

    for l2 in lambda2.keys():
        for i, model in enumerate(model_color.keys()):
            for j in range(10, 51, 10):
                dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+model+"-"+"1"+"-"+l2+"-"+str(j)+".csv")
                dom_dfs = combine_files(dom_filenames)
                dom_dfs = column_ops(dom_dfs)

                rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-"+model+"-"+"1"+"-"+l2+"-"+str(j)+".csv")
                rand_dfs = combine_files(rand_filenames)
                rand_dfs = column_ops(rand_dfs)

            print(l2)
            worst_handle = synthetic_plot_dfs(dom_dfs, ax[lambda2[l2],i], model_color[model], "solid")
            rand_handle = synthetic_plot_dfs(rand_dfs, ax[lambda2[l2],i], model_color[model], "dashed")

            ax[lambda2[l2],i].set_xticklabels([])
            ax[lambda2[l2],i].set_yticklabels([])
            if lambda2[l2]==3:
                ax[lambda2[l2],i].set_xticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"])
            if model=="small-world":
                ax[lambda2[l2],i].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])

            if lambda2[l2]==0:
                ax[lambda2[l2],i].set_xlabel(model_name[model]) # \u03BB is escaped unicode for lowercase lambda
                ax[lambda2[l2],i].xaxis.set_label_position("top")
            if model=="sbm":
                ax[lambda2[l2],i].set_ylabel("\u03BB2="+l2,rotation=-90,labelpad=15)
                ax[lambda2[l2],i].yaxis.set_label_position("right")

    leg = fig.legend([worst_handle[0], rand_handle[0]], ["worst-case", "random"], loc="upper center", ncol=2,)
    leg.legendHandles[0].set_color("black")
    leg.legendHandles[1].set_color("black")
    fig.supxlabel('Alter ego fraction of network')
    fig.supylabel('Normalized bias')
    plt.show()

def facebook_plot_dfs(dfs, ax, color, random):
    print(dfs)
    # mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    # num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")
    # ci = dfs.groupby("pt.adversaries")["diff.norm"].apply(lambda group: 1.96*(group.std()/np.sqrt(len(group))) )

    data = dfs.groupby("pt.adversaries")["diff.norm"].agg(list)
    print(len(dfs["diff.norm"]))
    print(len(dfs["pt.adversaries"]))
    if random:
        handle = ax.boxplot(data, whiskerprops=dict(color=color),
                            boxprops=dict(color=color),
                            capprops=dict(color=color),
                            medianprops=dict(color=color),
                            flierprops = dict(color=color, markeredgecolor=color),
                            positions=[1+0.15, 2+0.15, 3+0.15, 4+0.15, 5+0.15], 
                            labels=["0.0005", "0.0010", "0.0015", "0.0020", "0.0.0025"], widths=0.2)
    else:
        handle = ax.boxplot(data, whiskerprops=dict(color=color),
                            boxprops=dict(color=color),
                            capprops=dict(color=color),
                            medianprops=dict(color=color),
                            flierprops = dict(color=color, markeredgecolor=color),
                            positions=[1-0.15, 2-0.15, 3-0.15, 4-0.15, 5-0.15], labels=["0.0005", "0.0010", "0.0015", "0.0020", "0.0.0025"], widths=0.2) 
    ax.set_xticks([1, 2, 3, 4, 5])
    # ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(-0.01,0.06)
    ax.set_xlim(0.5, 5.5)
    # ax.fill_between(num_adversaries, mean-ci, mean+ci, facecolor=color, alpha=0.2)
    # print(mean)
    # print(ci)
    return handle

def facebook_total_plot():
    dom_filenames = []
    rand_filenames = []
    fig, ax = plt.subplots(nrows=4, ncols=1)
    fig.set_size_inches(3,6.5)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}

    for l2 in lambda2.keys():
            for j in range(10, 11, 10):
                dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-facebook-"+"1"+"-"+l2+"-"+str(j)+".csv")
                dom_dfs = combine_files(dom_filenames)
                dom_dfs["n"] = 4039
                dom_dfs = column_ops(dom_dfs)

                rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-facebook-"+"1"+"-"+l2+"-"+str(j)+".csv")
                rand_dfs = combine_files(rand_filenames)
                rand_dfs["n"] = 4039
                rand_dfs = column_ops(rand_dfs)

            worst_handle = facebook_plot_dfs(dom_dfs, ax[lambda2[l2]], "red", False)
            rand_handle = facebook_plot_dfs(rand_dfs, ax[lambda2[l2]], "blue", True)

            ax[lambda2[l2]].set_xticklabels([])
            # ax[lambda2[l2]].set_yticklabels([])
            if lambda2[l2]==3:
                ax[lambda2[l2]].set_xticklabels(["0.0005", "0.0010", "0.0015", "0.0020", "0.0025"])
            # ax[lambda2[l2]].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
            
            ax[lambda2[l2]].set_ylabel("\u03BB2="+l2,rotation=-90,labelpad=15) # \u03BB is escaped unicode for lowercase lambda
            ax[lambda2[l2]].yaxis.set_label_position("right")

    leg = fig.legend([worst_handle["boxes"][0], rand_handle["boxes"][0]], ["worst-case", "random"], loc="upper center", ncol=2)
    # leg.legendHandles[0].set_color("black")
    # leg.legendHandles[1].set_color("black")
    fig.supxlabel('Alter ego fraction of network')
    fig.supylabel('Normalized bias')
    plt.show()

if __name__=='__main__':
    facebook_total_plot()