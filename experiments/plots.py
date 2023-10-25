import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

def combine_files(filenames):
    dfs = pd.DataFrame([])
    for filename in filenames:
        df = pd.read_csv(filename) 
        dfs = pd.concat([dfs, df])
    return dfs

def column_ops(dfs):
    dfs["bias"] = dfs["ATE.true"] - dfs["ATE.adv.gui"]
    dfs["est.diff"] = dfs["nonadv.ATE"] - dfs["ATE.adv.gui"]
    dfs["diff.norm"] = abs(dfs["est.diff"] / dfs["nonadv.ATE"])
    dfs["pt.adversaries"] = (dfs["index"]*2)/dfs["n"]
    return dfs

def synthetic_plot_dfs(dfs, ax, color, linestyle):
    mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")

    num_adversaries = pd.concat([pd.Series([0]), num_adversaries])
    mean = pd.concat([pd.Series([0]), mean])

    # print(linestyle)
    # print(mean[mean.index>0.14])

    ci = dfs.groupby("pt.adversaries")["diff.norm"].apply(lambda group: 1.96*(group.std()/np.sqrt(len(group))) )
    handle = ax.plot(num_adversaries, mean, color=color, linestyle=linestyle)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(0,1)
    ax.set_xlim(-0.01,0.5)
    ax.fill_between(num_adversaries, mean-ci, mean+ci, facecolor=color, alpha=0.2)
    return handle

def synthetic_total_plot():
    dom_filenames = []
    rand_filenames = []
    fig, ax = plt.subplots(nrows=4, ncols=3)
    fig.set_size_inches(6,6.5)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}
    model_color = {"small-world": "blue", "forest-fire": "green", "sbm": "orange"}
    model_name = {"small-world": "small world", "forest-fire": "forest fire", "sbm": "SBM"}

    for l2 in lambda2.keys():
        for i, model in enumerate(model_color.keys()):
            dom_filenames = []
            rand_filenames = []
            for j in range(10, 51, 10):
                dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+model+"-"+"0.25"+"-"+l2+"-"+str(j)+".csv")
                dom_dfs = combine_files(dom_filenames)
                dom_dfs = column_ops(dom_dfs)

                rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-"+model+"-"+"0.25"+"-"+l2+"-"+str(j)+".csv")
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


def individual_total_plot(name):
    dom_filenames = []
    rand_filenames = []
    fig, ax = plt.subplots(ncols=4)
    fig.set_size_inches(7.5,2)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}

    for l2 in lambda2.keys():
        dom_filenames = []
        rand_filenames = []
        for j in range(10, 21, 5):
            dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+name+"-1-"+l2+"-"+str(j)+".csv")
            dom_dfs = combine_files(dom_filenames)
            if name=="facebook":
                dom_dfs["n"] = 4039
            dom_dfs = column_ops(dom_dfs)

        for j in range(10, 11, 10):
            rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-"+name+"-1-"+l2+"-"+str(j)+".csv")
            rand_dfs = combine_files(rand_filenames)
            if name=="facebook":
                rand_dfs["n"] = 4039
            rand_dfs = column_ops(rand_dfs)

        worst_handle = synthetic_plot_dfs(dom_dfs, ax[lambda2[l2]], "black", "solid")
        rand_handle = synthetic_plot_dfs(rand_dfs, ax[lambda2[l2]], "black", "dashed")

        ax[lambda2[l2]].set_xticklabels([])
        ax[lambda2[l2]].set_yticklabels([])
        # if lambda2[l2]==3:
        ax[lambda2[l2]].set_xticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"])
        if lambda2[l2]==0:
            ax[lambda2[l2]].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
        
        ax[lambda2[l2]].set_xlabel("\u03BB2="+l2,labelpad=15)
        ax[lambda2[l2]].xaxis.set_label_position("top")

    leg = fig.legend([worst_handle[0], rand_handle[0]], ["worst-case", "random"], loc="center right",)
    fig.supxlabel('Alter ego fraction of network')
    fig.supylabel('Normalized bias')
    plt.show()


def upper_bound_total_plot(name):
    dom_filenames = []
    fig, ax = plt.subplots(ncols=4)
    fig.set_size_inches(7.5,2)
    lambda2 = {"0":0, "0.1":1, "0.5":2, "1":3}

    for l2 in lambda2.keys():
        dom_filenames = []
        for j in range(10, 11, 10):
            dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+name+"-1-"+l2+"-"+str(j)+".csv")
            dom_dfs = combine_files(dom_filenames)
            if name=="facebook":
                dom_dfs["n"] = 4039
            dom_dfs = column_ops(dom_dfs)
        
        worst_handle = synthetic_plot_dfs(dom_dfs, ax[lambda2[l2]], "black", "solid")

        ax[lambda2[l2]].set_xticklabels([])
        ax[lambda2[l2]].set_yticklabels([])
        ax[lambda2[l2]].set_xticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"])
        if lambda2[l2]==0:
            ax[lambda2[l2]].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
        
        ax[lambda2[l2]].set_xlabel("\u03BB2="+l2,labelpad=15)
        ax[lambda2[l2]].xaxis.set_label_position("top")

    leg = fig.legend([worst_handle[0]], ["worst-case", "random"], loc="upper center", ncol=2,)
    fig.supxlabel('Alter ego fraction of network')
    fig.supylabel('Normalized bias')
    plt.show()


def spillover_plot_dfs(dfs, ax, color, linestyle):
    mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")
    ci = dfs.groupby("pt.adversaries")["diff.norm"].apply(lambda group: 1.96*(group.std()/np.sqrt(len(group))) )
    if linestyle == "solid":
        mean = (mean - num_adversaries)/mean
    else:
        mean = (mean - num_adversaries/2)/mean
    
    handle = ax.plot(num_adversaries, mean, color=color, linestyle=linestyle)
    ax.fill_between(num_adversaries, mean-ci, mean+ci, facecolor=color, alpha=0.2)
    return handle


def spillover_plot():
    fig, ax = plt.subplots(nrows=1)
    fig.set_size_inches(2.5,3.5)
    legend_handles = []

    model_name = {"small-world": "small world", "forest-fire": "forest fire", "sbm": "SBM"}
    model_color = {"small-world": "blue", "forest-fire": "green", "sbm": "orange"}
    model_id = {"small-world": 0, "forest-fire": 1, "sbm": 2}
    for model in model_name.keys():
        dom_filenames = []
        rand_filenames = []
        for j in range(10, 51, 10):
            dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+model+"-1-"+"1"+"-"+str(j)+".csv")
            dom_dfs = combine_files(dom_filenames)
            dom_dfs = column_ops(dom_dfs)

            rand_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-random-results-"+model+"-1-"+"1"+"-"+str(j)+".csv")
            rand_dfs = combine_files(rand_filenames)
            rand_dfs = column_ops(rand_dfs)
        
        worst_handle = spillover_plot_dfs(dom_dfs, ax, model_color[model], "solid")
        rand_handle = spillover_plot_dfs(rand_dfs, ax, model_color[model], "dashed")
        legend_handles.append(worst_handle[0])
        legend_handles.append(rand_handle[0])

    dom_filenames = []
    for j in range(10, 11, 10):
        dom_filenames.append("/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-"+"fan"+"-1-"+"1"+"-"+str(j)+".csv")
        dom_dfs = combine_files(dom_filenames)
        dom_dfs = column_ops(dom_dfs)

    worst_net_handle = spillover_plot_dfs(dom_dfs, ax, "black", "solid")
    legend_handles.append(worst_net_handle[0])

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(0,1)
    ax.set_xlim(0,0.5)
    ax.set_xticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"])
    ax.set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])

    leg = fig.legend(legend_handles, ["small world (worst-case)", "small world (random)",
                                      "forest fire (worst-case)", "forest fire (random)", "SBM (worst-case)", "SBM (random)",
                                      "worst-case network"], 
                     loc="upper center", ncol=2,)
    fig.supxlabel('Alter ego fraction of network')
    fig.supylabel('Fraction of normalized bias due to peer effects')
    plt.show()


if __name__=='__main__':
    # upper_bound_total_plot("fan")
    # spillover_plot()
    # synthetic_total_plot()
    individual_total_plot("facebook")
    # spillover_plot()