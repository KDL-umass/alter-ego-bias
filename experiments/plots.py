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

def plot_dfs(dfs, ax, color):
    mean = dfs.groupby("pt.adversaries")["diff.norm"].agg("mean")
    num_adversaries = dfs.groupby("pt.adversaries")["pt.adversaries"].agg("first")
    stdev = dfs.groupby("pt.adversaries")["diff.norm"].std()
    ax.plot(num_adversaries, mean, color=color)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(alpha=0.5)
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    ax.fill_between(num_adversaries, mean-stdev, mean+stdev, facecolor=color, alpha=0.2)


if __name__=='__main__':
    filenames = ["/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-forest-fire-1-0-10.csv",
                 "/Users/kavery/workspace/non-cooperative-spillover/results/paper_csvs/new-dominating-results-forest-fire-1-0-20.csv"]
    dfs = combine_files(filenames)
    dfs = column_ops(dfs)

    fig, ax = plt.subplots(nrows=2, ncols=2)
    plot_dfs(dfs, ax[0,0], "blue")
    ax[0,0].set_xticklabels([])
    ax[0,0].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
    plot_dfs(dfs, ax[0,1], "green")
    ax[0,1].set_xticklabels([])
    ax[0,1].set_yticklabels([])
    plot_dfs(dfs, ax[1,0], "red")
    ax[1,0].set_xticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
    ax[1,0].set_yticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
    plot_dfs(dfs, ax[1,1], "black")
    ax[1,1].set_xticklabels(["0.0", "0.25", "0.5", "0.75", "1.0"])
    ax[1,1].set_yticklabels([])
    plt.show()