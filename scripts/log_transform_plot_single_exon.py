import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from matplotlib.lines import Line2D
import seaborn as sns

animal_names = ["gibbon", "gorilla", "human", "bonobo", "chimpanzee", "borangutan", "sorangutan"]

for name in animal_names:
    df = pd.read_csv(f"../data/{name}_single_exon.tsv", sep="\t")

    refined_df = df[
        ['exon_splicing', 'exon_order',
         'upstream_len_in_scfr', 'downstream_len_in_scfr']
    ]

    upstream = refined_df[['upstream_len_in_scfr']]
    downstream = refined_df[['downstream_len_in_scfr']]

    # remove zeros
    upstream = upstream[upstream['upstream_len_in_scfr'] != 0]
    downstream = downstream[downstream['downstream_len_in_scfr'] != 0]

    # log10 transformation
    upstream['transformed_upstream'] = np.log10(upstream['upstream_len_in_scfr'])
    downstream['transformed_downstream'] = np.log10(downstream['downstream_len_in_scfr'])

    # create the exact dataframe names you want
    globals()[f"{name}_only_upstream"] = upstream
    globals()[f"{name}_only_downstream"] = downstream

conditions = []

for animal in animal_names:
    up_df = globals()[f"{animal}_only_upstream"]
    down_df = globals()[f"{animal}_only_downstream"]

    cond = {
        "name": animal.capitalize(),

        # ---- transformed values (for plotting) ----
        "up_trans": up_df["transformed_upstream"],
        "down_trans": down_df["transformed_downstream"],

        # ---- original values (for annotation only) ----
        "up_orig": up_df["upstream_len_in_scfr"],
        "down_orig": down_df["downstream_len_in_scfr"],
    }

    conditions.append(cond)

def plot_up_down_conditions_df(conditions, title="Upstream vs Downstream"):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    import pandas as pd

    # ---------------- Layout ----------------
    within_gap = 1.5      # gap between upstream & downstream of same animal
    group_gap = 4        # gap between animals
    box_width = 1.0
    side_padding = 1.5     # extra padding on left and right

    fig, ax = plt.subplots(figsize=(30, 10))  # SUPER big figure

    # ---------------- X positions ----------------
    x_positions = {}
    xticks = []
    xtick_labels = []

    x = side_padding
    for cond in conditions:
        name = cond["name"]
        up_x = x
        down_x = x + within_gap

        x_positions[(name, "Upstream")] = up_x
        x_positions[(name, "Downstream")] = down_x

        xticks.append((up_x + down_x) / 2)
        xtick_labels.append(name)

        x += group_gap

    # ---------------- Helper ----------------
    def get_whiskers(vals):
        q1 = vals.quantile(0.25)
        q3 = vals.quantile(0.75)
        iqr = q3 - q1
        lower = max(vals.min(), q1 - 1.5 * iqr)
        upper = min(vals.max(), q3 + 1.5 * iqr)
        return lower, upper, iqr

    # ---------------- Plot boxes + annotate ----------------
    for cond in conditions:
        name = cond["name"]
        data_map = {
            "Upstream": (cond.get("up_trans", pd.Series()).dropna(),
                         cond.get("up_orig", pd.Series()).dropna(),
                         "#7FB3D5"),
            "Downstream": (cond.get("down_trans", pd.Series()).dropna(),
                           cond.get("down_orig", pd.Series()).dropna(),
                           "#F5B7B1"),
        }

        for region, (vals_t, vals_o, color) in data_map.items():
            if len(vals_t) == 0:
                continue

            x_pos = x_positions[(name, region)]

            # Boxplot
            bp = ax.boxplot(
                vals_t,
                positions=[x_pos],
                widths=box_width,
                showfliers=False,
                patch_artist=True,
                boxprops=dict(facecolor=color, alpha=0.7),
                medianprops=dict(color="green", linewidth=2),
                whiskerprops=dict(color='black', linewidth=1.5),
                capprops=dict(color='red', linewidth=2)
            )

            # Stats
            median_t = vals_t.median()
            mean_t = vals_t.mean()
            wmin_t, wmax_t, iqr_t = get_whiskers(vals_t)

            median_o = vals_o.median() if len(vals_o) > 0 else 0
            mean_o = vals_o.mean() if len(vals_o) > 0 else 0
            min_o = vals_o.min() if len(vals_o) > 0 else 0
            max_o = vals_o.max() if len(vals_o) > 0 else 0

            offset = max(0.03 * iqr_t, 0.01 * (vals_t.max() - vals_t.min() + 1e-5))

            # Median annotation
            ax.text(x_pos, median_t + offset, f"{int(median_o)} bp",
                    ha="center", va="bottom", fontsize=12, fontweight="bold", color="green")

            # Mean line + dot + annotation
            ax.hlines(mean_t, x_pos - 0.5, x_pos + 0.5, colors="blue", linestyles="dotted", linewidth=2)
            ax.plot(x_pos, mean_t, "o", color="blue", zorder=5)
            ax.text(x_pos, mean_t - offset, f"{mean_o:.1f} bp", ha="center", va="top", fontsize=12, color="blue")

            # Min / Max annotation
            ax.text(x_pos, wmin_t - offset, f"{min_o:.0f} bp", ha="center", va="top", fontsize=10, color='red')
            ax.text(x_pos, wmax_t + offset, f"{max_o/1000:.2f} kb", ha="center", va="bottom", fontsize=10, color='red')

    # ---------------- Axis & Legend ----------------
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize=14, fontweight="bold")

    legend_elements = [
        Line2D([0], [0], color="#7FB3D5", lw=12, label="Upstream"),
        Line2D([0], [0], color="#F5B7B1", lw=12, label="Downstream"),
        Line2D([0], [0], color="green", lw=2, label="Median"),
        Line2D([0], [0], color="blue", lw=2, linestyle="dotted", label="Mean"),
        Line2D([0], [0], color='red', lw=2, label='Whisker caps (Min/Max)'),
    ]
    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    ax.set_ylabel("Transformed Counts (Yeo-Johnson Transform)", fontsize=14)
    ax.set_title(title, y=1.1, fontsize=20)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    # ---------------- Counts right above each animal ----------------
    y_max = ax.get_ylim()[1]
    y_count = y_max * 1.05  # slightly above plot top

    mapping = {"Upstream": "up_orig", "Downstream": "down_orig"}
    for cond in conditions:
        counts_text = []
        for region in ["Upstream", "Downstream"]:
            x_pos = x_positions[(cond["name"], region)]
            vals_orig = cond.get(mapping[region], pd.Series()).dropna()
            n = len(vals_orig)
            counts_text.append(f"{region[0]}: {n}")
        # place Up + Down counts in the same line
        ax.text((x_positions[(cond["name"], "Upstream")] + x_positions[(cond["name"], "Downstream")])/2,
                y_count,
                " | ".join(counts_text),
                ha="center", va="bottom", fontsize=14, fontweight="bold")

    # ---------------- Add horizontal padding ----------------
    all_x = list(x_positions.values())
    x_min, x_max = min(all_x), max(all_x)
    padding = side_padding
    ax.set_xlim(x_min - padding, x_max + padding)

    plt.savefig("../data/new_box_log_transform.png", dpi=300)
    plt.show()

plot_up_down_conditions_df(conditions)