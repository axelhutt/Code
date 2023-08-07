#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# A few helper functions:
import utils

# To illustrate examples
import numpy as np
from scipy.stats import mannwhitneyu, normaltest

dataset = pd.read_csv('kickstarter_projects.csv')
print(dataset.head())

print(list(dataset.Category.unique()))
#print(dataset.Category)

tech = dataset.loc[(dataset.Category=='Technology'), :]
print(tech)
rfs = tech.loc[(tech.Subcategory.isin(("Robots", "Flight", "Sound"))), :]
print(rfs)

subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
subcat_order = ['Robots', 'Flight', 'Sound']

with sns.plotting_context("notebook", font_scale=1.4):
    # Create new plot, setting a logarithmic scale for y
    ax = utils.get_log_ax()
    
    # Plot with seaborn
    sns.boxplot(ax=ax, data=rfs, x='Subcategory', y='Goal', palette=subcat_palette[1:])
    
    # Label (adds axes labels and title), and show
    label_plot_for_subcats(ax)
    plt.savefig("plot1.png")
    