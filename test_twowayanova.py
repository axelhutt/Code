#!/usr/bin/env python3

import numpy as np
import pandas as pd

#create data
df = pd.DataFrame({'water': np.repeat(['daily', 'weekly'], 15),
            'sun': np.tile(np.repeat(['low', 'med', 'high'], 5), 2),
            'height': [4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
                    6, 6, 6, 6, 6, 3, 3, 3, 3, 3,
                    4, 4, 4, 4, 4, 5, 5, 5, 5, 5]})
print(df)

import statsmodels.api as sm
from statsmodels.formula.api import ols

#perform two-way ANOVA
model = ols('height ~ C(water) + C(sun) + C(water):C(sun)', data=df).fit()
res = sm.stats.anova_lm(model, typ=2)
print(model)
print(res)

from statsmodels.stats.multicomp import pairwise_tukeyhsd
tukey = pairwise_tukeyhsd(endog=df['height'],groups=df['sun'], alpha=0.05)
print(tukey)