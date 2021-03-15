# Make Neut Curves and Analyze Mouse Neuts


```python
import math

import pandas as pd
import numpy as np
import seaborn as sns

import neutcurve

from plotnine import *

import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties

import scipy.stats
```


```python
CBP = neutcurve.colorschemes.CBPALETTE
CBM = neutcurve.colorschemes.CBMARKERS
theme_set(theme_seaborn(style='white', context='talk', font_scale=1))
```


```python
fidfs = {}
fidfs['rnd1'] = pd.read_csv('./fract_infect/201015_fractinfect.csv').drop(['Unnamed: 0'], axis=1)
fidfs['rnd2'] = pd.read_csv('./fract_infect/201016_fractinfect_rnd2.csv').drop(['Unnamed: 0'], axis=1)
fidfs['rnd3'] = pd.read_csv('./fract_infect/201017_fractinfect_rnd3.csv').drop(['Unnamed: 0'], axis=1)
fidfs['rnd4'] = pd.read_csv('./fract_infect/201020_fractinfect_rnd4.csv').drop(['Unnamed: 0'], axis=1)
fidfs['rnd5'] = pd.read_csv('./fract_infect/201021_fractinfect_rnd5.csv').drop(['Unnamed: 0'], axis=1)
fidfs['rnd6'] = pd.read_csv('./fract_infect/201022_fractinfect_rnd6.csv').drop(['Unnamed: 0'], axis=1)
```


```python
for fidf in fidfs.keys():
    display(fidfs[fidf].head())
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.988406</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.999059</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>1.013738</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.909787</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.900076</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.625716</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.710727</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.761430</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>1.159994</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.713713</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.565430</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.639826</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.635820</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.862117</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.814357</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>-0.002234</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>-0.002000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>-0.000372</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.021737</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000206</td>
      <td>0.135276</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.025000</td>
      <td>-0.001621</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.008333</td>
      <td>-0.001460</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.002778</td>
      <td>-0.001335</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000926</td>
      <td>-0.001225</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000309</td>
      <td>-0.000772</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.008333</td>
      <td>0.000013</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.002778</td>
      <td>0.006947</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000926</td>
      <td>0.030349</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000309</td>
      <td>0.209554</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000103</td>
      <td>0.642440</td>
    </tr>
  </tbody>
</table>
</div>



```python
smpls = {}
smpls['rnd1'] = pd.read_csv('./sample_maps/rnd1_sample_map.csv').rename(columns={'Sample': 'serum'})
smpls['rnd2'] = pd.read_csv('./sample_maps/rnd2_sample_map.csv').rename(columns={'Sample': 'serum'})
smpls['rnd3'] = pd.read_csv('./sample_maps/rnd3_sample_map.csv').rename(columns={'Sample': 'serum'})
smpls['rnd4'] = pd.read_csv('./sample_maps/rnd4_sample_map.csv').rename(columns={'Sample': 'serum'})
smpls['rnd5'] = pd.read_csv('./sample_maps/rnd5_sample_map.csv').rename(columns={'Sample': 'serum'})
smpls['rnd6'] = pd.read_csv('./sample_maps/rnd6_sample_map.csv').rename(columns={'Sample': 'serum'})
```


```python
for smpl in smpls.keys():
    display(smpls[smpl].head())
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201013</td>
      <td>0793-1</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201013</td>
      <td>0846-1</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201013</td>
      <td>0807-1</td>
      <td>Plate1</td>
      <td>3</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201013</td>
      <td>0813-1</td>
      <td>Plate1</td>
      <td>4</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201013</td>
      <td>0863-1</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201014</td>
      <td>0805-1</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201014</td>
      <td>0854-1</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201014</td>
      <td>0858-1</td>
      <td>Plate1</td>
      <td>3</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201014</td>
      <td>0817-1</td>
      <td>Plate1</td>
      <td>4</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201014</td>
      <td>0834-1</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201015</td>
      <td>0839-2</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201015</td>
      <td>0846-2</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201015</td>
      <td>0852-2</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201015</td>
      <td>0859-2</td>
      <td>Plate2</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201015</td>
      <td>0863-2</td>
      <td>Plate3</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.05</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201018</td>
      <td>0840-2</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.016667</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201018</td>
      <td>0850-2</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.016667</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201018</td>
      <td>0853-2</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.016667</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201018</td>
      <td>0856-2</td>
      <td>Plate2</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.016667</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201018</td>
      <td>0864-2</td>
      <td>Plate3</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.016667</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201019</td>
      <td>0795-2</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.025</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201019</td>
      <td>0803-2</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.025</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201019</td>
      <td>0807-2</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.025</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201019</td>
      <td>0812-2</td>
      <td>Plate2</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.025</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201019</td>
      <td>0818-2</td>
      <td>Plate3</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout2.csv</td>
      <td>0.025</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DateSetUp</th>
      <th>serum</th>
      <th>Plate</th>
      <th>SampleNum</th>
      <th>Virus</th>
      <th>PlateLayout</th>
      <th>StartDil</th>
      <th>DilFactor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>201020</td>
      <td>0854-1</td>
      <td>Plate1</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.050000</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>201020</td>
      <td>0856-1</td>
      <td>Plate1</td>
      <td>2</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.050000</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>201020</td>
      <td>0805-1</td>
      <td>Plate1</td>
      <td>3</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.050000</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>201020</td>
      <td>0818-1</td>
      <td>Plate1</td>
      <td>4</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.008333</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>201020</td>
      <td>0817-1</td>
      <td>Plate2</td>
      <td>1</td>
      <td>S-d21-D614G</td>
      <td>layout1.csv</td>
      <td>0.008333</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>


#### Add plate number to fraction infectivity df


```python
dfs={}
for rnd in fidfs.keys():
    dfs[rnd] = fidfs[rnd].merge(smpls[rnd], how='outer', on='serum').drop(['Virus', 'PlateLayout', 'StartDil', 'DilFactor'], axis=1)
```


```python
for rnd in dfs.keys():
    display(dfs[rnd].head())
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.988406</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.999059</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>1.013738</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.909787</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0840-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.900076</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.625716</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.710727</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.761430</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>1.159994</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.713713</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.565430</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.639826</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.635820</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.862117</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.814357</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>-0.002234</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>-0.002000</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>-0.000372</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.021737</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000206</td>
      <td>0.135276</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.025000</td>
      <td>-0.001621</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.008333</td>
      <td>-0.001460</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.002778</td>
      <td>-0.001335</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000926</td>
      <td>-0.001225</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000309</td>
      <td>-0.000772</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.008333</td>
      <td>0.000013</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.002778</td>
      <td>0.006947</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000926</td>
      <td>0.030349</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000309</td>
      <td>0.209554</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000103</td>
      <td>0.642440</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>


### Fix Plate 4 Issue for Rnd1

I forgot to discard the last 30 uL from the serum dilutions in row H of plate 4 in Rnd1, so I need to drop the final dilution for all plate 


```python
dfs['rnd1'].drop(dfs['rnd1'][(dfs['rnd1']['Plate']=='Plate4')&(dfs['rnd1']['concentration']==(0.05/(3**6)))].index, axis=0, inplace=True)
dfs['rnd1'] = dfs['rnd1'].sort_values(['Plate', 'SampleNum', 'replicate'])
display(dfs['rnd1'])
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>336</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.245596</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>337</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.396659</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>338</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.479400</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>339</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.822955</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>340</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>1.009480</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>51</th>
      <td>20/130 Ref</td>
      <td>S-d21-D614G</td>
      <td>2</td>
      <td>0.002778</td>
      <td>0.067901</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>5</td>
    </tr>
    <tr>
      <th>52</th>
      <td>20/130 Ref</td>
      <td>S-d21-D614G</td>
      <td>2</td>
      <td>0.000926</td>
      <td>0.184104</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>5</td>
    </tr>
    <tr>
      <th>53</th>
      <td>20/130 Ref</td>
      <td>S-d21-D614G</td>
      <td>2</td>
      <td>0.000309</td>
      <td>0.379211</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>5</td>
    </tr>
    <tr>
      <th>54</th>
      <td>20/130 Ref</td>
      <td>S-d21-D614G</td>
      <td>2</td>
      <td>0.000103</td>
      <td>0.605273</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>5</td>
    </tr>
    <tr>
      <th>55</th>
      <td>20/130 Ref</td>
      <td>S-d21-D614G</td>
      <td>2</td>
      <td>0.000034</td>
      <td>0.657573</td>
      <td>201013</td>
      <td>Plate7</td>
      <td>5</td>
    </tr>
  </tbody>
</table>
<p>384 rows × 8 columns</p>
</div>


## Calculate Fits using Neut Curve

For some reason must drop serum pool from rnd2.


```python
display(dfs['rnd2'][dfs['rnd2']['serum'] == '2017-2018 serum pool'])
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>322</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.050000</td>
      <td>0.860817</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>323</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.016667</td>
      <td>0.674258</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>324</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.005556</td>
      <td>0.980501</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>325</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.001852</td>
      <td>0.887813</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>326</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000617</td>
      <td>0.873441</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>327</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000206</td>
      <td>0.786647</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>328</th>
      <td>2017-2018 serum pool</td>
      <td>S-d21-D614G</td>
      <td>1</td>
      <td>0.000069</td>
      <td>0.830204</td>
      <td>201014</td>
      <td>Plate6</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



```python
dfs['rnd2'] = dfs['rnd2'][dfs['rnd2']['serum'] != '2017-2018 serum pool']
```


```python
fits = {}
fit_params = {}

for rnd in dfs.keys():
    fits[rnd] = neutcurve.CurveFits(dfs[rnd], fixtop=True)
    fit_params[rnd] = fits[rnd].fitParams(ics=[50, 90])
    fit_params[rnd]['nt50'] = 1/fit_params[rnd]['ic50']
    fit_params[rnd]['nt90'] = 1/fit_params[rnd]['ic90']
```

    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/neutcurve/hillcurve.py:689: RuntimeWarning: invalid value encountered in power
    /home/kdusenbu/.local/lib/python3.8/site-packages/scipy/optimize/minpack.py:828: OptimizeWarning: Covariance of the parameters could not be estimated



```python
for rnd in fit_params.keys():
    display(fit_params[rnd].head())

```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.010629</td>
      <td>interpolated</td>
      <td>0.0106</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.010629</td>
      <td>0.831619</td>
      <td>True</td>
      <td>0</td>
      <td>94.080392</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0846-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>4.132835</td>
      <td>0.265601</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0807-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.020702</td>
      <td>interpolated</td>
      <td>0.0207</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.020702</td>
      <td>0.903648</td>
      <td>True</td>
      <td>0</td>
      <td>48.303432</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0813-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>21279.859733</td>
      <td>0.126899</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0863-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000207</td>
      <td>interpolated</td>
      <td>0.000207</td>
      <td>0.000825</td>
      <td>interpolated</td>
      <td>0.000825</td>
      <td>0.000207</td>
      <td>1.588156</td>
      <td>True</td>
      <td>0</td>
      <td>4834.077623</td>
      <td>1211.892515</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>317328.751658</td>
      <td>0.079256</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0854-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>3.622431</td>
      <td>3.996576</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0858-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.002211</td>
      <td>interpolated</td>
      <td>0.00221</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.002211</td>
      <td>0.528155</td>
      <td>True</td>
      <td>0</td>
      <td>452.261297</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0817-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000075</td>
      <td>interpolated</td>
      <td>7.54e-05</td>
      <td>0.000418</td>
      <td>interpolated</td>
      <td>0.000418</td>
      <td>0.000075</td>
      <td>1.283698</td>
      <td>True</td>
      <td>0</td>
      <td>13258.708776</td>
      <td>2394.126118</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0834-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000096</td>
      <td>interpolated</td>
      <td>9.6e-05</td>
      <td>0.000757</td>
      <td>interpolated</td>
      <td>0.000757</td>
      <td>0.000096</td>
      <td>1.063690</td>
      <td>True</td>
      <td>0</td>
      <td>10419.562318</td>
      <td>1320.514686</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.324660</td>
      <td>0.264365</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2017-2018 Serum Pool</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.290436</td>
      <td>0.263139</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Ty1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.046310</td>
      <td>interpolated</td>
      <td>0.0463</td>
      <td>0.854942</td>
      <td>interpolated</td>
      <td>0.855</td>
      <td>0.046310</td>
      <td>0.753589</td>
      <td>True</td>
      <td>0</td>
      <td>21.593740</td>
      <td>1.169670</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Ref 20/130</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000249</td>
      <td>interpolated</td>
      <td>0.000249</td>
      <td>0.002299</td>
      <td>interpolated</td>
      <td>0.0023</td>
      <td>0.000249</td>
      <td>0.988148</td>
      <td>True</td>
      <td>0</td>
      <td>4019.844087</td>
      <td>435.031720</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0873-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000017</td>
      <td>interpolated</td>
      <td>1.66e-05</td>
      <td>0.000079</td>
      <td>interpolated</td>
      <td>7.86e-05</td>
      <td>0.000017</td>
      <td>1.415264</td>
      <td>True</td>
      <td>0</td>
      <td>60100.762592</td>
      <td>12724.170768</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000029</td>
      <td>interpolated</td>
      <td>2.86e-05</td>
      <td>0.000309</td>
      <td>interpolated</td>
      <td>0.000309</td>
      <td>0.000029</td>
      <td>0.922937</td>
      <td>True</td>
      <td>0</td>
      <td>34985.127556</td>
      <td>3235.667865</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0850-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000253</td>
      <td>interpolated</td>
      <td>0.000253</td>
      <td>0.000727</td>
      <td>interpolated</td>
      <td>0.000727</td>
      <td>0.000253</td>
      <td>2.084434</td>
      <td>True</td>
      <td>0</td>
      <td>3949.368926</td>
      <td>1376.363387</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0853-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.001809</td>
      <td>interpolated</td>
      <td>0.00181</td>
      <td>0.008289</td>
      <td>interpolated</td>
      <td>0.00829</td>
      <td>0.001809</td>
      <td>1.443564</td>
      <td>True</td>
      <td>0</td>
      <td>552.742140</td>
      <td>120.639716</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0856-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000318</td>
      <td>interpolated</td>
      <td>0.000318</td>
      <td>0.011380</td>
      <td>interpolated</td>
      <td>0.0114</td>
      <td>0.000318</td>
      <td>0.614245</td>
      <td>True</td>
      <td>0</td>
      <td>3143.337332</td>
      <td>87.875780</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0864-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.51e-06</td>
      <td>0.000038</td>
      <td>interpolated</td>
      <td>3.78e-05</td>
      <td>0.000009</td>
      <td>1.474670</td>
      <td>True</td>
      <td>0</td>
      <td>117451.650829</td>
      <td>26471.006109</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000017</td>
      <td>interpolated</td>
      <td>1.69e-05</td>
      <td>0.000080</td>
      <td>interpolated</td>
      <td>7.99e-05</td>
      <td>0.000017</td>
      <td>1.416749</td>
      <td>True</td>
      <td>0</td>
      <td>59027.745711</td>
      <td>12517.348896</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0824-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000006</td>
      <td>interpolated</td>
      <td>6.06e-06</td>
      <td>0.000027</td>
      <td>interpolated</td>
      <td>2.75e-05</td>
      <td>0.000006</td>
      <td>1.452346</td>
      <td>True</td>
      <td>0</td>
      <td>165125.537880</td>
      <td>36372.975212</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0818-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000007</td>
      <td>interpolated</td>
      <td>7.26e-06</td>
      <td>0.000039</td>
      <td>interpolated</td>
      <td>3.91e-05</td>
      <td>0.000007</td>
      <td>1.304890</td>
      <td>True</td>
      <td>0</td>
      <td>137831.787498</td>
      <td>25589.832347</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0831-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.73e-06</td>
      <td>0.000042</td>
      <td>interpolated</td>
      <td>4.23e-05</td>
      <td>0.000009</td>
      <td>1.392561</td>
      <td>True</td>
      <td>0</td>
      <td>114557.690454</td>
      <td>23647.284318</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0816-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.88e-06</td>
      <td>0.000038</td>
      <td>interpolated</td>
      <td>3.76e-05</td>
      <td>0.000009</td>
      <td>1.523549</td>
      <td>True</td>
      <td>0</td>
      <td>112644.478616</td>
      <td>26630.606659</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000148</td>
      <td>interpolated</td>
      <td>0.000148</td>
      <td>0.000527</td>
      <td>interpolated</td>
      <td>0.000527</td>
      <td>0.000148</td>
      <td>1.733238</td>
      <td>True</td>
      <td>0</td>
      <td>6740.868115</td>
      <td>1897.409580</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0834-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000177</td>
      <td>interpolated</td>
      <td>0.000177</td>
      <td>0.000807</td>
      <td>interpolated</td>
      <td>0.000807</td>
      <td>0.000177</td>
      <td>1.447257</td>
      <td>True</td>
      <td>0</td>
      <td>5657.925665</td>
      <td>1239.685829</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2017-2018 Serum Pool</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.077750</td>
      <td>2.549505</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Ty1-Fc</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050143</td>
      <td>interpolated</td>
      <td>0.0501</td>
      <td>0.337347</td>
      <td>interpolated</td>
      <td>0.337</td>
      <td>0.050143</td>
      <td>1.152656</td>
      <td>True</td>
      <td>0</td>
      <td>19.942813</td>
      <td>2.964304</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0817-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000126</td>
      <td>interpolated</td>
      <td>0.000126</td>
      <td>0.000506</td>
      <td>interpolated</td>
      <td>0.000506</td>
      <td>0.000126</td>
      <td>1.578931</td>
      <td>True</td>
      <td>0</td>
      <td>7943.119901</td>
      <td>1975.292048</td>
    </tr>
  </tbody>
</table>
</div>


### Merge fit parameter data and sample data


```python
fits_dfs = {}

for rnd in fit_params.keys():
    fits_dfs[rnd] = pd.merge(fit_params[rnd], smpls[rnd], how='outer', on='serum')
    fits_dfs[rnd].drop(['PlateLayout', 'Virus', 'DilFactor', 'StartDil'], axis=1, inplace=True)
```


```python
for rnd in fits_dfs.keys():
    display(fits_dfs[rnd].head())
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0793-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.010629</td>
      <td>interpolated</td>
      <td>0.0106</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.010629</td>
      <td>0.831619</td>
      <td>True</td>
      <td>0</td>
      <td>94.080392</td>
      <td>20.000000</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0846-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>4.132835</td>
      <td>0.265601</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0807-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.020702</td>
      <td>interpolated</td>
      <td>0.0207</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.020702</td>
      <td>0.903648</td>
      <td>True</td>
      <td>0</td>
      <td>48.303432</td>
      <td>20.000000</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0813-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>21279.859733</td>
      <td>0.126899</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201013</td>
      <td>Plate1</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0863-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000207</td>
      <td>interpolated</td>
      <td>0.000207</td>
      <td>0.000825</td>
      <td>interpolated</td>
      <td>0.000825</td>
      <td>0.000207</td>
      <td>1.588156</td>
      <td>True</td>
      <td>0</td>
      <td>4834.077623</td>
      <td>1211.892515</td>
      <td>201013</td>
      <td>Plate2</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0805-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>317328.751658</td>
      <td>0.079256</td>
      <td>True</td>
      <td>0.0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0854-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>3.622431</td>
      <td>3.996576</td>
      <td>True</td>
      <td>0.0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0858-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.002211</td>
      <td>interpolated</td>
      <td>0.00221</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.002211</td>
      <td>0.528155</td>
      <td>True</td>
      <td>0.0</td>
      <td>452.261297</td>
      <td>20.000000</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0817-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000075</td>
      <td>interpolated</td>
      <td>7.54e-05</td>
      <td>0.000418</td>
      <td>interpolated</td>
      <td>0.000418</td>
      <td>0.000075</td>
      <td>1.283698</td>
      <td>True</td>
      <td>0.0</td>
      <td>13258.708776</td>
      <td>2394.126118</td>
      <td>201014</td>
      <td>Plate1</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0834-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000096</td>
      <td>interpolated</td>
      <td>9.6e-05</td>
      <td>0.000757</td>
      <td>interpolated</td>
      <td>0.000757</td>
      <td>0.000096</td>
      <td>1.063690</td>
      <td>True</td>
      <td>0.0</td>
      <td>10419.562318</td>
      <td>1320.514686</td>
      <td>201014</td>
      <td>Plate2</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GF-8</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.324660</td>
      <td>0.264365</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2017-2018 Serum Pool</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.290436</td>
      <td>0.263139</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Ty1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.046310</td>
      <td>interpolated</td>
      <td>0.0463</td>
      <td>0.854942</td>
      <td>interpolated</td>
      <td>0.855</td>
      <td>0.046310</td>
      <td>0.753589</td>
      <td>True</td>
      <td>0</td>
      <td>21.593740</td>
      <td>1.169670</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Ref 20/130</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000249</td>
      <td>interpolated</td>
      <td>0.000249</td>
      <td>0.002299</td>
      <td>interpolated</td>
      <td>0.0023</td>
      <td>0.000249</td>
      <td>0.988148</td>
      <td>True</td>
      <td>0</td>
      <td>4019.844087</td>
      <td>435.031720</td>
      <td>201015</td>
      <td>Plate9</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0873-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000017</td>
      <td>interpolated</td>
      <td>1.66e-05</td>
      <td>0.000079</td>
      <td>interpolated</td>
      <td>7.86e-05</td>
      <td>0.000017</td>
      <td>1.415264</td>
      <td>True</td>
      <td>0</td>
      <td>60100.762592</td>
      <td>12724.170768</td>
      <td>201015</td>
      <td>Plate8</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0840-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000029</td>
      <td>interpolated</td>
      <td>2.86e-05</td>
      <td>0.000309</td>
      <td>interpolated</td>
      <td>0.000309</td>
      <td>0.000029</td>
      <td>0.922937</td>
      <td>True</td>
      <td>0</td>
      <td>34985.127556</td>
      <td>3235.667865</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0850-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000253</td>
      <td>interpolated</td>
      <td>0.000253</td>
      <td>0.000727</td>
      <td>interpolated</td>
      <td>0.000727</td>
      <td>0.000253</td>
      <td>2.084434</td>
      <td>True</td>
      <td>0</td>
      <td>3949.368926</td>
      <td>1376.363387</td>
      <td>201018</td>
      <td>Plate1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0853-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.001809</td>
      <td>interpolated</td>
      <td>0.00181</td>
      <td>0.008289</td>
      <td>interpolated</td>
      <td>0.00829</td>
      <td>0.001809</td>
      <td>1.443564</td>
      <td>True</td>
      <td>0</td>
      <td>552.742140</td>
      <td>120.639716</td>
      <td>201018</td>
      <td>Plate2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0856-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000318</td>
      <td>interpolated</td>
      <td>0.000318</td>
      <td>0.011380</td>
      <td>interpolated</td>
      <td>0.0114</td>
      <td>0.000318</td>
      <td>0.614245</td>
      <td>True</td>
      <td>0</td>
      <td>3143.337332</td>
      <td>87.875780</td>
      <td>201018</td>
      <td>Plate2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0864-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.51e-06</td>
      <td>0.000038</td>
      <td>interpolated</td>
      <td>3.78e-05</td>
      <td>0.000009</td>
      <td>1.474670</td>
      <td>True</td>
      <td>0</td>
      <td>117451.650829</td>
      <td>26471.006109</td>
      <td>201018</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0828-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000017</td>
      <td>interpolated</td>
      <td>1.69e-05</td>
      <td>0.000080</td>
      <td>interpolated</td>
      <td>7.99e-05</td>
      <td>0.000017</td>
      <td>1.416749</td>
      <td>True</td>
      <td>0</td>
      <td>59027.745711</td>
      <td>12517.348896</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0824-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000006</td>
      <td>interpolated</td>
      <td>6.06e-06</td>
      <td>0.000027</td>
      <td>interpolated</td>
      <td>2.75e-05</td>
      <td>0.000006</td>
      <td>1.452346</td>
      <td>True</td>
      <td>0</td>
      <td>165125.537880</td>
      <td>36372.975212</td>
      <td>201019</td>
      <td>Plate4</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0818-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000007</td>
      <td>interpolated</td>
      <td>7.26e-06</td>
      <td>0.000039</td>
      <td>interpolated</td>
      <td>3.91e-05</td>
      <td>0.000007</td>
      <td>1.304890</td>
      <td>True</td>
      <td>0</td>
      <td>137831.787498</td>
      <td>25589.832347</td>
      <td>201019</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0831-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.73e-06</td>
      <td>0.000042</td>
      <td>interpolated</td>
      <td>4.23e-05</td>
      <td>0.000009</td>
      <td>1.392561</td>
      <td>True</td>
      <td>0</td>
      <td>114557.690454</td>
      <td>23647.284318</td>
      <td>201019</td>
      <td>Plate3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0816-2</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000009</td>
      <td>interpolated</td>
      <td>8.88e-06</td>
      <td>0.000038</td>
      <td>interpolated</td>
      <td>3.76e-05</td>
      <td>0.000009</td>
      <td>1.523549</td>
      <td>True</td>
      <td>0</td>
      <td>112644.478616</td>
      <td>26630.606659</td>
      <td>201019</td>
      <td>Plate7</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0824-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000148</td>
      <td>interpolated</td>
      <td>0.000148</td>
      <td>0.000527</td>
      <td>interpolated</td>
      <td>0.000527</td>
      <td>0.000148</td>
      <td>1.733238</td>
      <td>True</td>
      <td>0</td>
      <td>6740.868115</td>
      <td>1897.409580</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0834-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000177</td>
      <td>interpolated</td>
      <td>0.000177</td>
      <td>0.000807</td>
      <td>interpolated</td>
      <td>0.000807</td>
      <td>0.000177</td>
      <td>1.447257</td>
      <td>True</td>
      <td>0</td>
      <td>5657.925665</td>
      <td>1239.685829</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2017-2018 Serum Pool</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.050000</td>
      <td>lower</td>
      <td>&gt;0.05</td>
      <td>0.077750</td>
      <td>2.549505</td>
      <td>True</td>
      <td>0</td>
      <td>20.000000</td>
      <td>20.000000</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Ty1-Fc</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.050143</td>
      <td>interpolated</td>
      <td>0.0501</td>
      <td>0.337347</td>
      <td>interpolated</td>
      <td>0.337</td>
      <td>0.050143</td>
      <td>1.152656</td>
      <td>True</td>
      <td>0</td>
      <td>19.942813</td>
      <td>2.964304</td>
      <td>201020</td>
      <td>Plate3</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0817-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000126</td>
      <td>interpolated</td>
      <td>0.000126</td>
      <td>0.000506</td>
      <td>interpolated</td>
      <td>0.000506</td>
      <td>0.000126</td>
      <td>1.578931</td>
      <td>True</td>
      <td>0</td>
      <td>7943.119901</td>
      <td>1975.292048</td>
      <td>201020</td>
      <td>Plate2</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>


### Find any samples with interpolated IC50s

Note that these samples were re-run in later rounds of neuts


```python
for rnd in fits_dfs.keys():
    print(rnd)
    display(fits_dfs[rnd][fits_dfs[rnd]['ic50_bound']=='upper'])
```

    rnd1



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>14</th>
      <td>0829-1</td>
      <td>S-d21-D614G</td>
      <td>average</td>
      <td>2</td>
      <td>0.000206</td>
      <td>upper</td>
      <td>&lt;0.000206</td>
      <td>0.000451</td>
      <td>interpolated</td>
      <td>0.000451</td>
      <td>0.000179</td>
      <td>2.376361</td>
      <td>True</td>
      <td>0</td>
      <td>4860.0</td>
      <td>2218.104626</td>
      <td>201013</td>
      <td>Plate4</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>


    rnd2



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


    rnd3



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


    rnd4



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


    rnd5



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


    rnd6



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>replicate</th>
      <th>nreplicates</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>ic50_str</th>
      <th>ic90</th>
      <th>ic90_bound</th>
      <th>ic90_str</th>
      <th>midpoint</th>
      <th>slope</th>
      <th>top</th>
      <th>bottom</th>
      <th>nt50</th>
      <th>nt90</th>
      <th>DateSetUp</th>
      <th>Plate</th>
      <th>SampleNum</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


## Plot all neut curves


```python
for rnd in fits.keys():
    fig, axes = fits[rnd].plotSera(xlabel='concentration/dilution')
```

## Plot Ty1-Fc control neut curves


```python
curve = fits['rnd1'].getCurve(serum='Ty1-Fc', virus='S-d21-D614G', replicate='1')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd2'].getCurve(serum='Ty1-FC', virus='S-d21-D614G', replicate='1')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd3'].getCurve(serum='Ty1', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd4'].getCurve(serum='Ty1', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd5'].getCurve(serum='Ty1', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd6'].getCurve(serum='Ty1-Fc', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd1'].getCurve(serum='20/130 Ref', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {1/curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd3'].getCurve(serum='Ref 20/130', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {1/curve.ic50():.3g}")
fig, ax = curve.plot()
```


```python
curve = fits['rnd5'].getCurve(serum='Ref 20/130', virus='S-d21-D614G', replicate='average')
print(f"The IC50 is {1/curve.ic50():.3g}")
fig, ax = curve.plot()
```

## Combine data into one df



```python
all_data_df = pd.concat(fits_dfs)
```

## Drop Earlier Data for Samples that were re-run

Drop first sample and only keep re-run for samples that were re-run.


```python
rerun = ['0829-1', '0879-1', '0858-1', '0854-1', '0856-1', '0805-1', '0818-1',
         '0817-1', '0835-1', '0873-1', '0866-1', '0824-1', '0834-1']
```


```python
print(len(all_data_df))
drop_idxs = []
for sample in rerun:
    sample_data = all_data_df[all_data_df['serum']==sample]
    drop_idxs.append(sample_data['DateSetUp'].idxmin()) # drop earlier run

cleaned_data_df = all_data_df.drop(drop_idxs).reset_index(drop=True)
display(cleaned_data_df.head())
print(len(cleaned_data_df))
```

### Output relevant data to csv

Output columns: `serum`, `ic50`, `ic90`, `nt50`, `nt90`, `DateSetUp`, `Plate`, `SampleNum`


```python
export_df_all = all_data_df[['serum', 'ic50', 'ic90', 'nt50', 'nt90', 'DateSetUp', 'Plate', 'SampleNum']]
export_df_clean = cleaned_data_df[['serum', 'ic50', 'ic90', 'nt50', 'nt90', 'DateSetUp', 'Plate', 'SampleNum']]

# `all_neut_results` includes all neuts, including the first neuts
# from samples I reran and all controls.
export_df_all.to_csv('./all_neut_results.csv') 
# `mouse_plus_ctrls_neuts` only includes the later neuts for the
# samples I re-ran and all controls.
export_df_clean.to_csv('./mouse_plus_ctrls_neuts.csv')
```

## Initial Analyses

Only including re-runs for samples I ran twice.
Not including controls.

Export csv after adding timepoint and mouse data.


```python
drop_samples = ['GF-8', 'Ref 20/130', '2017-2018 Serum Pool', '2017-2018 Pool', '20/130 Ref', '2017-2018 serum pool', 'Ty1', 'Ty1-FC', 'Ty1-Fc']
cleaned_data_df = cleaned_data_df[~cleaned_data_df['serum'].isin(drop_samples)].copy()
cleaned_data_df['Timepoint'] = cleaned_data_df['serum'].apply(lambda x: 'Prime' if '-1' in x else 'Boost')
cleaned_data_df['Mouse'] = cleaned_data_df['serum'].apply(lambda x: x[:-2])
```


```python
display(cleaned_data_df.sort_values('Mouse').reset_index(drop=True).head())
```


```python
groups = {'Group 1': ['0840', '0837', '0839', '0794', '0795', '0793'], 
          'Group 2': ['0801', '0803', '0805', '0850', '0848', '0846'],
          'Group 3': ['0808', '0807', '0810', '0853', '0854', '0852'],
          'Group 4': ['0813', '0812', '0815', '0856', '0858', '0859'],
          'Group 5': ['0816', '0818', '0817', '0864', '0861', '0863'],
          'Group 6': ['0835', '0831', '0834', '0867', '0866', '0870'],
          'Group 7': ['0827', '0828', '0829', '0874', '0873', '0872'],
          'Group 8': ['0823', '0824', '0825', '0879', '0877', '0876']}
```


```python
cleaned_data_df["Group"] = cleaned_data_df["Mouse"].apply(lambda x: [group for group in groups.keys() if x in groups[group]][0])
```


```python
display(cleaned_data_df.head())
```

### Export csv of cleaned mouse neut data

Only includes latest run for samples I re-ran and doesn't include human naive serum or Ty1-FC controls.


```python
export_df_clean_noctrls = cleaned_data_df[['serum', 'ic50', 'ic90', 'nt50', 'nt90', 'DateSetUp', 'Plate', 'SampleNum', 'Timepoint', 'Mouse', 'Group']]
# `mouse_neuts` only includes the mouse neuts (no controls) and 
# only includes the later data for the samples I re-ran
export_df_clean_noctrls.to_csv('./mouse_neuts.csv')
```

### Initial plotting


```python
group_list = ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5',
              'Group 6', 'Group 7', 'Group 8']
group_cat = pd.Categorical(cleaned_data_df['Group'], categories=group_list)

tp_list = ['Prime', 'Boost']
tp_cat = pd.Categorical(cleaned_data_df['Timepoint'], categories=tp_list)

# assign to a new column in the DataFrame
cleaned_data_df = cleaned_data_df.assign(group_order = group_cat)
cleaned_data_df = cleaned_data_df.assign(timepoint_order = tp_cat)
```


```python
ic50s_plot = (ggplot(cleaned_data_df, aes('group_order', 'nt50', color='tp_cat')) +
              geom_boxplot(outlier_alpha=0) +
              geom_point(size=2, alpha=0.5, position=position_dodge(width=0.75)) +
              scale_color_manual(values=CBP) +
              theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    figure_size=(8, 6)) +
#               facet_wrap('~ timepoint_order') +
              geom_hline(yintercept=20, color='grey', linetype='dashed') +
              scale_y_continuous(trans='log10') +
              labs(color='Timepoint')
             )

_ = ic50s_plot.draw()
```

### Convert to Markdown


```python
!jupyter nbconvert MouseNeuts.ipynb --to markdown
```


```python

```


```python

```
