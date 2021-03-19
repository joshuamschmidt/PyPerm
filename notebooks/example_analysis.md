### PySetPerm design
The pysetperm.py module includes a number of classes that provide simple building blocks for testing set enrichments.
Features can be anything: genes, regulatory elements etc. as long as they have chr, start (1-based!), end(1-based) and name columns:


```bash
%%bash
head -n3 data/genes.txt
```

    chr	start	end	gene
    1	904115	905037	HES4
    1	921857	922761	ISG15


Annotations are also simply specified:


```bash
%%bash
head -n3 data/kegg.txt
```

    id	feature	name
    hsa00010	ACSS1	Glycolysis / Gluconeogenesis
    hsa00010	ACSS2	Glycolysis / Gluconeogenesis


### An example analysis
We import features and annotaions via respective classes. Features can be altered with a distance (i.e. genes +- 2000 bp). Annotations can also be filtered to have a minimum set size (i.e. at least 5 genes)


```python
import pysetperm as psp
features = psp.Features('data/genes.txt', 2000)
annotations = psp.AnnotationSet('data/kegg.txt', features.features_user_def, 5)
n_perms = 200000
cores = 10
```

Initiate test groups using the Input class:


```python
e_input = psp.Input('data/eastern_candidates.txt',
                    'data/eastern_background.txt.gz',
                    features,
                    annotations)

c_input = psp.Input('data/central_candidates.txt',
                    'data/central_background.txt.gz',
                    features,
                    annotations)

i_input = psp.Input('data/internal_candidates.txt',
                    'data/internal_background.txt.gz',
                    features,
                    annotations)


```

A Permutation class holds the permuted datasets.


```python
e_permutations = psp.Permutation(e_input, n_perms, cores)
c_permutations = psp.Permutation(c_input, n_perms, cores)
i_permutations = psp.Permutation(i_input, n_perms, cores)
```

Once permutions are completed, we determine the distribution of the Pr. X of genes belonging to Set1...n, using the SetPerPerm class. This structure enables the easy generation of joint distributions.


```python
e_per_set = psp.SetPerPerm(e_permutations,
                           annotations,
                           e_input,
                           cores)
c_per_set = psp.SetPerPerm(c_permutations,
                           annotations,
                           c_input,
                           cores)
i_per_set = psp.SetPerPerm(i_permutations,
                           annotations,
                           i_input,
                           cores)
```

Here, we can use join_objects() methods for both Imput and SetPerPerm objects, to get the joint distribution of two or more indpendent tests.


```python
# combine sims
ec_input = psp.Input.join_objects(e_input, c_input)
ec_per_set = psp.SetPerPerm.join_objects(e_per_set, c_per_set)
ei_input = psp.Input.join_objects(e_input, i_input)
ei_per_set = psp.SetPerPerm.join_objects(e_per_set, i_per_set)
ci_input = psp.Input.join_objects(c_input, i_input)
ci_per_set = psp.SetPerPerm.join_objects(c_per_set, i_per_set)
eci_input = psp.Input.join_objects(ec_input, i_input)
eci_per_set = psp.SetPerPerm.join_objects(ec_per_set, i_per_set)
```

Call the make_results_table function to generate a pandas format results table.


```python
# results
e_results = psp.make_results_table(e_input, annotations, e_per_set)
c_results = psp.make_results_table(c_input, annotations, c_per_set)
i_results = psp.make_results_table(i_input, annotations, i_per_set)
ec_results = psp.make_results_table(ec_input, annotations, ec_per_set)
ei_results = psp.make_results_table(ei_input, annotations, ei_per_set)
ci_results = psp.make_results_table(ci_input, annotations, ci_per_set)
eci_results = psp.make_results_table(eci_input, annotations, eci_per_set)
```


```python
from itables import show
from IPython.display import display
from ipywidgets import HBox, VBox
import ipywidgets as widgets
display(e_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>226</th>
      <td>hsa04658</td>
      <td>Th1 and Th2 cell differentiation</td>
      <td>[CD3D, CD3G, IL12RB1, IL13, IL4, MAML3, MAPK14...</td>
      <td>11</td>
      <td>3.698510</td>
      <td>0.000495</td>
      <td>0.999880</td>
      <td>0.092385</td>
      <td>1.0</td>
      <td>0.133954</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>4</th>
      <td>hsa00051</td>
      <td>Fructose and mannose metabolism</td>
      <td>[FUK, GMDS, HKDC1, MPI, PMM1, SORD]</td>
      <td>6</td>
      <td>1.466265</td>
      <td>0.000820</td>
      <td>0.999925</td>
      <td>0.092385</td>
      <td>1.0</td>
      <td>0.133954</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>47</th>
      <td>hsa00520</td>
      <td>Amino sugar and nucleotide sugar metabolism</td>
      <td>[CYB5RL, FUK, GFPT2, GMDS, HKDC1, MPI, PMM1]</td>
      <td>7</td>
      <td>1.929355</td>
      <td>0.001095</td>
      <td>0.999865</td>
      <td>0.092385</td>
      <td>1.0</td>
      <td>0.133954</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>334</th>
      <td>hsa05169</td>
      <td>Epstein-Barr virus infection</td>
      <td>[AKT3, B2M, CD3D, CD3G, HLA-A, MAPK14, NFKBIB,...</td>
      <td>15</td>
      <td>7.124730</td>
      <td>0.003280</td>
      <td>0.998825</td>
      <td>0.164768</td>
      <td>1.0</td>
      <td>0.300938</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>109</th>
      <td>hsa03009</td>
      <td>Ribosome biogenesis</td>
      <td>[DBR1, HSPA8, MDN1, NIP7, RBM19, REXO1, REXO4,...</td>
      <td>14</td>
      <td>6.734525</td>
      <td>0.004700</td>
      <td>0.998345</td>
      <td>0.188352</td>
      <td>1.0</td>
      <td>0.344978</td>
      <td>0.999925</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>339</th>
      <td>hsa05204</td>
      <td>Chemical carcinogenesis</td>
      <td>[]</td>
      <td>0</td>
      <td>1.549155</td>
      <td>1.000000</td>
      <td>0.202404</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>33</th>
      <td>hsa00410</td>
      <td>beta-Alanine metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.931795</td>
      <td>1.000000</td>
      <td>0.101269</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>31</th>
      <td>hsa00380</td>
      <td>Tryptophan metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.743145</td>
      <td>1.000000</td>
      <td>0.159219</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>349</th>
      <td>hsa05217</td>
      <td>Basal cell carcinoma</td>
      <td>[]</td>
      <td>0</td>
      <td>2.590350</td>
      <td>1.000000</td>
      <td>0.063135</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>0.999925</td>
    </tr>
    <tr>
      <th>69</th>
      <td>hsa00630</td>
      <td>Glyoxylate and dicarboxylate metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.237520</td>
      <td>1.000000</td>
      <td>0.264189</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>0.999925</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(c_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>153</th>
      <td>hsa04060</td>
      <td>Cytokine-cytokine receptor interaction</td>
      <td>[ACVR1, BMP6, BMP7, CCL24, CCR3, CCR9, CD4, CX...</td>
      <td>22</td>
      <td>6.795135</td>
      <td>0.000005</td>
      <td>1.000000</td>
      <td>0.000000</td>
      <td>1.0</td>
      <td>0.001835</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>369</th>
      <td>hsa05340</td>
      <td>Primary immunodeficiency</td>
      <td>[ADA, AICDA, BLNK, CD4, RFX5, TNFRSF13B]</td>
      <td>6</td>
      <td>0.824675</td>
      <td>0.000120</td>
      <td>0.999990</td>
      <td>0.010405</td>
      <td>1.0</td>
      <td>0.022020</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>317</th>
      <td>hsa05140</td>
      <td>Leishmaniasis</td>
      <td>[IFNGR2, IRAK4, ITGAM, MAPK12, MAPK13, NFKBIB,...</td>
      <td>8</td>
      <td>2.049320</td>
      <td>0.000455</td>
      <td>0.999925</td>
      <td>0.027655</td>
      <td>1.0</td>
      <td>0.055661</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>156</th>
      <td>hsa04064</td>
      <td>NF-kappa B signaling pathway</td>
      <td>[BLNK, EDARADD, ERC1, IL1R1, IRAK4, LYN, PLCG2...</td>
      <td>10</td>
      <td>4.053805</td>
      <td>0.004135</td>
      <td>0.998785</td>
      <td>0.196703</td>
      <td>1.0</td>
      <td>0.301142</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>150</th>
      <td>hsa04050</td>
      <td>Cytokine receptors</td>
      <td>[CCR3, CCR9, CXCR6, IFNGR2, IL1R1, IL20RA, IL3...</td>
      <td>10</td>
      <td>3.904215</td>
      <td>0.004380</td>
      <td>0.998690</td>
      <td>0.196703</td>
      <td>1.0</td>
      <td>0.301142</td>
      <td>1.0</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>110</th>
      <td>hsa03010</td>
      <td>Ribosome</td>
      <td>[]</td>
      <td>0</td>
      <td>1.253405</td>
      <td>1.000000</td>
      <td>0.280364</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>73</th>
      <td>hsa00730</td>
      <td>Thiamine metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.111140</td>
      <td>1.000000</td>
      <td>0.283889</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>128</th>
      <td>hsa03051</td>
      <td>Proteasome</td>
      <td>[]</td>
      <td>0</td>
      <td>1.280525</td>
      <td>1.000000</td>
      <td>0.267244</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>58</th>
      <td>hsa00563</td>
      <td>Glycosylphosphatidylinositol (GPI)-anchor bios...</td>
      <td>[]</td>
      <td>0</td>
      <td>0.842490</td>
      <td>1.000000</td>
      <td>0.412613</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>74</th>
      <td>hsa00740</td>
      <td>Riboflavin metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>0.161505</td>
      <td>1.000000</td>
      <td>0.848071</td>
      <td>1.000000</td>
      <td>1.0</td>
      <td>1.000000</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(i_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>124</th>
      <td>hsa03036</td>
      <td>Chromosome and associated proteins</td>
      <td>[AHDC1, AKAP9, ALDOC, ANAPC7, ANKRD17, ARID1A,...</td>
      <td>84</td>
      <td>57.924460</td>
      <td>0.000220</td>
      <td>0.999865</td>
      <td>0.040880</td>
      <td>1.000000</td>
      <td>0.080740</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>119</th>
      <td>hsa03021</td>
      <td>Transcription machinery</td>
      <td>[AFF1, ARID1A, ARID2, ATXN7, BRD4, CCNT1, CHD1...</td>
      <td>21</td>
      <td>10.080570</td>
      <td>0.000640</td>
      <td>0.999785</td>
      <td>0.058760</td>
      <td>1.000000</td>
      <td>0.117439</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>26</th>
      <td>hsa00310</td>
      <td>Lysine degradation</td>
      <td>[ALDH2, ALDH3A2, ASH1L, GCDH, HADHA, KMT2D, KM...</td>
      <td>9</td>
      <td>3.038860</td>
      <td>0.001385</td>
      <td>0.999820</td>
      <td>0.087762</td>
      <td>1.000000</td>
      <td>0.156433</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>113</th>
      <td>hsa03013</td>
      <td>RNA transport</td>
      <td>[AAAS, EIF2B1, EIF5B, NDC1, NUP155, NUP214, PY...</td>
      <td>11</td>
      <td>3.972560</td>
      <td>0.001705</td>
      <td>0.999535</td>
      <td>0.087762</td>
      <td>1.000000</td>
      <td>0.156433</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>374</th>
      <td>hsa05418</td>
      <td>Fluid shear stress and atherosclerosis</td>
      <td>[ACVR2A, ACVR2B, AKT2, CHUK, MAP2K5, MAPK7, NF...</td>
      <td>12</td>
      <td>4.984255</td>
      <td>0.002435</td>
      <td>0.999180</td>
      <td>0.092680</td>
      <td>1.000000</td>
      <td>0.178728</td>
      <td>0.999865</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>52</th>
      <td>hsa00534</td>
      <td>Glycosaminoglycan biosynthesis - heparan sulfa...</td>
      <td>[]</td>
      <td>0</td>
      <td>2.715120</td>
      <td>1.000000</td>
      <td>0.031885</td>
      <td>1.000000</td>
      <td>0.174461</td>
      <td>1.000000</td>
      <td>0.486884</td>
    </tr>
    <tr>
      <th>152</th>
      <td>hsa04054</td>
      <td>Pattern recognition receptors</td>
      <td>[]</td>
      <td>0</td>
      <td>1.622385</td>
      <td>1.000000</td>
      <td>0.190054</td>
      <td>1.000000</td>
      <td>0.518005</td>
      <td>1.000000</td>
      <td>0.996872</td>
    </tr>
    <tr>
      <th>51</th>
      <td>hsa00533</td>
      <td>Glycosaminoglycan biosynthesis - keratan sulfate</td>
      <td>[]</td>
      <td>0</td>
      <td>0.661145</td>
      <td>1.000000</td>
      <td>0.479118</td>
      <td>1.000000</td>
      <td>0.954243</td>
      <td>1.000000</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>73</th>
      <td>hsa00730</td>
      <td>Thiamine metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.075285</td>
      <td>1.000000</td>
      <td>0.291179</td>
      <td>1.000000</td>
      <td>0.712622</td>
      <td>1.000000</td>
      <td>0.999865</td>
    </tr>
    <tr>
      <th>55</th>
      <td>hsa00537</td>
      <td>Glycosylphosphatidylinositol (GPI)-anchored pr...</td>
      <td>[]</td>
      <td>0</td>
      <td>3.945295</td>
      <td>1.000000</td>
      <td>0.002980</td>
      <td>1.000000</td>
      <td>0.021366</td>
      <td>1.000000</td>
      <td>0.099423</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(ec_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>317</th>
      <td>hsa05140</td>
      <td>Leishmaniasis</td>
      <td>[IL4, MAPK14, NFKBIB, PRKCB, STAT1, TAB2, IFNG...</td>
      <td>14</td>
      <td>4.159845</td>
      <td>0.000005</td>
      <td>1.000000</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>226</th>
      <td>hsa04658</td>
      <td>Th1 and Th2 cell differentiation</td>
      <td>[CD3D, CD3G, IL12RB1, IL13, IL4, MAML3, MAPK14...</td>
      <td>20</td>
      <td>7.298305</td>
      <td>0.000010</td>
      <td>0.999995</td>
      <td>0.000470</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>153</th>
      <td>hsa04060</td>
      <td>Cytokine-cytokine receptor interaction</td>
      <td>[ACKR3, CCR9, IL12RB1, IL13, IL31, IL34, IL4, ...</td>
      <td>31</td>
      <td>13.442275</td>
      <td>0.000015</td>
      <td>0.999995</td>
      <td>0.000658</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>hsa00051</td>
      <td>Fructose and mannose metabolism</td>
      <td>[FUK, GMDS, HKDC1, MPI, PMM1, SORD, FUK, GMDS,...</td>
      <td>11</td>
      <td>2.958685</td>
      <td>0.000030</td>
      <td>1.000000</td>
      <td>0.001225</td>
      <td>1.000000</td>
      <td>0.002752</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>47</th>
      <td>hsa00520</td>
      <td>Amino sugar and nucleotide sugar metabolism</td>
      <td>[CYB5RL, FUK, GFPT2, GMDS, HKDC1, MPI, PMM1, F...</td>
      <td>12</td>
      <td>4.055615</td>
      <td>0.000170</td>
      <td>0.999980</td>
      <td>0.006955</td>
      <td>1.000000</td>
      <td>0.010748</td>
      <td>1.000000</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>81</th>
      <td>hsa00830</td>
      <td>Retinol metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>2.617205</td>
      <td>1.000000</td>
      <td>0.066540</td>
      <td>1.000000</td>
      <td>0.774800</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>82</th>
      <td>hsa00860</td>
      <td>Porphyrin and chlorophyll metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.553795</td>
      <td>1.000000</td>
      <td>0.200909</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>34</th>
      <td>hsa00430</td>
      <td>Taurine and hypotaurine metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>0.890655</td>
      <td>1.000000</td>
      <td>0.380428</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>349</th>
      <td>hsa05217</td>
      <td>Basal cell carcinoma</td>
      <td>[]</td>
      <td>0</td>
      <td>5.069310</td>
      <td>1.000000</td>
      <td>0.004680</td>
      <td>1.000000</td>
      <td>0.321138</td>
      <td>1.000000</td>
      <td>0.825542</td>
    </tr>
    <tr>
      <th>104</th>
      <td>hsa02042</td>
      <td>Bacterial toxins</td>
      <td>[]</td>
      <td>0</td>
      <td>0.084765</td>
      <td>1.000000</td>
      <td>0.921640</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(ei_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>226</th>
      <td>hsa04658</td>
      <td>Th1 and Th2 cell differentiation</td>
      <td>[CD3D, CD3G, IL12RB1, IL13, IL4, MAML3, MAPK14...</td>
      <td>19</td>
      <td>7.307075</td>
      <td>0.000050</td>
      <td>0.999975</td>
      <td>0.009475</td>
      <td>1.000000</td>
      <td>0.018350</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>26</th>
      <td>hsa00310</td>
      <td>Lysine degradation</td>
      <td>[ASH1L, EHMT1, KMT2A, KMT5A, SMYD1, SMYD3, WHS...</td>
      <td>16</td>
      <td>6.353000</td>
      <td>0.000235</td>
      <td>0.999945</td>
      <td>0.023565</td>
      <td>1.000000</td>
      <td>0.039146</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>334</th>
      <td>hsa05169</td>
      <td>Epstein-Barr virus infection</td>
      <td>[AKT3, B2M, CD3D, CD3G, HLA-A, MAPK14, NFKBIB,...</td>
      <td>27</td>
      <td>13.345100</td>
      <td>0.000320</td>
      <td>0.999900</td>
      <td>0.023565</td>
      <td>1.000000</td>
      <td>0.039146</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>124</th>
      <td>hsa03036</td>
      <td>Chromosome and associated proteins</td>
      <td>[AKAP9, ANAPC7, ANKS4B, ARHGEF10, ARID1A, ARMC...</td>
      <td>152</td>
      <td>118.538310</td>
      <td>0.000525</td>
      <td>0.999650</td>
      <td>0.027963</td>
      <td>1.000000</td>
      <td>0.048169</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>361</th>
      <td>hsa05235</td>
      <td>PD-L1 expression and PD-1 checkpoint pathway i...</td>
      <td>[AKT3, CD274, CD3D, CD3G, MAPK14, NFATC2, NFKB...</td>
      <td>19</td>
      <td>9.191205</td>
      <td>0.001015</td>
      <td>0.999660</td>
      <td>0.043419</td>
      <td>1.000000</td>
      <td>0.071259</td>
      <td>0.999975</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>104</th>
      <td>hsa02042</td>
      <td>Bacterial toxins</td>
      <td>[]</td>
      <td>0</td>
      <td>0.130180</td>
      <td>1.000000</td>
      <td>0.880446</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>81</th>
      <td>hsa00830</td>
      <td>Retinol metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>2.387655</td>
      <td>1.000000</td>
      <td>0.084105</td>
      <td>1.000000</td>
      <td>0.420211</td>
      <td>1.000000</td>
      <td>0.734914</td>
    </tr>
    <tr>
      <th>77</th>
      <td>hsa00770</td>
      <td>Pantothenate and CoA biosynthesis</td>
      <td>[]</td>
      <td>0</td>
      <td>3.101705</td>
      <td>1.000000</td>
      <td>0.021000</td>
      <td>1.000000</td>
      <td>0.245089</td>
      <td>1.000000</td>
      <td>0.428165</td>
    </tr>
    <tr>
      <th>250</th>
      <td>hsa04744</td>
      <td>Phototransduction</td>
      <td>[]</td>
      <td>0</td>
      <td>1.338225</td>
      <td>1.000000</td>
      <td>0.249809</td>
      <td>1.000000</td>
      <td>0.770701</td>
      <td>1.000000</td>
      <td>0.999975</td>
    </tr>
    <tr>
      <th>82</th>
      <td>hsa00860</td>
      <td>Porphyrin and chlorophyll metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>1.482415</td>
      <td>1.000000</td>
      <td>0.216759</td>
      <td>1.000000</td>
      <td>0.736863</td>
      <td>1.000000</td>
      <td>0.999975</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(ci_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>317</th>
      <td>hsa05140</td>
      <td>Leishmaniasis</td>
      <td>[IFNGR2, IRAK4, ITGAM, MAPK12, MAPK13, NFKBIB,...</td>
      <td>13</td>
      <td>4.095365</td>
      <td>0.000075</td>
      <td>0.999985</td>
      <td>0.014505</td>
      <td>1.000000</td>
      <td>0.014068</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>153</th>
      <td>hsa04060</td>
      <td>Cytokine-cytokine receptor interaction</td>
      <td>[ACVR1, BMP6, BMP7, CCL24, CCR3, CCR9, CD4, CX...</td>
      <td>29</td>
      <td>13.523935</td>
      <td>0.000110</td>
      <td>0.999945</td>
      <td>0.014505</td>
      <td>1.000000</td>
      <td>0.014068</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>369</th>
      <td>hsa05340</td>
      <td>Primary immunodeficiency</td>
      <td>[ADA, AICDA, BLNK, CD4, RFX5, TNFRSF13B, BLNK,...</td>
      <td>8</td>
      <td>1.576295</td>
      <td>0.000115</td>
      <td>0.999980</td>
      <td>0.014505</td>
      <td>1.000000</td>
      <td>0.014068</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>198</th>
      <td>hsa04350</td>
      <td>TGF-beta signaling pathway</td>
      <td>[ACVR1, BMP6, BMP7, FBN1, GREM1, INHBA, LTBP1,...</td>
      <td>21</td>
      <td>9.623055</td>
      <td>0.000200</td>
      <td>0.999940</td>
      <td>0.014505</td>
      <td>1.000000</td>
      <td>0.018350</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>226</th>
      <td>hsa04658</td>
      <td>Th1 and Th2 cell differentiation</td>
      <td>[CD4, DLL1, IFNGR2, IL13, IL4R, MAPK12, MAPK13...</td>
      <td>17</td>
      <td>7.208360</td>
      <td>0.000565</td>
      <td>0.999830</td>
      <td>0.024924</td>
      <td>1.000000</td>
      <td>0.041471</td>
      <td>0.999985</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>73</th>
      <td>hsa00730</td>
      <td>Thiamine metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>2.186425</td>
      <td>1.000000</td>
      <td>0.083050</td>
      <td>1.000000</td>
      <td>0.554494</td>
      <td>1.000000</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>51</th>
      <td>hsa00533</td>
      <td>Glycosaminoglycan biosynthesis - keratan sulfate</td>
      <td>[]</td>
      <td>0</td>
      <td>1.286265</td>
      <td>1.000000</td>
      <td>0.242634</td>
      <td>1.000000</td>
      <td>0.941195</td>
      <td>1.000000</td>
      <td>0.999985</td>
    </tr>
    <tr>
      <th>46</th>
      <td>hsa00515</td>
      <td>Mannose type O-glycan biosynthesis</td>
      <td>[]</td>
      <td>0</td>
      <td>2.905925</td>
      <td>1.000000</td>
      <td>0.029310</td>
      <td>1.000000</td>
      <td>0.357775</td>
      <td>1.000000</td>
      <td>0.705176</td>
    </tr>
    <tr>
      <th>81</th>
      <td>hsa00830</td>
      <td>Retinol metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>2.518320</td>
      <td>1.000000</td>
      <td>0.074920</td>
      <td>1.000000</td>
      <td>0.532925</td>
      <td>1.000000</td>
      <td>0.948121</td>
    </tr>
    <tr>
      <th>116</th>
      <td>hsa03018</td>
      <td>RNA degradation</td>
      <td>[]</td>
      <td>0</td>
      <td>4.608260</td>
      <td>1.000000</td>
      <td>0.008305</td>
      <td>1.000000</td>
      <td>0.132108</td>
      <td>1.000000</td>
      <td>0.319839</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python
display(eci_results)
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
      <th>id</th>
      <th>name</th>
      <th>candidate_features</th>
      <th>n_candidates_in_set</th>
      <th>mean_n_resample</th>
      <th>emp_p_e</th>
      <th>emp_p_d</th>
      <th>fdr_e</th>
      <th>fdr_d</th>
      <th>BH_fdr_e</th>
      <th>BH_fdr_d</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>226</th>
      <td>hsa04658</td>
      <td>Th1 and Th2 cell differentiation</td>
      <td>[CD3D, CD3G, IL12RB1, IL13, IL4, MAML3, MAPK14...</td>
      <td>28</td>
      <td>10.906870</td>
      <td>0.000010</td>
      <td>1.000000</td>
      <td>0.000545</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>317</th>
      <td>hsa05140</td>
      <td>Leishmaniasis</td>
      <td>[IL4, MAPK14, NFKBIB, PRKCB, STAT1, TAB2, IFNG...</td>
      <td>19</td>
      <td>6.205890</td>
      <td>0.000010</td>
      <td>1.000000</td>
      <td>0.000545</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>334</th>
      <td>hsa05169</td>
      <td>Epstein-Barr virus infection</td>
      <td>[AKT3, B2M, CD3D, CD3G, HLA-A, MAPK14, NFKBIB,...</td>
      <td>39</td>
      <td>19.527630</td>
      <td>0.000015</td>
      <td>1.000000</td>
      <td>0.000687</td>
      <td>1.000000</td>
      <td>0.001835</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>369</th>
      <td>hsa05340</td>
      <td>Primary immunodeficiency</td>
      <td>[CD3D, TNFRSF13B, ADA, AICDA, BLNK, CD4, RFX5,...</td>
      <td>10</td>
      <td>2.458895</td>
      <td>0.000130</td>
      <td>0.999970</td>
      <td>0.005809</td>
      <td>1.000000</td>
      <td>0.009542</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>153</th>
      <td>hsa04060</td>
      <td>Cytokine-cytokine receptor interaction</td>
      <td>[ACKR3, CCR9, IL12RB1, IL13, IL31, IL34, IL4, ...</td>
      <td>38</td>
      <td>20.171075</td>
      <td>0.000130</td>
      <td>0.999945</td>
      <td>0.005809</td>
      <td>1.000000</td>
      <td>0.009542</td>
      <td>1.000000</td>
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
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>75</th>
      <td>hsa00750</td>
      <td>Vitamin B6 metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>0.612195</td>
      <td>1.000000</td>
      <td>0.530487</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>104</th>
      <td>hsa02042</td>
      <td>Bacterial toxins</td>
      <td>[]</td>
      <td>0</td>
      <td>0.154295</td>
      <td>1.000000</td>
      <td>0.861841</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>118</th>
      <td>hsa03020</td>
      <td>RNA polymerase</td>
      <td>[]</td>
      <td>0</td>
      <td>2.140525</td>
      <td>1.000000</td>
      <td>0.109329</td>
      <td>1.000000</td>
      <td>0.623935</td>
      <td>1.000000</td>
      <td>0.866453</td>
    </tr>
    <tr>
      <th>74</th>
      <td>hsa00740</td>
      <td>Riboflavin metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>0.555795</td>
      <td>1.000000</td>
      <td>0.564242</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>82</th>
      <td>hsa00860</td>
      <td>Porphyrin and chlorophyll metabolism</td>
      <td>[]</td>
      <td>0</td>
      <td>2.239385</td>
      <td>1.000000</td>
      <td>0.098825</td>
      <td>1.000000</td>
      <td>0.609766</td>
      <td>1.000000</td>
      <td>0.866453</td>
    </tr>
  </tbody>
</table>
<p>367 rows × 11 columns</p>
</div>



```python

```
