Structure-sum
========================================

This is a script to summarize the results from several structure runs. Works for outputs of Structure 2.2 using the admixture, no admixture, correlated and uncorrelated allele frequency models. It works for datasets containing recessive alleles (e.g. AFLP) and for datasets including information about populations of origin (but not for results from the USEPOPINFO option). It does not work for the linkage model.

# Original author

**Dorothee Ehrich**

* <https://en.uit.no/ansatte/person?p_document_id=41186>
* <https://www.researchgate.net/profile/Dorothee-Ehrich>

# Updated and maintained by

Vojtěch Zeisek, <https://trapa.cz/>.

# Homepage and reporting issues

<https://github.com/V-Z/structure-sum>

# License

GNU General Public License 3.0, see `LICENSE.md` and <https://www.gnu.org/licenses/gpl-3.0.html>.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# Note

Script and documentation is right now in state from 2011. It'll be updated to be fully compatible with R 4.3+ and more. You can also look at <https://github.com/MarekSlenker/structureSum>, which is another attempt to do similar job.

# Original documentation for Structure-sum

Version 2011

A series of R functions for summarizing the outputs of the program Structure ver. 2.3.3 (and some previous versions)

Dorothée Ehrich

Last revised 1.02.11

Structure-sum-2009 is a series of R functions which summarizes the output of the software Structure ver. 2.2 and 2.3 (<http://pritch.bsd.uchicago.edu/structure.html>, Pritchard et al. 2000), usually consisting of a series of runs.

If you did not use R before, you first need to install it. It can be downloaded from <http://www.r-project.org/> . R is a free software environment for statistical computing and graphics, which is very useful for many things. Numerous additional packages can be downloaded for different specialized analyses, notably also in population genetics. It is definitely worth discovering…

Second, download the file `Structure-sum-2009.R` from <http://www.nhm.uio.no/ncb> and copy it to any directory. Contrary to previous versions of the script, this version works for result files from the pc version as well as from the Bioportal of the University of Oslo.

The exact layout of the output files of Structure varies somewhat with the options chosen to run the program and the information contained in the input file. The Structure-sum2009.R script in its present version is compatible with outputs from files containing individual names as only information in addition to the data matrix, and files containing both individual and population names. Marker names, popflag and other additional information should not be used (might be added in a future version). The script works for runs using the admixture and non-admixture models, the RECESSIVE ALLELES option and correlated or uncorrelated allele frequencies, as well as the LOCPRIOR option, which is new in Structure 2.3 (Hubisz et al. 2009). For all other options, default values are expected. Structure-sum in its present version does not work for runs using the linkage model.

In order to use Structure-sum-2009.R, start R.

Choose "Source R code …" in the File Menu and open the file

Structure-sum-2009.R.

The next thing to do is to group all structure result files into one folder. It does not matter if they were produced by runs with different lengths, but they should be results for the same model (e.g. only no admixture model with uncorrelated allele frequencies), and of course from the same input file.

Then you have to make a list of these result files and the K they correspond to (e.g.in Excel or any text editor). The list should be saved as a text file into the same folder as the structure result files. It should be **tab separated**. When saving as text file from excel, chose tab separated. Sometimes it may be easier to rename the result files.

For example:

```
1	noAdunCor500000_run_1_f
1	noAdunCor500000_run_2_f
2	noAdunCor500000_run_3_f
2	noAdunCor300000_run_1_f
3	noAdunCor300000_run_2_f
```

etc…

Next, use "Change dir…" in the File Menu of R to choose the folder where your output files from Structure are. Then you are ready to use the functions.

The rapid way to calculate similarity coefficients in Structure-sum-2009.R is described in Ehrich D, Gaudeul M, Assefa A, Koch M, Mummenhoff K, Nemomissa S, Intrabiodiv Consortium, Brochmann C. 2007. Genetic consequences of Pleistocene range shifts: Contrast between the Arctic, the Alps and the East African mountains. Molecular Ecology, 16: 2542-2559.

This paper can be cited as a reference for the script.

Alternatively, it can be cited as a part of

Ehrich D. 2006. AFLPdat: a collection of R functions for convenient handling of AFLP data. Molecular Ecology Notes 6: 603-604.

Or cite it as Ehrich (2011) and refer to <http://tiny.cc/dorothee_ehrich>.

## Functions

### Structure.table

This function reads all result files after the list of files you made, and makes a table similar to the simulation summary in the front end version of Structure 2.3. The table contains also a column which records whether there were empty groups for this run. In addition, the function makes a plot of ln P(Data) in function of K. As input the function needs a list of result files and, if your input file contained a column with the population to which individuals belong, the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information. In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (If you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum).

```R
Structure.table ("list.txt", 0)
```

or

```R
Structure.table ("list.txt", pop = 5, locprior = T)
```

The output file is a table called lnPData.txt, which can be opened in excel, and a figure.

### Structure.simil

This function calculates a similarity coefficient among each pair of structure runs according to Nordborg et al. (2005). It makes a plot showing the average similarity coefficient for each K with standard deviations. In the files simil_coefficients.txt the average coefficient of similarity with standard deviation is listed for each K. In addition a file called similK (where K has the different values of K used) is produced for each K. These files contain the individual coefficients of similarity for the comparison between each pair of runs. Structure-sum-2009.R uses a slightly different approach to align the matrices than the approach described in Norborg et al. 2005. The matrices are aligned column by column, finding for each column of matrix 1 the column of matrix 2 which is most similar. This approach is much faster than taking the minimum over all possible permutations, especially for runs with a large number of groups. The approach is described in the appendix of Ehrich et al. (2007).

As input the function needs a list of result files and, if your input file contained a column with the population to which individuals belong, the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information.

```R
Structure.simil ("list.txt", 0)
```

In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (if you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum). If the number of localities used in the analysis was not the same as the number of populations used in the input file, this needs also to be entered into the function with the parameter `nbloc`.

```R
Structure.simil ("list.txt", pop = 16, locprior = T, nbloc = 5)
```

For this example a file was used where 16 populations were collected in 5 different regions, and the 5 different regions were used as locality prior in the Structure analysis.

The output files are tables called simil_coefficients.txt and similK.txt, which can be opened in Excel, and a figure.

### Structure.deltaK

This function makes four plots which can be used to determine the optimal K after the method of Evano et al. (2005). Read the paper to see what he does and why ☺

As input the function needs a list of result files and, if your input file contained a column with the population to which individuals belong, the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information. . In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (If you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum).

```R
Structure.deltaK ("list.txt", 0)
```

The output is a figure and a file called DeltaK.txt (can be opened in Excel).

### Structure.order

This function orders the output of several runs of Structure obtained for the same value of K in order to place the same cluster in the same column for all the runs. The ordering is done according to the simplified approach described for the function Structure.simil, and not according to the permutation approach of Rosenberg et al. (2002) or Nordborg et al. (2005). The input file for the function is a list of runs similar to that used by the previous functions, except that all runs should be carried out for the same value of K. In addition, as for the previous functions, if your input file contained a column with the population to which individuals belong, it needs the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information.

```R
Structure.order ("list.txt", 0)
```

In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (if you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum). If the number of localities used in the analysis was not the same as the number of populations used in the input file, this needs also to be entered into the function with the parameter `nbloc`.

```R
Structure.simil ("list.txt", pop = 16, locprior = T, nbloc = 5)
```

As output, the function writes a text file with the matrix part of each result file. The files have the same name as the result files, but –ord.txt is added at the end of the file name. In addition the function produces a text file with a matrix of the average q values for each individual in each cluster, resulting from averaging over the compared runs (as used for example by Nordborg et al. 2005). This function is thus similar to the software CLUMPP (Jakobsson and Rosenberg 2007).

### Structure.cluster

This function calculates a coefficient of clusteredness for Structure results after Rosenberg et al. (2005). This coefficient measures the extent to which individuals were estimated to belong to a single cluster rather than to a combination of clusters.

As input the function needs a list of result files and, if your input file contained a column with the population to which individuals belong, the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information.

```R
Structure.cluster ("list.txt", 0)
```

In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (if you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum). If the number of localities used in the analysis was not the same as the number of populations used in the input file, this needs also to be entered into the function with the parameter `nbloc`.

```R
Structure.simil ("list.txt", pop = 16, locprior = T, nbloc = 5)
```

The output of this function is a text file containing the clusteredness value for each run, named clustered_coefficients.txt.

Structure.groups

This function adds a column to outputs of Structure which contains a number for the group each individual was assigned to. The output contains only the matrix without individual or population names and the additional column. This is convenient for example to plot the results in a GIS program.

As input the function needs the name of the file to which you want to add a column with cluster numbers, and, if your input file contained a column with the population to which individuals belong, the number of populations (this is not K!). The default for the number of populations is 0, corresponding to a file without population information.

```R
Structure.groups ("noAdunCor500000_run_1_f", 0)
```

In addition, if you used the LOCPRIOR option, you need to indicate this by typing `locprior = T` in the parenthesis specifying the arguments to the function (if you did not use LOCPRIOR, it is not necessary to add anything, as the "old" model is the default in Structure-sum). If the number of localities used in the analysis was not the same as the number of populations used in the input file, this needs also to be entered into the function with the parameter `nbloc`.

Ex: Structure.simil ("noAdunCor500000_run_1_f", pop = 16, locprior = T, nbloc = 5)

The output is a table containing the part of the Structure output with the inferred ancestry of each individual and the new column with a different number for each group (can be opened in Excel).

## References

* Evano G, Regnaut S and Goudet J (2005) Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Mol. Ecol. 14: 2611-2620.
* Jakobsson M and Rosenberg NA (2007) CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics 23:1801-1806.
* Hubisz KJ, Falush D, Stephens M and Pritchard JK (2009) Inferring weak population structure with the assistance of sample group information. Molecular Ecology Ressources 9: 1322-1332.
* Nordborg M, Hu TT, Ishino Y, Jhaveri J, Toomajian C, et al. (2005) The pattern of polymorphism in Arabidopsis thaliana. PLoS Biol 3(7): e196.
* Pritchard JK, Stephens M and Donnelly PJ (2000) Inference of population structure using multilocus genotype data. Genetics 155: 945-959.
* Rosenberg NA, Mahajan S, Ramachandran S, Zhao C, Pritchard JK, et al. (2005) Clines, clusters, and the effect of study design on the inference of human population structure. PLoS Genet 1(6): e70.
* Rosenberg NA, Pritchard JK, Weber JL, Cann HM, Kidd KK, Zhivotovsky LA and Feldman MW (2002) The genetic structure of human populations. Science, 298: 23812385.

Good luck!

Dorothee Ehrich

<dorothee.ehrich@ib.uit.no>

(feel free to complain, make suggestions and ask questions)

