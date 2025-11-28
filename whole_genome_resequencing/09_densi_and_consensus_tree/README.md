# Consensus tree and visualisation

The previous step produced a large collection of phylogenetic trees, each representing the relationships among individuals for a different genomic region. Two approaches can provide an overall view of the genetic relationships among these individuals: (i) constructing a consensus tree that summarizes information across all individual trees, and (ii) visualizing the entire set of trees using a DensiTree graph.

### 1. Consensus tree
We can generate a consensus species tree that integrates information from all individual gene trees using **ASTRAL**. Using the `project.trees` file containing all ML trees, the ASTRAL species tree can be estimated with the following command:

```bash
#!bin/bash

cd ~/project

# Astral consensus tree
java -jar  astral.5.6.3.jar -i ./trees/project.trees -o./astral/project.astral.contree´
```

This consensus tree can now be visualized with your program of choice, such as [FigTree](https://tree.bio.ed.ac.uk/software/figtree/).

Note that the branch lengths in the ASTRAL output tree are not standard evolutionary distances like in a typical phylogenetic tree inferred from sequence data, but are instead measured in coalescence units and reflect gene tree concordance. An informative measure to infer the concordance at which individuals are within the same clade is the calculation of gene concordance factors (gCFs), which indicate the percentage of gene trees supporting each branch. Using the project.trees file and the ASTRAL species tree, gCFs can be calculated in IQ-TREE using the following script:

```bash
#!bin/bash

cd ~/project

# Define in and output trees
TREE_ASTRAL="~/project/astral/project.astral.contree"
TREES_WINDOWS="~/project/trees/project.trees"
OUT="~/project/astral"

# Calculate gCF
iqtree2 -t "$TREE_ASTRAL" --gcf "$TREES_WINDOWS" -pre "$OUT"
```



