# DMPlex dependencies for DMNetwork
This document lists the DMPlex dependencies for DMNetwork

DMNetwork holds a pointer to the DMPlex data structurre and uses it for mesh creation, querying sizes and topology information, partitioning and data (component) distribution. It uses the following functions from DMPlex

### Mesh creation

1. [DMPlexCreateFromCellList](https://www.mcs.anl.gov/petsc/petsc-3.6/docs/manualpages/DM/DMPlexCreateFromCellList.html)
2. [DMPlexCreateFromCellListParallel](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexCreateFromCellListParallel.html)

### Mesh sizes

1. [DMPlexGetChart](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetChart.html)
2. [DMPlexGetHeightStratum](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetHeightStratum.html)

### Partitioning and data distribution
1. [DMPlexGetPartitioner](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetPartitioner.html)
2. [DMPlexDistribute](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexDistribute.html)
3. [DMPlexDistributeData](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexDistributeData.html)

### Mesh topology query
1. [DMPlexGetSupportSize](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetSupportSize.html)
2. [DMPlexGetSupport](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetSupport.html)
3. [DMPlexGetCone](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMPLEX/DMPlexGetCone.html)
