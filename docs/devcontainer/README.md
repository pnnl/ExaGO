# Using ExaGO inside a VSCode Dev Container

## Prerequisites

1. install Docker Desktop and launch the app
2. install the "Remote Development" extension in VSCode
3. open your local clone of exago in VSCode
4. copy the CoinHSL 2019 tarball to the top level of your exago directory
5. if connected, disconnect from the PNNL VPN
6. launch the container build  
    * `cmd shift p` to open the command pallette in vscode
    * click `> Dev Container: rebuild and reopen container`
    * this will start building the container, taking about 40 minutes
    * click on the pop up with `(show log)` to view the progress
7. once the conatiner has build and launched, open `docs/devcontainer/tutorial.ipynb`
8. select the (newly built) existing jupyter kernel "ExaGO"
9. run all cells!
