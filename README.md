# Annotation of V1, V2, and V3 in the Human Connectome Project

This repository is a tool for annotating the iso-eccentricity and visual area boundary contours in V1, V2, and V3 that were initially annotated by [Benson **et al**. (2021)](https://doi.org/10.1101/2020.12.30.424856).

## How to use this repository

### Before you start

Using this repository requires signing up for GitHub and downloading/running some software, but all of these are free.

1. You must have a GitHub account ([join here](https://github.com/join)) and you must have `git` installed ([instructions here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)).
2. You must have Docker installed ([instructions here](https://docs.docker.com/get-docker/)). Docker is a virtual-machine tool that allows you to execute code packaged up in docker images by other authors without having to share an OS or software dependencies.
3. These instructions are written for use with a BASH shell (which is usually the default for Mac OS and many flavors of Linux). You may have to change them slightly if you are using a different shell.

### Instructions

1. Make a fork of this repository into your own GitHub profile. This is done using the "Fork" button in the upper-right corner of this repository's GitHub page.
2. In a BASH shell, clone your own version of this repository.

   ```bash
   $ git clone https://github.com/<username>/hcp-annot-v123
   ```
   
   Then `cd` into the new directory.
   
   ```bash
   $ cd hcp-annot-v123
   ```
3. Run `docker-compose` to start the virtual-machine.

   ```bash
   $ docker compose up
   ```
   
   This command will likely take a long time to run the first time you run it, and it will produce a lot of output. Near the end of the output will be something that looks like this:
   
   ```
   neuropythy_1  |     To access the notebook, open this file in a browser:
   neuropythy_1  |         file:///home/jovyan/.local/share/jupyter/runtime/nbserver-9-open.html
   neuropythy_1  |     Or copy and paste one of these URLs:
   neuropythy_1  |         http://ba0d0b5346d7:8888/?token=3e27b357b6a63752d130ae22ac0c160d8e26ec252c4d3c66
   neuropythy_1  |      or http://127.0.0.1:8888/?token=3e27b357b6a63752d130ae22ac0c160d8e26ec252c4d3c66
   ```
   
   This indicates that the virtual machine is running.
4. Copy-and-paste the web address starting with `http://127.0.0.1:8888/` into your web-browser. This should bring up a [Jupyter](https://jupyter.org/) interface.
5. Click to open the file `annotate.ipynb`. This should open a new window containing a Jupyter notebook.
6. Select the first cell of the notebook and press Control+Enter. This should create an output cell below the original cell that contains the annotation tool. The use of this tool is described elsewhere. The tool automatically saves your progress as you go.
7. Once you have finished annotating (you do not need to be completely finished---just finished for the time being), you can close the Jupyter browser tabs and press Control+C in the terminal window that is running the virtual machine to shut it down.
8. To commit your work to GitHub you can use the following command (from within the root of the `hcp-annot-v123` repository:

   ```bash
   $ git commit -a -m '<some message>'
   ```
   
   The `<some message>` can be any message you wish to attach to the batch of contours worked on since your last commit. It's fine for this project to commit an empty message (`''`)---comments on individual contours should go in the Notes section of the contour editor.
   
   **Important**: You should commit your work often, ideally after every time you finish working on a set of contours. Committing the work mostly just prevents work from being lost.
