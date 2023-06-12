# GreenEx_Py installation manual for Windows

Currently (update 10-06-23), the python module can only be locally installed by cloning the GitHub repository and installing the corresponding python environment. Please follow the steps below;

**Cloning GitHub repository using Git Bash** 
<br>For additional information, please check this [link](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository?tool=webui)

Prerequisites:
- [Git Bash](https://git-scm.com/) installed

Steps:
1.	On GitHub.com, navigate to the main page of the repository.
<br>```https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/GreenEx_Py```
2.	Above the list of files, click Code.
3.	Copy the URL for the repository.
4.	Open Git Bash.
5.	Change the current working directory to the location where you want the cloned directory.
6.	Type git clone, and then paste the URL you copied earlier.
<br>```$ git clone https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/GreenEx_Py.git```
7.	Press Enter to create your local clone.

<br><br>**Cloning GitHub repository using GitHub Desktop app** 
<br>For additional information, please check this [link](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository?tool=desktop)

Prerequisites:
 - [GitHub desktop](https://desktop.github.com/) installed

Steps:
1.	On GitHub.com, navigate to the main page of the repository.
<br>```https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/GreenEx_Py```
2.	Above the list of files, click Code.
3.  To clone and open the repository with GitHub Desktop, click  Open with GitHub Desktop.
4.  Follow the prompts in GitHub Desktop to complete the clone.

<br><br>**Installing and using conda environment** 
<br>For additional information, please check this [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

Prerequisites:
- Either [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed
- GitHub repository is already cloned

Steps:
1.	Open anaconda or miniconda prompt
2.	Navigate to the cloned GitHub repository directory using ‘cd’ command, example;
<br>```cd Documents/GitHub/GreenEx_Py```
3.	Create the environment from the GreenExpPy.yml file:
<br>```conda env create -f GreenExpPy.yml```
4.	Activate the new environment: 
<br>```conda activate Greenex_py```
5.	Open jupyter notebook by typing:
<br>```jupyter notebook```
6.	When jupyter notebook is opened, navigate to cloned GitHub repository 
7.	Use example.ipynb as an example to run the code
