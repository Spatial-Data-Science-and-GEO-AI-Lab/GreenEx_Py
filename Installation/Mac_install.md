# GreenEx_Py installation manual for Windows

Currently (update 10-06-23), the python module can only be locally installed by cloning the GitHub repository and installing the corresponding python environment. Please follow the steps below;

**Cloning GitHub repository**
<br>For additional information, please check this [link](https://blogs.sap.com/2019/07/12/how-to-clone-a-github-repository-to-local-mac-computer/)

Steps:
1. Open the main page of the repository in browser. click Clone or download.
<br>```https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/GreenEx_Py```
2. Above the list of files, click Code. 
3. Copy the URL for the repository. 
4. Open Terminal on your mac. 
5. Type “cd” and the directory where you want the cloned directory to be made. You can right-click the folder in Finder and choose “Copy 'the folder name'” to copy the path into clipboard. Then by pressing “Command” and “v” on your keyboard to paste the path into terminal.
6. Type “git clone”, and then paste the URL you copied in step 2. Press Enter. The local clone will be created.

<br><br>**Installing and using conda environment** 
<br>For additional information, please check this [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

Prerequisites:
- Either [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed
- GitHub repository is already cloned

Steps: 
1. Open Terminal 
2. Navigate to the cloned GitHub repository directory using ‘cd’ command, example; 
<br>```cd Documents/GitHub/GreenEx_Py```
3. Create the environment from the GreenExpPy_Mac.yml file: 
<br>```conda env create -f GreenExpPy_Mac.yml```
4. Activate the new environment: 
<br>```conda activate Greenex_py```
5. Open jupyter notebook by typing: 
<br>```jupyter notebook```
6. When jupyter notebook is opened, navigate to cloned GitHub repository 
7. Use example.ipynb as an example to run the code


**NOTE: Local issues may be encountered when running the module locally on a Mac device, as was the case during testing. The module is developed on a Windows device and a separate yml file, containing the python environment, had to be created for Mac.** 
