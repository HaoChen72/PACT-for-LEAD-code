# **Feature extraction code for PACT**



##### **Description**

###### This code accompanies manuscript entitled "Functional Photoacoustic Foot Imaging Enables Ischemia Detection and Lesion Localization in Lower Extremity Arterial Disease." These codes are used to automatically extract the feature values. All extracted feature values for all participants can be found in the manuscript associated data (which will be available after peer review). This repo serves as a platform for running 2D and 3D feature extraction using an example participant's data -"**exampleData.zip**". "**Main\_calculate\_2D.m**" code is used to extract vessel density, number of vessels, vessel total area and unit area, vessel width, and relative blood volume inferred from PA intensity. "**Main\_calculate\_3D.m**" code is used to calculate number of branches, number of branch points, distance metrics, sum of angle means, max vessel lengths, and average branch vessel lengths. Details see manuscript. We did not provide analysis for the oxygen saturation code as that is associated with an ongoing manuscript where the used spectrum unmixing method is central to.



##### **System requirement**

###### All software and data were tested on MATLAB 2025b and MATLAB 2023b. All additional toolboxes were installed. No other non-standard hardware is required.



##### **Installation guide**

###### Download all the files from Github to one folder. No specific installation is required.



##### **Demo and instructions for use**

###### To run the data, download the "exampleData.zip", unzip the file, and copy the path to "root\_path" variable for both Main\_calculate\_2D.m and Main\_calculate\_3D.m. Make sure the path is until the /exampleData/ folder path level, and only until /exampleData/ folder path. Run Main\_calculate\_2D.m will generate an .xlsx file named as vessel\_analysis\_2D.xlsx, and run Main\_calcualte\_3D.m will generate an .xlsx file named as vessel\_analysis\_3D.xlsx for the given participant's data. The calculation time should not be over 2 minutes for a modern computer.



##### **Data description**

###### In the shared example data, each foot contains two .mat files, i.e., resize\_divide.mat and top\_view\_seg.mat.



###### resize\_divide.mat:

###### The dataset contains three variables for 2D analysis, resized\_figs consists of a 640×1536×30 double-precision array containing the original 2D cross sectional PA images, resized\_labels is a 640×1536×30 logical array containing the binarized labels after smart skin removal, resized\_masks is a 640×1536×30 logical array containing the deep learning model extracted foot vasculature. In all three variables, the first two dimensions correspond to the image size, while the last dimension represents the number of 2D cross-sectional PA images.



###### top\_view\_seg.mat:

###### The dataset contains a 901×1600 double-precision array representing the top-view maximum intensity projection image of the foot for 3D analysis.

