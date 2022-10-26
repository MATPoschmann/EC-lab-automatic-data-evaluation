These codelines automatically evaluate data related to Oxygen Evolution Reaction (OER) from biologic electrochemistry setup.

Packages that need to be installed beside typical ones are: eclabfiles and galvani. Those can be pip-installed.
To run properly you need to have PEIS measurements done in the beginning of the method and in every loop.

It opens every *.mpr-file in its own folder and does automated data evaluation as follows:
From PEIS measurement the program determines Re(Z) at absolute minimum of -Im(Z) and uses these values to correct determined potentials in following methods.
It IR-corrects Ewe data by Re(Z) * I and the potential of the reference electrode if entered.
It creates Tafel-plots and determines Tafelslope from Modular Potentio measurements as lowest fittable slope to 4 points within the measurement range.
It generates a lot of graphics to investigate the data to compare the electrode performance in every loop.
It saves the data of CP, CV and MP in CSV format and arrangement which is easily processable in plotting programs like Qti-Plot, Origin ...etc.

In the methodlist there can but mustn't be OCV, LSV, CV, CP and Modular Potentio (MP) present.  

The evaluation just needs two inputs: potential of the reference elektrode and area of the working electrode.

Have a look in the code to find out details.
