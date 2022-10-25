These codelines automatically evaluate data related to Oxygen Evolution Reaction (OER) from biologic electrochemistry setup.

Packages that need to be installed beside typical ones are: eclabfiles and galvani. Those can be pip-installed.

To run properly you need to have PEIS measurements done in the beginning of the method and in every loop.

From PEIS measurement the program determines Re(Z) at absolute minimum of -Im(Z) and uses these values to correct determined potentials in following methods.

In the methodlist there can but mustn't be OCV, ZIR, LSV, CV, CP and Modular Potential present. ZIR is redundant with PEIS.  

The evaluation just needs two inputs: potential of the reference elektrode and area of the working electrode.

Have a look in the code to find out details.