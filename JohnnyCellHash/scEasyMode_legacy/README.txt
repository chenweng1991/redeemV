this package implements:

- multiseq correction in python using a within-barcode zscore correction
- plotting for stacked barplots in your dataset
- mouse cell filtering/separation from mixed dataset
- scanpy wrapper that simplifies the workflow 

to use this:

- clone the package into your directory by running:
	git clone https://github.com/johnnyUCSF/scEasyMode
- then load modules:
	from scEasyMode import mousefilter
	from scEasyMode import clusterplot
	from scEasyMode import pymulti
	from scEasyMode import sceasy
- then call functions from modules, example:
	pymulti.pymulti(R1,R2,bcsmulti,bcs10x,split=True)

** note that parts of this is taken from various parts of the internet
