Basic usage
------------------------------
From the command line:
	python3.2 run.py
will execute the program and use a basic CLI interface to configure the problem

Options:
	python3.2 run.py (y/n) [(mode) (mode opts)]
y/n: selects whether to auto execute all modules
mode:	0 = Guderley
	1 = LILAC
	2 = HYADES
mode opts: see below

To run a HYADES implosion:
	python3.2 run.py y 2 {fname}
{fname} -> name of HYADES netCDF file

Old functionality not currently supported
------------------------------
To run a Guderley implosion:
	python3.2 run.py y 0 {R} {tc} {xi} {D2} {3He} {ei}
replace:
{R} -> capsule radius (um)
{tc} -> shock collapse time (ns)
{xi} -> shock strength (um/ns^a)
{D2} -> atm of D2
{3He} -> atm of 3He
{ei} -> Rygg-style ei coupling.

To run a LILAC implosion:
	python3.2 run.py y 1 {fname}
{fname} -> name of LILAC file