from setuptools import setup

setup(
   name='chembl_grounder',
   version='0.6',
   description='Grounding of a DailyMed ingredient (with a UNII and a name) to ChEMBL using G-SRS to fetch names and structures for a UNII. These are matched to ChEMBL using query indexes.',
   author='Tim Rozday',
   author_email='timrozday@ebi.ac.uk',
   packages=['chembl_grounder'],  #same as name
)
