
Dope_compare.py 

- Script to automate single, multi and loop model generation and DOPE comparison in MODELLER.
- Produces final plot to compare DOPE scores of different modelling methods: single template, multi-template, automatic loop refinement on multi-template model.
<br /> <br /> <br /> 

Requirements:
- Python modules: modeller, pylab, os, sys 
- Python 2.7



      usage:  template1 template2 chainid_template1 chainid_template2 query_ID

<br /> <br /> <br /> 
Notes:
- template1 should be the template with the highest sequence % id to query.
- The query_ID needs to be the same as the ID used to generate the .ali file. 
- chainid needs to be capitalised. 
- Depending on compute power n models could be changed - look for starting_model and ending_model. <br /> <br /> <br /> 



* to look out for; best DOPE score not necessarily best model - particularly for loop models, hint: view it in pymol.
