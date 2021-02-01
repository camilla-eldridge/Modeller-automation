# Automating protein homology modelling with Modeller 

Scripts to automate protein homology modelling with `Modeller`.             
<br /> <br /> <br /> 


`fasta2ali.py`

- Converts fasta file to ali format for `Modeller`**

            usage: fasta2ali.py ID fasta_file
                     
            
            
            
<br /> <br /> <br /> 

`Dope_compare.py` 

- Script to automate single, multi and loop model generation and DOPE comparison in `Modeller`.
- Produces final plot to compare DOPE scores of different modelling methods: single template, multi-template, automatic loop refinement on multi-template model.
<br /> <br /> <br /> 



![alt text](https://github.com/camilla-eldridge/Modeller-scripts/blob/master/DOPE_compare/dope_example.png)

<br /> <br /> <br /> 

Requirements:
- Python modules: `modeller`, `pylab`, `os`, `sys` 
- `Python 2.7`

<br /> <br /> <br /> 



      usage:  template1  template2   chainid_template1  chainid_template2  query_ID

<br /> <br /> <br /> 
Notes:
- template1 should be the template with the highest sequence % id to query.
- The query_ID needs to be the same as the ID used to generate the .ali file**
- chainid needs to be capitalised. 
- Depending on compute power n models could be changed - look for starting_model and ending_model. 
- To look out for; best DOPE score not necessarily best model - particularly for loop models, hint: **view it in pymol.**
