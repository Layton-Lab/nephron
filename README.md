# About
This is latest version of sex-specific mathematical models for epithelial transport along the nephron implemented in Python 3. This includes the normal pregnancy study. Related research papers are listed below. Please cite appropriately. To speed up the computation time, parallel computation of different types of nephron is implemented for the multiple nephron model. 

# Instructions
To run the parallel simulation code use command: **python3 parallel_simulate.py --sex [option] --species [option] --type [option] --diabetes [option] --inhibition [option] --pregnant [option]**

The options here are:

sex: **Male, Female** (required);

species: **human, rat** (required);

type: **superficial, multiple** (required);

diabetes: **Severe, Moderate, Non** (optional, default: Non);

pregnant: **mid, late** (optional, default: non, only for female rat);

inhibition: **ACE, SGLT2, NHE3-50, NHE3-80, NKCC2-70, NKCC2-100, NCC-70, NCC-100, ENaC-70, ENaC-100, SNB-70, SNB-100** (optional, default: None).

unx: **N, Y** (optional, default: N)

Notes:
* Human only have ACE and SGLT2 inhibition cases. The others are for rats.
* pregnancy: only has been characterized for normal pregnant rat superficial nephron at this time (i.e., not done for humans and for diabetes, also multiple nephron)

### Understanding output

All the output files' names are in following structure: 'sex_species_segment_concentration/flow_of_solute_in_compartment.txt'. 

Here is an example: female_rat_ccd_con_of_Cl_in_Bath.txt. It contains interstitial concentration of Chloride along cortical collecting duct in female rat.

Another example: male_hum_pt_flow_of_Na_in_Lumen.txt. It contains luminal flow of Sodium along proximal convolute tubule in male human.

These results are scaled per nephron.

The unit of concentration from outputs is **mmol/L (mM)**.

The unit of volume is **nl/min**.

The unit of flow is **pmol/min**.

**/plot/** contains scripts for plotting output

## Research Papers
Please cite appropriate paper(s) when using this model.
Published papers related to/using model:

* **superficial nephron (sex-specific):** [2019 Hu et al. "Functional implications of the sex differences in transporter abundance along the rat nephron: modeling and analysis"](https://journals.physiology.org/doi/full/10.1152/ajprenal.00352.2019)
* **multiple nephron (sex-specific):** [2020 Hu et al. "Sex differences in solute transport along the nephrons: effects of Na+ transport inhibition"](https://journals.physiology.org/doi/abs/10.1152/ajprenal.00240.2020?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org)
* **human multiple nephron (male only):** [2019 Layton and Layton "A computational model of epithelial solute and water transport along a human nephron"](https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1006108)
* **human multiple nephron (sex-specific):** [2021 Hu et al. "Sex differences in solute and water handling in the human kidney: Modeling and functional implications"](https://www.sciencedirect.com/science/article/pii/S2589004221006350)
* **diabetic human (sex-specific):** [2021 Hu et al. "A Computational Model of Kidney Function in a Patient with Diabetes"](https://www.mdpi.com/1422-0067/22/11/5819)
* **pregnant rat:** (under review)

## Previous versions
Previous versions of this model are on [this](https://github.com/uwrhu) page.
