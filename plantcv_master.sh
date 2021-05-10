###
# @author: tvanbutselaar
# Current version: 20201217
# Changes: HPC compatibility

# Working script for PlantCV image analysis tool to measure rosette area in multiplant pictures. Information on https://plantcv.readthedocs.io/. This script encompasses 5 EOF files that form the main functions, and a final for loop that calls the first EOF bash script for all pictures in input folder.

# It is important to consider the limitations of this pipeline. The pipeline will have difficulties in the following cases:
	# When separate plants in the picture are touching each other or overlapping each other
	# When there is too much moss on the pots
	# When there are still plant labels present (especially those that light up in the green-magenta channel)
	# When the plants are not distributed in a regular grid and this grid is not similar across pictures (USE A TRAY FOR THIS)
	# When there are plants of multiple categorical variables (genotype/treatment) in one picture

# The script will analyze all the files that are in an input folder within the current working directory
# Execute script with cwd as directory where your input folder with images is. Execute the script as follows: bash plantcvbash.sh [name of input folder in pwd] [name of output folder in pwd]

# Dependencies of this script: python3 (tested and used with Python 3.6.7) with following modules installed. All can be installed through pip3 install --user [module]. pip3 should also install necessary dependencies of PIL and plantcv
	# PIL
	# plantcv
	# pandas
	# Seaborn

# Subsequent analysis will prep a figure with bars for every category in the experiment. This part is dependent on a notes.txt file, which should contain: 
	# picture names and their corresponding genotypes (also include a stylized genotype if you want), this all tab-separated. The stylized genotype should be the last part of this line. Each picture on new line. 
	# a line describing the squared pixel to squared mm ratio. the line must contain the word 'Pixel ratio' (Case-sensitive) and must be ending with a colon followed by the squared pixel to squared mm ratio.
		# Protips: 
			# flank your text with $ if you want the text to be printed in italics in the plot
			# use mathregular-based expressions if you want to use subscript, superscript, or special characters in the plot: https://matplotlib.org/3.3.3/tutorials/text/mathtext.html
		
# Use a representative picture and analyze in Photoshop to estimate the parameters in the variable list below. These parameters will be called on in the image analysis pipeline. These are the primary parameters that need to be adjusted for your images. The main variables to be altered are those for ROI, Cluster Nrow/Ncol, Rotation, and Shift numbers. Run the program locally and use output from the first image to further estimate other parameters

###

mpwd=$(pwd)

#Use a representative picture and analyze in Photoshop to estimate the parameters in the variable list below. These parameters will be called on in the image analysis pipeline. These are the primary parameters that need to be adjusted for your images. The main variables to be altered are those for ROI, Cluster Nrow/Ncol, Rotation, and Shift numbers. Run the program locally and use output from the first image to further estimate other parameters

cat <<\EOF > "$mpwd"/$2_variables.txt
0 ##			Script 1+2 ROI x_adj			the top-left pixel defining the rectangle of ROI
0 ##			Script 1+2	ROI y_adj			the top-left pixel defining the rectangle of ROI
1700 ##	Script 1+2	ROI h_adj			value for how far along y-axis from top-left pixel the ROI stretches
3400 ##	Script 1+2	ROI w_adj			value for how far along x-axis from top-left pixel the ROI stretches
4 ##		Script 1		cluster nrow		value for number of rows of cluster. if plants parts are identified okay, but are not declustered correctly, this is the main parameter to tweak
6 ##		Script 1		cluster ncol		value for number of columns to cluster. if plants parts are identified okay, but are not declustered correctly, this is the main parameter to tweak
-0.3 ##		Script 1		rotation_deg		rotation angle in degrees, can be a negative number, positive values move counter clockwise
100 ##		Script 1		shift number		value for shifting image up x pixels, must be greater than 1
150 ##		Script 1		shift number		value for shifting image left x pixels, must be greater than 1
120 ##		Script 1		bin treshold		threshold value (0-255) for binary tresholding, higher values will generate more background. if missing a lot of plant parts, or you have too much background, this is the main parameter to tweak
200 ##		Script 1		fill size				minimum pixel size for objects, those smaller will be filled
1 ##			Script 1		dilation kernel 	an odd integer that is used to build a ksize x ksize matrix using np.ones. Must be greater than 1 to have an effect. Greater values will ensure all plant parts are included, but also will overestimate size of plant
1 ##			Script 1		dilation iter		number of consecutive filtering passes for dilation. Greater values will ensure all plant parts are included, but also will overestimate size of plant
EOF

#Generate main python script for declustering pictures
cat <<\EOF > "$mpwd"/$2_python1.py
#!/usr/bin/python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys, traceback
import os
import re
import numpy as np
import argparse
import string
import glob
import shutil
from PIL import Image
from plantcv import plantcv as pcv

### Parse command-line arguments
def options():
	parser = argparse.ArgumentParser(description="Imaging processing with opencv")
	parser.add_argument("-i", "--input", help="Input image directory.", required=True)
	parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=True)
	parser.add_argument("-n", "--names", help="path to txt file with names of genotypes to split images into", required =False)
	parser.add_argument("-r","--result", help="result file.", required= False )
	parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action=None)
	parser.add_argument("-v", "--variables",help="List of variables to be passed to pipeline", required=True)
	args = parser.parse_args()
	return args

##Main declustering function
def main(img8,tray, args, pwd):
	# Get options
	args = options()
	if args.names!= None:
		args.names="../../names/%s.txt" %(str(tray))
		print("Names file: %s" %(args.names))
	else:
		print("No names file submitted for %s" %(str(tray)))
		
	#Get variables from variable arguments
	with open(args.variables) as f:
		roi1=[line.split(' ')[0] for line in f]

	# Read image
	img, path, filename = pcv.readimage(img8)
	pcv.params.debug=args.debug #set debug mode

	## STEP 2: Normalize the white color so you can later compare color between images. Inputs
	img1 = pcv.white_balance(img,roi=None)
	
	## STEP 3: Rotate the image
	rotate_img = pcv.rotate(img1, float(roi1[6]), True)
	
	## STEP 4: Shift image. This step is important for clustering later on. The resulting image is the same size as the original.
	shift1 = pcv.shift_img(rotate_img, int(roi1[7]), "bottom")
	shift2 = pcv.shift_img(shift1, int(roi1[8]), "right")
	img1 = shift2

	## STEP 5: Convert image from RGB colorspace to LAB colorspace. Keep only the green-magenta channel (grayscale).
	a = pcv.rgb2gray_lab(img1, "a")

	## STEP 6: Set a binary threshold on the saturation channel image.
	img_binary = pcv.threshold.binary(a, int(roi1[9]), 255, "dark")

	## STEP 7: Fill in small objects (speckles)
	fill_image = pcv.fill(img_binary, int(roi1[10]))

	## STEP 8: Dilate so that you don't lose leaves (just in case). Only do this if dilation parameters allow this.
	if int(roi1[11]) > 1 and int(roi1[12]) > 1:
		dilated = pcv.dilate(fill_image, int(roi1[11]), int(roi1[12]))
	else:
		print("Dilation kernal and iter <= 1, skipping dilation...")
		dilated=fill_image

	## STEP 9: Find objects (contours: black-white boundaries)
	id_objects, obj_hierarchy = pcv.find_objects(img1, dilated)

	## STEP 10: Define region of interest (ROI)
	roi_contour, roi_hierarchy = pcv.roi.rectangle(img1, int(roi1[0]), int(roi1[1]), int(roi1[2]), int(roi1[3]))

	## STEP 11: Keep objects that overlap with the ROI
	roi_objects, roi_obj_hierarchy, kept_mask, obj_area = pcv.roi_objects(img1, roi_contour, roi_hierarchy, id_objects, obj_hierarchy, roi_type='partial')
	
    ## NEW STEP: automatically crop an image to a contour 
	cropped = pcv.auto_crop(img1, np.vstack(id_objects), padding_x=0, padding_y=0, color='black')

	## STEP 12: This function take a image with multiple contours and clusters them based on user input of rows and columns
	clusters_i, contours, hierarchies = pcv.cluster_contours(img1, roi_objects, roi_obj_hierarchy, int(roi1[4]), int(roi1[5]), show_grid=True)
	
	## STEP 13: This function takes the outputs from cluster_contours and plots on it the labels of each object. This helps for better visualisation of what objects declustered to what output
	pcv.params.debug = "print"
	clusterfuck = pcv.visualize.clustered_contours(img1, clusters_i, roi_objects, roi_obj_hierarchy, nrow=int(roi1[4]),ncol=int(roi1[5]))
	listedpic=glob.glob(os.path.join(args.outdir, "*_clusters.png"))
	latest_file=max(listedpic, key=os.path.getctime)
	shutil.copy(latest_file, os.path.join(os.path.dirname(args.outdir), "{0}_clusters_contours.jpg".format(tray)))
	#shutil.copy(os.path.join(args.outdir, str(sorted(os.listdir(args.outdir))[-1])), os.path.join(os.path.dirname(args.outdir), "{0}_clusters_contours.jpg".format(tray)))
	pcv.params.debug = args.debug
	
	## STEP 13: This function takes clustered contours and splits them into multiple images
	out = "./"
	names = args.names
	output_path = pcv.cluster_contour_splitimg(img1, clusters_i, contours, hierarchies, out, file=filename, filenames=names)

#Get the necessary parameters
args=options()
pcwd=os.getcwd()
ima=args.input
imb=ima.rsplit("/",1)[-1]
picname=ima.rsplit("/",1)[-1][:-4]

#And run the declustering on the file
main(ima, picname, args, pcwd)

EOF

#Generate second python script for output image measurements
cat <<\EOF > "$mpwd"/$2_python2.py
#!/usr/bin/python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys, traceback
import os
import re
import numpy as np
import argparse
import string
import glob
import shutil
import matplotlib
from PIL import Image
from plantcv import plantcv as pcv

### Parse command-line arguments
def options():
	parser = argparse.ArgumentParser(description="Imaging processing with opencv")
	parser.add_argument("-i", "--input", help="Input image directory.", required=True)
	parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=False)
	parser.add_argument("-n", "--names", help="path to txt file with names of genotypes to split images into", required =False)
	parser.add_argument("-r","--result", help="result file.", required= True )
	parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action=None)
	parser.add_argument("-v", "--variables",help="List of variables to be passed to pipeline", required=True)
	args = parser.parse_args()
	return args

#Functuin to analyze the rosette area
def analysis(image1, binary1, plantnr, tray):
	args = options()
	
	#Get variable info from variables arguments
	with open(args.variables) as f:
		roi1=[line.split(' ')[0] for line in f]

	# Convert images to np arrays for analysis, something I have added
	binarpic=Image.open(binary1)
	imagpic=Image.open(image1)
	binar=np.array(binarpic)
	ima=np.array(imagpic)

	# Find objects again within ROI, see also in main function!, from PlantCV
	id_objects, obj_hierarchy = pcv.find_objects(ima, binar)
	roi_contour, roi_hierarchy = pcv.roi.rectangle(ima, int(roi1[0]), int(roi1[1]), int(roi1[2]), int(roi1[3]))

	#From image inputs, generate np arrays and calculate shape area, generate error if analyze_object could fail due to vertices errors, from PlantCV
	roi_objects, hierarchy3, kept_mask, obj_area = pcv.roi_objects(ima, roi_contour, roi_hierarchy, id_objects, obj_hierarchy, roi_type='partial')
	obj, mask = pcv.object_composition(ima, roi_objects, hierarchy3)
	try:
		shape_image = pcv.analyze_object(ima, obj, mask)
		shape_area = pcv.outputs.observations['area']['value']
		data=' '.join((tray, str(plantnr), str(shape_area)))
		out=open(args.result, "a+")
		out.write((data + '\n'))
		out.close
	except:
		print("%s has not enough vertices, omitting from data.." %(str(image1)))
		data=' '.join((tray,str(plantnr),"vertix_err"))
		out=open(args.result, "a+")
		out.write((data + '\n')) 
		out.close

#Get necessary parameters
args=options()
imga=args.input
imgb=imga[:-9] + ".png"
picname=imga.rsplit('/')[-1].rsplit('_')[-3]
trayname=imga.rsplit('/')[-2]

#And analyze each picture
analysis(imgb, imga, picname, trayname)

EOF

#And third python script to summarize all data in a figure
cat <<\EOF > "$mpwd"/$2_python3.py
#Post analysis script, to analyze results of plantcv output and to generate figures 
#NOTES file should contain: 
#    - Tab-separated picture name, genotype name, each on a new line. The order of genotypes will used as key to the order in the final graph
#    - A line describing the Pixel ratio as "Pixel ratio:00.00", the proogram will take any number after the last colon in the line containing Pixel ratio.

import os
import pandas as pd
import argparse
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import textwrap
import scipy.stats as stats
from statsmodels.stats.multicomp import (pairwise_tukeyhsd, MultiComparison)
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib

try:
    __IPYTHON__
except NameError:
    mwd=os.getcwd()
    print("Running from command line")
else:
    mwd="C:/Users/3630188/OneDrive - Universiteit Utrecht/dmr6dlo1 EMS/20201119 Complementation assay" #If using in Spyder/IPython
    print("Running from Spyder/IPython")

parser = argparse.ArgumentParser(description="Imaging postprocessing")
parser.add_argument("-o", "--output", help="output folder name", required=True)
args = parser.parse_args()
inres=os.path.join(mwd, "{0}.totalresults.txt".format(args.output))
notes=os.path.join(mwd, "notes.txt")

#Dependency: pairwise comparison letter generator
"""
'Python pairwise comparison letter generator'

Github: PhilPlantMan
"""
import string
import random
import statsmodels.stats.multicomp as multi
import scikit_posthocs as sp


def multi_comparisons_letter_df_generator(comparisons_df, letter_ordering_series = None, 
                                          primary_optimisation_parameter = "Number of different letters", 
                                          monte_carlo_cycles = 5, letter_separator = '', ): 
    """
    Function takes a df listing pairwise comparisons with a cols labelled 'group1' and 'group2' for the two groups being compared 
    and another column labelled 'reject' with boolean values corresponding to whether the null hypothesis should be rejected 
    i.e. True: Both treatments are significantly different
    
    letter_ordering_series (default = None): In theory, which letters are assigned to each non-significance grouping is
    arbitrary and therefor the order can be changed. Offering letter_ordering_series a series with the same index as the output
    will make sure that the order that letters are assigned will follow letter_ordering_series from max to min. For boxplots,
    letter_ordering_series with median values is a good choice.
    
    monte_carlo_cycles (default = 5): Function will always return correct letter representation however it may be suboptimal. 
    Within each monte carlo cycle, a random letter is deleted until the representation breaks. The result with the optimum
    number layout of letters after n monte_carlo_cycles is returned. 
    
    The optimum letter layout is set by primary_optimisation_parameter (default = "Number of different letters"):
        'Number of different letter' optimises for fewest different letters
        "Min letters per row" optimises for the fewest letters assigned per treatment
        "Letter total" optimises for the fewest total letters of the treatments combined
        
    letter_separator (default = ''): Separator for each letter in string assigned to each treatment
    
    
    Letter representation is determined by the method described by Piepho 2004: An Algorithm for a Letter-Based Representation
    of All-Pairwise Comparisons
    """
    #'Insert' stage
    #make df with all unique groups as index
    letters_df = comparisons_df['group1'].append(comparisons_df['group2']).drop_duplicates().to_frame().set_index(0)
    
    letters_df[letters_df.shape[1]] = 1
    for pos_result in comparisons_df.loc[comparisons_df['reject']==True].index:
        group1 = comparisons_df.loc[pos_result, 'group1']
        group2 = comparisons_df.loc[pos_result, 'group2']
        for letter_col in letters_df:
            group1_val = letters_df.loc[group1,letter_col]
            group2_val = letters_df.loc[group2,letter_col]
            if group1_val == 1 and group2_val == 1:
                #duplicate column
                new_col = letters_df.shape[1]
                letters_df[new_col] = letters_df[letter_col]
                #del val at group1 first col and at group2 new col
                letters_df.loc[group1,letter_col] = 0
                letters_df.loc[group2,new_col] = 0
    #'Absorb' stage          
    for col in letters_df:
       other_cols_list = list(letters_df)
       other_cols_list.remove(col)
       col_total = letters_df[col].sum()
       for other_col in other_cols_list:
           matched_total = 0
           for row in letters_df.index:
               if letters_df.loc[row, col] == 1 and letters_df.loc[row, other_col]: matched_total +=1
           if col_total == matched_total:
               letters_df.drop(col, axis = 1, inplace = True)  
               break
        
    def check_letters_against_tests(test_df, letters_df):
        if letters_df.sum(axis = 1).min() == 0: return False
        for result_row in test_df.index:
            group1 = test_df.loc[result_row, 'group1']
            group2 = test_df.loc[result_row, 'group2']
            reject = bool(test_df.loc[result_row, 'reject'])
            count_of_true_trues = 0
            count_of_true_falses = 0
            for letter_col in letters_df:
                group1_val = letters_df.loc[group1,letter_col]
                group2_val = letters_df.loc[group2,letter_col]
                if reject:
                    if group1_val != group2_val: count_of_true_trues += 1
                    if group1_val == 1 and group2_val == 1: 
                        return False
                if reject == False:
                    if group1_val == 1 and group2_val == 1: count_of_true_falses += 1
            if reject and count_of_true_trues == 0: 
                return False
            if reject == False and count_of_true_falses == 0: 
                return False
        return True

    #'Sweep stage' with monte carlo optimisation
    for i in range(monte_carlo_cycles):
        num_of_letters = letters_df.sum().sum()
        num_list = list(np.arange(start = 1, stop = 1+ num_of_letters))
        letters_df_monte_order = letters_df.copy()
        for row in letters_df_monte_order.index:
            for col in letters_df_monte_order:
                if letters_df_monte_order.loc[row,col] == 0: continue
                random_num = random.sample(num_list, 1)[0]
                letters_df_monte_order.loc[row,col] = random_num
                num_list.remove(random_num)
        
        current_letters_df = letters_df.copy()
        for pos in range(num_of_letters + 1):     
            mask = letters_df_monte_order.isin([pos])
            zero_df = letters_df.copy().loc[:] = 0
            letters_df_copy = current_letters_df.copy()
            letters_df_copy.mask(mask, other = zero_df, inplace = True)
            if check_letters_against_tests(comparisons_df,letters_df_copy):
                current_letters_df = letters_df_copy
        
        for col in letters_df:
            if current_letters_df[col].sum() == 0: current_letters_df.drop(col, axis = 1, inplace = True)
            
        # determine fitness parameters for optimisation
        current_fitness_parameter_vals = {"Min letters per row":current_letters_df.sum(axis = 1).max(),
                                          "Number of different letters": current_letters_df.shape[1],
                                          "Letter total": current_letters_df.sum().sum()}
        if i == 0: 
            best_fitness_parameter_vals = current_fitness_parameter_vals
            best_letters_df = current_letters_df
            continue
        
        if current_fitness_parameter_vals[primary_optimisation_parameter] > best_fitness_parameter_vals[primary_optimisation_parameter]:
            continue
        if current_fitness_parameter_vals[primary_optimisation_parameter] < best_fitness_parameter_vals[primary_optimisation_parameter]:
            best_letters_df = current_letters_df.copy()
            best_fitness_parameter_vals = current_fitness_parameter_vals
            
        if sum(current_fitness_parameter_vals.values()) < sum(best_fitness_parameter_vals.values()):
            best_letters_df = current_letters_df.copy()
            best_fitness_parameter_vals = current_fitness_parameter_vals
    
    #order cols
    if isinstance(letter_ordering_series, pd.Series):
        scoring_df = pd.DataFrame(index = best_letters_df.index)
        for row in best_letters_df.index:
            for col in best_letters_df:
                scoring_df.loc[row, col] = best_letters_df.loc[row, col] * letter_ordering_series[row]
        scoring_df = scoring_df.replace(0, np.NaN)
        scoring_means = scoring_df.mean(axis = 0).sort_values(ascending = False)
        best_letters_df = best_letters_df[scoring_means.index]
    # letter the cols     
    for col_name, col_num in zip(best_letters_df, range(len(best_letters_df.columns))):
        letter = string.ascii_lowercase[col_num]
        best_letters_df.loc[best_letters_df[col_name] == 1, col_name] = letter
    # make df with strings ready for presentation
    best_string_df = pd.DataFrame(index = best_letters_df.index)
    best_string_df.loc[:,'string'] = ""
    for row in best_letters_df.index:
        for col in best_letters_df:
            if best_letters_df.loc[row, col] != 0:
                letter_string = best_string_df.loc[row, 'string']
                letter = best_letters_df.loc[row, col]
                if letter_string == "": best_string_df.loc[row, 'string'] = letter
                else: best_string_df.loc[row, 'string'] = letter_separator.join((letter_string, letter))
                
    return best_string_df

print("Environment loaded. Processing data...")

#Now load the input as a pandas file and the notes file up until the first blank line and then input to dict
dictio=dict(zip(["pic", "plantnr", "SqPix"],['str','int','float']))
indf=pd.read_csv(inres, sep=" ", header=None, names=["pic", "plantnr", "SqPix"])
indf=indf[indf['SqPix'].apply(np.isreal)]
indf=indf.astype(dictio)
indf['plantnr']=indf.index

picnames=list(indf.pic.unique())
catdict={}
catorder=[]
for line in open(notes, 'r').readlines():
    if any(name in line for name in picnames):
        cat=line.strip().split("\t")[-1]
        name=line.split("\t")[0].split(".")[0]
        if cat not in catorder:
            catorder.append(cat)
            catdict[name] = cat
    if "Pixel ratio" in line:
        pixratio=float(line.strip().split(":")[-1])

# catorder=sorted(set(catorder), key=catorder.index)

indf['cat'] = indf.pic.replace(catdict)
indf['Sqmm']=indf['SqPix']/pixratio

#General statistics
try:
	shap=stats.shapiro(indf.groupby(by="cat")['Sqmm'].agg(stats.sem))
except ValueError as err:
	shap="No shapiro test performed: {0}".format(err)
f=ols("Sqmm ~ cat", data = indf).fit()
aov=sm.stats.anova_lm(f, typ=2)
mc=MultiComparison(indf['Sqmm'], indf['cat'], group_order=catorder)
thsd=mc.tukeyhsd()
thdf = pd.DataFrame(data=thsd._results_table.data[1:], columns=thsd._results_table.data[0])
letterdf=multi_comparisons_letter_df_generator(thdf)

#Now plot the dataseries per column
sb.set_style("darkgrid")
matplotlib.rcParams['ps.fonttype'] = 42 #Needed for Illustrator-editable files
matplotlib.rcParams['pdf.fonttype'] = 42 #Needed for Illustrator-editable files
sb.set(font_scale=1.2)
fig, ax = plt.subplots(figsize=(1+0.75*len(catorder),7))
bx=sb.boxplot(
    data=indf.reset_index(),
    x="cat",
    y="Sqmm",
    order=catorder,
    showfliers=False)
sw=sb.swarmplot(
    data=indf,
    x="cat",
    y="Sqmm",
    order=catorder,
    linewidth=1,
    color='.2',
    size=4)

maxes=indf.groupby(['cat'])['Sqmm'].max()
soneg=dict(zip(range(0, len(catorder),1), catorder))
pos=range(len(letterdf))

ydif=np.abs(ax.get_yticks()[0]-ax.get_yticks()[1])
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]+ydif)

for tick, label in zip(pos, ax.get_xticklabels()):
    ax.text(pos[tick], maxes[soneg.get(tick)] + 0.4*ydif,letterdf.iloc[tick]['string'], horizontalalignment="center", color="black", size="medium", weight="semibold")

ax.set_xticks(ax.get_xticks())
if any("mathregular" in s for s in catorder):
    ax.set_xticklabels([textwrap.fill(i, 30) for i in catorder],rotation=45, ha="right")
else:
    ax.set_xticklabels([textwrap.fill(i, 15) for i in catorder],rotation=45, ha="right")
ax.set_ylabel("Rosette area ($\mathregular{mm^2}$)")
ax.set_xlabel("")
plt.tight_layout()
plt.savefig(os.path.join(mwd, "{0}_boxplot.png".format(args.output)))
plt.savefig(os.path.join(mwd, "{0}_boxplot.svg".format(args.output)))

#And save the statistics file
with open(os.path.join(mwd, "{0}_statistics.txt".format(args.output)), "w") as text_file:
    text_file.write("""Shapiro-Wilk on SEM:
{0}\n\n\n
ANOVA:
{1}\n\n\n
Post-Hoc:
{2}""".format(str(shap), str(aov), str(thsd)))

print("All done. boxplot.png/.svg and statistics.txt saved to {0}".format(mwd))

EOF

#First bash script that calls Python1 and Python2 per input picture
cat <<\EOF > "$mpwd"/$2_bash1.sh
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=t.vanbutselaar@uu.nl
#SBATCH --job-name=decluster

#Prep environment
mpwd="$3"
file="$4"
cd "$mpwd"
picname=$(basename "$file")
picname="${picname%.*}"
mkdir "$mpwd"/$2/$picname 
cd "$mpwd"/$2/$picname

#Decluster the image
echo "Declustering image " $file #A shout-out to stdout what image it's currently working on, so you have an idea where in your set-up it is.
python3 "$mpwd"/$2_python1.py -i "$file" -o "$mpwd"/$2/$picname -D print -v "$mpwd"/$2_variables.txt #Start the declustering of multi-plant pictures. The -D argument will ensure that also intermediary pictures of each processing step are saved. It is important to inspect these intermediary pictures after the pipeline is finished to get an idea of how the declustering behaved, and if parameters need tweaking.

#And per declustered output file, measure the leaf area
for file2 in "$mpwd/$2/$picname/$picname"*_mask.png
do
	python3 "$mpwd"/$2_python2.py -i "$file2" -D print -r "$mpwd"/$2/$picname.results.txt -v "$mpwd"/$2_variables.txt #Now analyze every single (supposed) plant from the declustered output
done

EOF

#Second bash script that summarizes output and calls Python3
cat <<\EOF > "$mpwd"/$2_bash2.sh
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=t.vanbutselaar@uu.nl
#SBATCH --job-name=finalize

#Prep environment
mpwd="$3"

#Concatenate all single image results into one major results file, containing per line, the original picture name, the plant number, and the size in squared pixels
cat "$mpwd"/$2/*.results.txt | sort > "$mpwd"/$2.totalresults.txt

#Now do post-analysis summarizing of all data into a boxplot and a statistics file
python3 "$mpwd/$2_python3.py" -o $2
 
EOF

#Check if user wants to run the script in parallel on the HPC or sequentially on local PC
printf "\n"
testcmd () {
    command -v "$1" >/dev/null
}
if testcmd sbatch; then
	echo "sbatch program detected, starting HPC analysis"
	hpc=true
else
	echo "sbatch program NOT detected, starting local analysis"
	hpc=false
fi

sleep 2
printf "\nStarting image analysis script using PlantCV.\n\nBe sure to check out PlantCV main workflows for the image declustering and image analysis pipelines on https://plantcv.readthedocs.io/en/stable/multi-plant_tutorial/ and https://plantcv.readthedocs.io/en/stable/vis_tutorial/. On these sites, the authors of the PlantCV explain all processing steps and what each output image is from.\n\n"
sleep 2

# Prep the environment and perform necessary checks
#Check that the notes.txt file is present, terminate script if not running
set -e
if ! [ -f "$mpwd/notes.txt" ]; then
	echo "notes.txt does not exist, terminating..."
	exit 1
fi
#Check that input folder exists and contains files
if [ -n "$(find "$mpwd/$1" -maxdepth 0 -type d -empty 2>/dev/null)" ]; then
	echo "input folder does not exist or is empty, terminating..."
	exit 1
fi

if ! [ -d "$mpwd/$2" ]; then
	mkdir "$mpwd"/$2
fi

#And now call all scripts per input picture and submit jobs
if $hpc; then
	for file in "$mpwd"/$1/*; do #Call all files that are in the specified input folder
		sbatch "$mpwd"/$2_bash1.sh $1 $2 "$mpwd" "$file"
		sleep 0.5
	done
	sbatch --dependency=singleton --job-name=decluster "$mpwd"/$2_bash2.sh $1 $2 "$mpwd"
	watch squeue -u $USER
else
	for file in "$mpwd"/$1/*; do #Call all files that are in the specified input folder
		bash "$mpwd"/$2_bash1.sh $1 $2 "$mpwd" "$file"
	done
	bash "$mpwd"/$2_bash2.sh $1 $2 "$mpwd"
fi
