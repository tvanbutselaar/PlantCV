###
# @author: tvanbutselaar
# Current version: 20211123
# Changes: Aesthetic changes to plot, bugfix to catdict/catorder

# This version compatible with PlantCV v3.12.1 (and possibly later)

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
			# Protips: 
			# flank your text with $ if you want the text to be printed in italics in the plot
			# use mathregular-based expressions if you want to use subscript, superscript, or special characters in the plot: https://matplotlib.org/3.3.3/tutorials/text/mathtext.html
	# all variables for this script to run separated on newlines


# Before use, update the notes.txt file to optimized parameters for your batch of pictures. Use a representative picture and analyze in Photoshop to estimate the parameters in the variable list. 
###

mpwd=$(pwd)



#Generate main python script for declustering pictures
cat <<\EOF > "$mpwd"/"$2"_python1.py
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
	# if args.names!= None:
		# args.names="../../names/%s.txt" %(str(tray))
		# print("Names file: %s" %(args.names))
	# else:
		# print("No names file submitted for %s" %(str(tray)))
		
	#Get variables from variable arguments
	for line in open(args.variables, 'r').readlines():
		if line.startswith("shift_up"):
			shift_up=int(line.strip().split(":")[-1])
		if line.startswith("shift_left"):
			shift_left=int(line.strip().split(":")[-1])
		if line.startswith("ROI_x_adj"):
			ROI_x_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_y_adj"):
			ROI_y_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_h_adj"):
			ROI_h_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_w_adj"):
			ROI_w_adj=int(line.strip().split(":")[-1])
		if line.startswith("rotation_deg"):
			rotation_deg=float(line.strip().split(":")[-1])
		if line.startswith("cluster_nrow"):
			cluster_nrow=int(line.strip().split(":")[-1])
		if line.startswith("cluster_ncol"):
			cluster_ncol=int(line.strip().split(":")[-1])
		if line.startswith("bin_tresh"):
			bin_tresh=int(line.strip().split(":")[-1])
		if line.startswith("fill_size"):
			fill_size=int(line.strip().split(":")[-1])
		if line.startswith("dilation_kernel"):
			dilation_kernel=int(line.strip().split(":")[-1])
		if line.startswith("dilation_iter"):
			dilation_iter=int(line.strip().split(":")[-1])
		
	# Read image
	img, path, filename = pcv.readimage(img8)
	pcv.params.debug=args.debug #set debug mode

	## STEP 2: Normalize the white color so you can later compare color between images. Inputs
	img1 = pcv.white_balance(img,roi=None)
	
	## STEP 3: Rotate the image
	if rotation_deg != 0:
		rotate_img = pcv.transform.rotate(img1, float(rotation_deg), True)
	else:
		rotate_img = img1
	
	## STEP 4: Shift image. This step is important for clustering later on. The resulting image is the same size as the original.
	if shift_up > 0:
		shift1=pcv.shift_img(rotate_img, shift_up, "bottom")
	else:
		shift1=rotate_img
	if shift_left > 0:
		shift2=pcv.shift_img(shift1, shift_left, "right")
		img1=shift2
	else:
		img1=shift1

	## STEP 5: Convert image from RGB colorspace to LAB colorspace. Keep only the green-magenta channel (grayscale).
	a = pcv.rgb2gray_lab(img1, "a")

	## STEP 6: Set a binary threshold on the saturation channel image.
	img_binary = pcv.threshold.binary(a, bin_tresh, 255, "dark")

	## STEP 7: Fill in small objects (speckles)
	fill_image = pcv.fill(img_binary, fill_size)

	## STEP 8: Dilate so that you don't lose leaves (just in case). Only do this if dilation parameters allow this.
	if dilation_kernel > 1 and dilation_iter > 1:
		dilated = pcv.dilate(fill_image, dilation_kernel, dilation_iter)
	else:
		print("Dilation kernal and iter <= 1, skipping dilation...")
		dilated=fill_image

	## STEP 9: Find objects (contours: black-white boundaries)
	id_objects, obj_hierarchy = pcv.find_objects(img1, dilated)

	## STEP 10: Define region of interest (ROI)
	roi_contour, roi_hierarchy = pcv.roi.rectangle(img1, ROI_x_adj, ROI_y_adj, ROI_h_adj, ROI_w_adj)

	## STEP 11: Keep objects that overlap with the ROI
	roi_objects, roi_obj_hierarchy, kept_mask, obj_area = pcv.roi_objects(img1, roi_contour, roi_hierarchy, id_objects, obj_hierarchy, roi_type='partial')
	
    ## NEW STEP: automatically crop an image to a contour 
	cropped = pcv.auto_crop(img1, np.vstack(id_objects), padding_x=0, padding_y=0, color='black')

	## STEP 12: This function take a image with multiple contours and clusters them based on user input of rows and columns
	clusters_i, contours, hierarchies = pcv.cluster_contours(img1, roi_objects, roi_obj_hierarchy, cluster_nrow, cluster_ncol, show_grid=True)
	
	## STEP 13: This function takes the outputs from cluster_contours and plots on it the labels of each object. This helps for better visualisation of what objects declustered to what output
	pcv.params.debug = "print"
	clusterfuck = pcv.visualize.clustered_contours(img1, clusters_i, roi_objects, roi_obj_hierarchy, nrow=cluster_nrow,ncol=cluster_ncol)
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
cat <<\EOF > "$mpwd"/"$2"_python2.py
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
	for line in open(args.variables, 'r').readlines():
		if line.startswith("ROI_x_adj"):
			ROI_x_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_y_adj"):
			ROI_y_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_h_adj"):
			ROI_h_adj=int(line.strip().split(":")[-1])
		if line.startswith("ROI_w_adj"):
			ROI_w_adj=int(line.strip().split(":")[-1])

	# Convert images to np arrays for analysis, something I have added
	binarpic=Image.open(binary1)
	imagpic=Image.open(image1)
	binar=np.array(binarpic)
	ima=np.array(imagpic)

	# Find objects again within ROI, see also in main function!, from PlantCV
	id_objects, obj_hierarchy = pcv.find_objects(ima, binar)
	roi_contour, roi_hierarchy = pcv.roi.rectangle(ima, ROI_x_adj, ROI_y_adj, ROI_h_adj, ROI_w_adj)

	#From image inputs, generate np arrays and calculate shape area, generate error if analyze_object could fail due to vertices errors, from PlantCV
	roi_objects, hierarchy3, kept_mask, obj_area = pcv.roi_objects(ima, roi_contour, roi_hierarchy, id_objects, obj_hierarchy, roi_type='partial')
	obj, mask = pcv.object_composition(ima, roi_objects, hierarchy3)
	try:
		shape_image = pcv.analyze_object(ima, obj, mask)
		shape_area = pcv.outputs.observations['default']['area']['value']
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
cat <<\EOF > "$mpwd"/"$2"_python3.py
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

def pletters(thsd):
    #this is a function to do Piepho method.  AN Alogrithm for a letter based representation of al-pairwise comparisons.  
    tot=len(thsd.groupsunique)
    #make an empty dataframe that is a square matrix of size of the groups. #set first column to 1
    df_ltr=pd.DataFrame(np.nan, index=np.arange(tot),columns=np.arange(tot))
    df_ltr.iloc[:,0]=1
    count=0
    df_nms = pd.DataFrame('', index=np.arange(tot), columns=['names'])  # I make a dummy dataframe to put axis labels into.  sd stands for signifcant difference
    
    for i in np.arange(tot):   #I loop through and make all pairwise comparisons. 
        for j in np.arange(i+1,tot):
            #print('i=',i,'j=',j,thsd.reject[count])
            if thsd.reject[count]==True:
                for cn in np.arange(tot):
                    if df_ltr.iloc[i,cn]==1 and df_ltr.iloc[j,cn]==1: #If the column contains both i and j shift and duplicat
                        df_ltr=pd.concat([df_ltr.iloc[:,:cn+1],df_ltr.iloc[:,cn+1:].T.shift().T],axis=1)
                        df_ltr.iloc[:,cn+1]=df_ltr.iloc[:,cn]
                        df_ltr.iloc[i,cn]=0
                        df_ltr.iloc[j,cn+1]=0
                    #Now we need to check all columns for abosortpion.
                    for cleft in np.arange(len(df_ltr.columns)-1):
                        for cright in np.arange(cleft+1,len(df_ltr.columns)):
                            if (df_ltr.iloc[:,cleft].isna()).all()==False and (df_ltr.iloc[:,cright].isna()).all()==False: 
                                if (df_ltr.iloc[:,cleft]>=df_ltr.iloc[:,cright]).all()==True:  
                                    df_ltr.iloc[:,cright]=0
                                    df_ltr=pd.concat([df_ltr.iloc[:,:cright],df_ltr.iloc[:,cright:].T.shift(-1).T],axis=1)
                                if (df_ltr.iloc[:,cleft]<=df_ltr.iloc[:,cright]).all()==True:
                                    df_ltr.iloc[:,cleft]=0
                                    df_ltr=pd.concat([df_ltr.iloc[:,:cleft],df_ltr.iloc[:,cleft:].T.shift(-1).T],axis=1)
    
            count+=1
    
    #I sort so that the first column becomes A        
    df_ltr=df_ltr.sort_values(by=list(df_ltr.columns),axis=1,ascending=False)
    
    # I assign letters to each column
    for cn in np.arange(len(df_ltr.columns)):
        df_ltr.iloc[:,cn]=df_ltr.iloc[:,cn].replace(1,chr(97+cn)) 
        df_ltr.iloc[:,cn]=df_ltr.iloc[:,cn].replace(0,'')
        df_ltr.iloc[:,cn]=df_ltr.iloc[:,cn].replace(np.nan,'') 
    
    #I put all the letters into one string
    df_ltr=df_ltr.astype(str)
    #print(df_ltr)
    #print('\n')
    #print(df_ltr.sum(axis=1))
    return df_ltr.sum(axis=1)

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
        catdict[name] = cat
        if cat not in catorder:
            catorder.append(cat)
    if line.startswith("pixel_ratio"):
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
letterdf=pletters(thsd)

#Now plot the dataseries per column
plt.style.use("default")
# sb.set_style("whitegrid")
matplotlib.rcParams['ps.fonttype'] = 42 #Needed for Illustrator-editable files
matplotlib.rcParams['pdf.fonttype'] = 42 #Needed for Illustrator-editable files
plt.rc('font', size=15)
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.axisbelow'] = True
matplotlib.rcParams['axes.grid.axis'] = 'both'

fig, ax = plt.subplots(figsize=(1+0.75*len(catorder),7))
bx=sb.boxplot(
    data=indf.reset_index(),
    x="cat",
    y="Sqmm",
    order=catorder,
    palette="Dark2",
    showfliers=False)
sw=sb.swarmplot(
    data=indf,
    x="cat",
    y="Sqmm",
    order=catorder,
    linewidth=1,
    palette="Dark2",
    alpha=0.4,
    size=4)
sb.despine()

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

##And the dataframes
indf.to_csv(os.path.join(mwd, "{0}_df.tsv".format(args.output)), sep="\t")
indf.to_pickle(os.path.join(mwd, "{0}_df.pkl".format(args.output)))

print("All done. {1}_ boxplot.png/.svg, statistics.txt, dataframe.tsv/.pkl saved to {0}".format(mwd, args.output))

EOF

#First bash script that calls Python1 and Python2 per input picture
cat <<\EOF > "$mpwd"/"$2"_bash1.sh
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
python3 "$mpwd"/$2_python1.py -i "$file" -o "$mpwd"/$2/$picname -D print -v "$mpwd"/notes.txt #Start the declustering of multi-plant pictures. The -D argument will ensure that also intermediary pictures of each processing step are saved. It is important to inspect these intermediary pictures after the pipeline is finished to get an idea of how the declustering behaved, and if parameters need tweaking.

#And per declustered output file, measure the leaf area
for file2 in "$mpwd/$2/$picname/$picname"*_mask.png
do
	python3 "$mpwd"/$2_python2.py -i "$file2" -D print -r "$mpwd"/$2/$picname.results.txt -v "$mpwd"/notes.txt #Now analyze every single (supposed) plant from the declustered output
done

EOF

#Second bash script that summarizes output and calls Python3
cat <<\EOF > "$mpwd"/"$2"_bash2.sh
#!/bin/bash
#SBATCH --time=00:30:00
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
	mkdir "$mpwd"/"$2"
fi

#And now call all scripts per input picture and submit jobs
if $hpc; then
	for file in "$mpwd"/"$1"/*; do #Call all files that are in the specified input folder
		sbatch "$mpwd"/"$2"_bash1.sh "$1" "$2" "$mpwd" "$file"
		sleep 0.5
	done
	sbatch --dependency=singleton --job-name=decluster "$mpwd"/"$2"_bash2.sh "$1" "$2" "$mpwd"
	watch squeue -u $USER
else
	for file in "$mpwd"/"$1"/*; do #Call all files that are in the specified input folder
		bash "$mpwd"/"$2"_bash1.sh "$1" "$2" "$mpwd" "$file"
	done
	bash "$mpwd"/"$2"_bash2.sh "$1" "$2" "$mpwd"
fi

rm "$mpwd"/"$2"_bash*.sh
rm "$mpwd"/"$2"_python*.sh