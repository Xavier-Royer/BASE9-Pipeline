#!/usr/bin/python3

# this code will take as an input from the user a .res file 
#   and a .phot file. It will then run makeCMD for the 
#   average values from the .res file, and overplot the isochrone 
#   on the user-specified filters for a CMD.
#
#   USAGE (example):
#          calc_plot_iso.py myPhot.phot V B-V file.res 

# check the CMD for a BASE-9 phot file

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# read the files in
#photfile = sys.argv[1]
#mag     = sys.argv[2]
#color     = sys.argv[3].split('-')
#resfile   = sys.argv[4]

#photfile = '../hyades_UBVRIKs_noGiants.phot'
#mag = 'V'
#color = 'V-Ks'.split('-')
#resfile = 'hyades_UBVRIKs_noGiants_mistRot_all.res'


def create_iso_plot(mag, color,clusterName, version = ""):
    root_dir = clusterName + "/"
    color     = color.split('-')
    resfile = root_dir + clusterName + version +".res"
    photfile = root_dir + clusterName + version +".phot"
    model_file = root_dir + clusterName + version + ".ms"
    #model_file = rootname+'.ms'

    ########################################
    #### Read in .res file,    #############
    #### calculate statistics. #############
    ########################################

    # load the data into arrays
    data    = np.loadtxt(resfile,skiprows=1,usecols=[0,1,2,3,5])
    age_all = data[:,0]
    feh_all = data[:,1]
    dis_all = data[:,2]
    av_all  = data[:,3]
    stage   = data[:,4]

    # only read from main run (don't plot burn in)
    ind    = np.where(stage == 3)
    age    = age_all[ind]
    feh    = feh_all[ind]
    dis    = dis_all[ind]
    av     = av_all[ind]
    i      = list(range(0,len(age),1)) # iteration number

    # calculate stats
    avg_age = np.average(age)     ; std_age = np.std(age)
    avg_age_l = np.average(10**age)/1e9     ; std_age_l = np.std(10**age)/1e9
    avg_feh = np.average(feh)     ; std_feh = np.std(feh)
    avg_dis = np.average(dis)     ; std_dis = np.std(dis)
    avg_av  = np.average(av)      ; std_av  = np.std(av)

    # make plot labels, mean +/- stdev
    a = '%.3f' % avg_age ; s = '%.4f' % std_age ; age_lab = a+' +/- '+s
    a = '%.3f' % avg_age_l ; s = '%.4f' % std_age_l ; age_lab_l = a+' +/- '+s
    a = '%.3f' % avg_feh ; s = '%.4f' % std_feh ; feh_lab = a+' +/- '+s
    a = '%.3f' % avg_dis ; s = '%.4f' % std_dis ; dis_lab = a+' +/- '+s
    #avg_av1 = avg_av * 1000 ; std_av1 = std_av * 1000
    #a = '%.2f' % avg_av1  ; s = '%.2f' % std_av1  ; av_lab  = a+' +/- '+s
    a = '%.4f' % avg_av  ; s = '%.4f' % std_av  ; av_lab  = a+' +/- '+s

    ########################################
    #### Run makeCMD using mean ############
    #### values of .res file.   ############
    ########################################

    #cmd = 'makeCMD --outputFileBase '+'"'+rootname+'" '+'--priorLogAge '+str(avg_age)+' --startingLogAge '+str(avg_age)+' --priorFe_H '+str(avg_feh)+' --startingFe_H '+str(avg_feh)+' --priorDistMod '+str(avg_dis)+' --startingDistMod '+str(avg_dis)+' --sigmaDistMod '+str(std_dis)+' --priorAv '+str(avg_av)+' --startingAv '+str(avg_av)

    # run makeCMD; delete .wd file because I don't need it.
    #os.system(cmd)
    #os.system('rm '+rootname+'.wd')

    ########################################
    #### Thin makeCMD output    ############
    ########################################

    #cmd = 'thin_makeCMD.py '+rootname+'.ms'
    make_thinCMD(model_file)
    model_file += ".thin"

    # run thin_makeCMD
    #os.system(cmd)

    # figure out which columns from the model file
    # are the filters the user wants to plot
    header = np.loadtxt(model_file,max_rows=1,dtype=str) #+'.thin'

    yplot     = np.where(header == mag)[0][0]
    color1    = np.where(header == color[0])[0][0]
    color2    = np.where(header == color[1])[0][0]

    # read in model file, set variables, use only user-specified mags
    data = np.loadtxt(model_file,skiprows=1)
    y_model   = data[:,yplot]
    color_mod = data[:,color1]-data[:,color2]

    ########################################
    #### Read in .res file,    #############
    #### calculate statistics. #############
    ########################################
    # below is unchanged from the regular check_cmd.py file; it will plot the photometry
    # only difference is that I've added code to overplot the newly made isochrone and 
    # display the average values on the CMD.

    header = np.loadtxt(photfile,max_rows=1,dtype=str)

    yplot     = np.where(header == mag)[0][0]
    color1    = np.where(header == color[0])[0][0]
    color2    = np.where(header == color[1])[0][0]

    ylabel = header[int(yplot)]
    xlabel = header[int(color1)]+' - '+header[int(color2)]

    data = np.loadtxt(photfile,skiprows=1)
    nstars = len(data[:,0])
    nfilts = int((data.shape[1] - 6)/2)


    # include only stars for which phot in all 3 filters exists
    good   = np.where((data[:,yplot+nfilts] > 0) & (data[:,color1+nfilts] > 0) & (data[:,color2+nfilts] > 0))

    ymag    = data[:,yplot][good]
    color   = data[:,color1][good] - data[:,color2][good]
    mag_err = data[:,yplot+nfilts][good]
    col_err = np.sqrt((data[:,color1+nfilts][good])**2 + (data[:,color2+nfilts][good])**2)

    print("\n  Plotting "+ylabel+" vs. "+xlabel+" CMD for "+str(len(ymag))+" / "+str(nstars)+" stars from "+photfile+".\n  File is cmd_iso.png.\n")

    N = 0.6
    #print(np.max(ymag),np.min(ymag))
    if ((np.min(ymag) < 0) & (abs(np.min(ymag)) < 0.1)):
        N = 20
    if ((np.min(ymag) < 0) & (abs(np.min(ymag)) > 0.1)):
        N = 2.0
    if ((np.min(ymag) > 0) & (abs(np.min(ymag)) < 0.1)):
        N = 0.01
    #print(1.1*np.max(ymag),N*np.min(ymag))

    mag_range = abs(N*np.min(ymag)-1.1*np.max(ymag))
    lab_start_y = mag_range/2.1
    dy = mag_range/40

    plt.figure(figsize=(6,10))
    plt.errorbar(color,ymag,xerr=col_err,yerr=mag_err,fmt='o',markersize=1,color='k')
    plt.ylim(1.1*np.max(ymag),N*np.min(ymag))
    plt.xlim([np.min(color)-0.5,np.max(color)+2.5])
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(clusterName)
    plt.plot(color_mod,y_model,'r--')
    plt.text(np.max(color),lab_start_y-4*dy,'log(Age) = '+age_lab)
    plt.text(np.max(color),lab_start_y-3*dy,'Age (Gyr) = '+age_lab_l)
    plt.text(np.max(color),lab_start_y-2*dy,'[Fe/H] = '+feh_lab)
    plt.text(np.max(color),lab_start_y-dy,'(m - M)_V = '+dis_lab)
    plt.text(np.max(color),lab_start_y,'A_V = '+av_lab)

    plt.savefig(root_dir + 'cmd_iso.png')
    plt.show()
    plt.cla()
    plt.clf()
    plt.close()

    # display results to screen
    print("log(Age):  "+age_lab)
    print("Age (Gry): "+age_lab_l)
    print("[Fe/H]:    "+feh_lab)
    print("dist:      "+dis_lab)
    print("Av:        "+av_lab)

    # end




def make_thinCMD(model_file):
    output_file = model_file+'.thin'

    # read in the data
    fo = open(model_file,"r")
    header=fo.readline()
    fo.close()
    #header = np.loadtxt(model_file,dtype=str,max_rows=1,delimiter=',')
    data   = np.loadtxt(model_file,dtype=str,skiprows=1)

    # read in line-by-line
    with open(model_file,"r") as ins:
        array = []
        for line in ins:
            array.append(line)

    mass = np.around(data[1:,0].astype(float),4)

    # set up initital values
    mass_prev = 0.0
    i         = []

    for idx,val in enumerate(mass):
        mass_diff = val - mass_prev
        if mass_diff == 0:
            continue
        #print(idx,val)
        i.append(idx)   # save the values that are non-repeats
        mass_prev = val

    outfile = open(output_file,'w')
    outfile.write(array[0])

    # write the file
    for j in i:
        outfile.write(array[j+1])

    outfile.close()


def createAllCMDS(mag,color,root_dir):
    files = os.listdir(root_dir)
    for f in files:
        if f.removesuffix(".ms") != f and f.removesuffix(".thin.ms") == f:
            #file is a simple ms file 
            version = f.removesuffix(".ms")
            version = version.removeprefix(root_dir)
            print(f"CMD Diagram for version: {version}")
            create_iso_plot(mag=mag,color=color, clusterName= root_dir, version = version)
        


