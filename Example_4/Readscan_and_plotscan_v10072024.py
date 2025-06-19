import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import kaiser
from matplotlib.widgets import Slider, Button, TextBox  #, CheckButtons
from tkinter import Tk # in Python 2, use "Tkinter" instead
import scipy.constants as const
import OpenBinary_module_v11072024 as obm
# Created by Sven Herbers 10.07.2024 - Sven_herbers@web.de
##### INPUT

label    = "a_"  # Label before index # if JUG data, idx is not used and single file is assumed.
JUGDATA=False # If False, Ftmw++ .dat file is loaded, if True, JUG file is loaded.
PowerSpectrum=False #Shows the square of the usual amplitude spectrum - like the HAUPTPROGRAMM.exe does for JUGDATA. If False the usual Amplitude spectrum will be shown.

masklines=[]       # lines to be removed
maskwidth=0.1       # +-width will be removed by mask in MHz.
startidx =  2     # Starting index, you can ommit or add the leading zeros, it works either way. the file indices must consist of 5 digits.
endidx   = 5      # End index, again, the leading zeros are not nessecary but can be added for fun
forcecalc= 1       # if zero, the script checks if it already ran once and produced a spectrum file and plots that instead of recalculating. If 1, recalculation is forced.
rangeMHz = 501    # plus and minus this value around probe frequency is what will be added from each scan step. About HWHM=FWHM/2=0.5MHz should be good. If your mode ran away a bit, you might want to increase this values
                   # For chirps: about chirp width/2 should be good. e.g. for 1GHz chirp rangeMHz= 500.
##### ADVANCED INPUT

Kaiser = 0.0      # application of a kaiserwindow with kaiserparameter Kaiser=beta. If zero, no kaiserwindow is used. If you change this, you might want to set forcecalc=1 to overwrite previous run.
Kaiser_Tail= 0#     # tailcut applied in mus
Kaiser_Head = 0#   # Headcut applied in mus
rebin = False          # In case you want to reduce resolution, you can now with the rebin option. rebin=0 omits rebinning. The produced line lists from stickspectrum, are then also inaccurate.
rebinsize = 0.250
stickspectrum = True #prdouce a stickspectrum with stickthreshold
stick_dopplerread = True # if this is turned on an approach is made to find doppler pairs based on stick_exampledoppler
##### END OF INPUT


######## Advanced control parameters
rezero = 1.0           # substract noisefloor by the inverse squared value weighted average of each section.
stickthreshold = 0.5   # threshold in terms of global spectrum. local threshold hardwired to 2.5 for three datapoints at and around maximum.
stickmaxtest = 5       # checks for every point, stickmaxtest left and right of it, if it is the largest, it is identified as a local maximum and subjected to interpolation
stick_exampledoppler = [[ 13341.22734,             0.03578],[  9866.81114,             0.02523]] #  - uses two dopplersplittings (+- from centre frequency)  with two frequencies, to calculate the expected dopplersplit throughout the spectrum
stick_dopplerthresh = 0.011# defines how precise in MHz the doppler must be to identify a pair - depends on how much fine structure you have but 15 kHz should be fine in most cases
# exp window not tested:
exponentialwindow=ew = 0 # exponentialwindow parameter , 0 deactivates it
def window(fid, tailcut, headcut, timestep, beta,exponentialwindow):
    'The Kaiserwindow uses Besselfunctions - luckily it is already availabe from scipy. FID is cut in half assuming 2 channels'
    tailing=int(tailcut*1000/timestep)#asumming 0.8ns timestep, the tail is getting cut
    heading=int(headcut*1000/timestep)
    L = len(fid) // 2
    time=np.linspace(heading,L-1-tailing-heading,L-tailing)


    ew=exponentialwindow
    if ew==0:
        expw=1
    else:
        tau=len(time)/(ew)
        expw=np.exp(-timestep/1000*abs(time-(len(time)-1))*(1/tau))
    if tailing==0 and heading==0:
        k=kaiser(L,beta)
        side1=np.array(fid[0:L])*k*expw
        side2=np.array(fid[L:])*k*expw
        fid[0:L]=side1
        fid[L:]=side2
    if tailing!=0 or heading !=0:
        k=kaiser(L-tailing-heading,beta)
        side1=np.array(fid[0+heading:L-tailing])*k*expw
        side2=np.array(fid[L+heading:2*L-tailing])*k*expw
        fidnew=np.zeros((L-tailing-heading)*2)
        fidnew[0:L-tailing-heading]=side1
        fidnew[L-tailing-heading:]=side2
        fid=fidnew
    return fid
def addtobins(frequencies, freq_bins, bin_counts, spec_freq, spec_int, probe, rangeMHz):
    'takes frequency array and intensity array of spectrum and puts it into the bins'
    startidx = 0
    stopidx = -1
    for j in range(len(spec_freq)):
        if startidx == 0 and spec_freq[j] > (probe - rangeMHz):
            startidx = j
        elif stopidx == -1 and spec_freq[j] > (probe + rangeMHz):
            stopidx = j
            break
    spec_freq_new = spec_freq[startidx:stopidx]
    #print(spec_freq_new)
    spec_int_new = spec_int[startidx:stopidx]

    minfreq = spec_freq_new[0]
    maxfreq = spec_freq_new[-1]
    #centrefreq = probe
    idxstart = 0
    idxend = 0
    for j in range(len(frequencies)):
        if idxstart == 0 and frequencies[j] > minfreq:
            idxstart = j
        elif idxend == 0 and frequencies[j] > maxfreq:
            idxend = j
            break
    factor = (1.0 * (idxend - idxstart)) / (1.0 * len(spec_freq_new))
    for j in range(len(spec_freq_new)):
        idxguess = idxstart + factor * j - 10
        for k in range(int(idxguess), len(frequencies)):
            if frequencies[k] > spec_freq_new[j]:
                fup = frequencies[k]
                fdown = frequencies[k - 1]
                dif = fup - fdown
                dif2 = spec_freq_new[j] - fdown
                #print(dif)
                #print(dif2)
                factor_up = dif2 / dif
                factor_down = 1 - factor_up
                if factor_up > 1 or factor_down > 1:
                    print("Problem in binning! can not assign bins properly!")
                idx_bin = k
                break
        freq_bins[k] += factor_up * spec_int_new[j]
        freq_bins[k - 1] += factor_down * spec_int_new[j]
        bin_counts[k] += factor_up
        bin_counts[k-1] += factor_down

    return freq_bins, bin_counts

def evaluate_scan(Leadstring,startindex,endindex,overwrite):
    "simple function to run through all scanfiles and returns frequencies, and intensities.If overwrite=0 it will check if a spectrum is already available"
    if overwrite==0:
        try:

            a=np.loadtxt(Leadstring+"{:05d}_to_".format(startindex)+"{:05d}".format(endindex)+"compiledSpectrum.txt")
            frequencies, freq_bins=a[:,0],a[:,1]
            overwrite=-1
            print("FILE FOUND")
            print(Leadstring + "{:05d}_to_".format(startindex) + "{:05d}".format(endindex) + "compiledSpectrum.txt")
            print("CONTINUE PLOTTING, WITHOUT REPROCESSING SCAN FILES")
        except:
            pass
    if overwrite!=-1:
        informations0,numbers0=obm.convertfile(Leadstring+"{:05d}.dat".format(startindex),outputnameview="view.txt",outputnamefid="fid0.txt",writefid=0,writeview=0)
        print(informations0[3])
        freq0,int0=obm.createspectrum(fidtype=0,fidfile=[informations0,numbers0],paddingfactor=3,outname="Freq_spectrum",writespectrum=0,writefid=0)
        probe0=informations0[6]

        informationsend,numbersend=obm.convertfile(Leadstring+"{:05d}.dat".format(endindex),outputnameview="view.txt",outputnamefid="fid1.txt",writefid=0,writeview=0)
        #freqend,intend=obm.createspectrum(fidfile="fid.txt",paddingfactor=3,outname="Freq_spectrum",writefid=1)
        probeend=informationsend[6]
        #print(probe0)
        #print(probeend)
        therange = range(startindex,endindex+1)
        rangetobin = rangeMHz#5 # 0.5 means +- means using 1MHz of each measurment - the size of the program window
        binsize = (freq0[1] - freq0[0]) * 2
        if rebin==1:
            binsize=rebinsize

        if probeend<probe0:
            frequencies = np.arange(probeend - 2 * rangetobin, probe0 + 2 * rangetobin, binsize)
        else:
            frequencies = np.arange(probe0 - 2 * rangetobin, probeend + 2 * rangetobin, binsize)

        freq_bins = np.copy(frequencies) * 0
        bin_counts = np.copy(frequencies) * 0
        #rangefreq=abs(probeend-probe0)
        oldeprobe=probe0
        for k in therange:
            informations,numbers = obm.convertfile(Leadstring+"{:05d}.dat".format(k), outputnameview="view.txt",
                                              outputnamefid="fid{:d}.txt".format(k%2),writefid=0,writeview=0) # i alternate between 0 and 1 here, because sometimes, computer is a bit slow with read/write. For speedup, avoid save and directly pass to next function can be implemented...
            if Kaiser != 0 or Kaiser_Tail!=0 or exponentialwindow !=0 or Kaiser_Head!=0:
                timestep=informations[3]*1000

                numbers=window(numbers,Kaiser_Tail,Kaiser_Head,timestep,Kaiser,exponentialwindow)
            freqk,intk=obm.createspectrum(fidtype=0,fidfile=[informations,numbers],paddingfactor=3,outname="Freq_spectrum",writespectrum=0,writefid=0)
            probek = informations[6]
            stepsize=probek-oldeprobe
            oldeprobe=probek
            print("Step {:12.5f} MHz - index {:05d}".format(stepsize,k))
            if rezero==1:
                newbins=np.copy(freq_bins)*0
                newbincounts=np.copy(bin_counts)*0
                newbins, newbincounts = addtobins(frequencies, newbins, newbincounts, freqk, intk, probek, rangetobin)
                idxrange=np.nonzero(newbins)
                idxrange=np.sort(idxrange)[0]
                print(np.shape(idxrange))
                weights=np.copy(newbins)*0
                weights[idxrange[1:-1]]=newbins[idxrange[1:-1]]**-2
                av=np.sum(weights*newbins)/np.sum(weights)
                print(av)
                newbins[idxrange]=newbins[idxrange]-av
                newbins[idxrange[0]]=0
                newbins[idxrange[-1]]=0
                freq_bins=freq_bins+newbins
                bin_counts=bin_counts+newbincounts
            else:
                freq_bins,bin_counts=addtobins(frequencies,freq_bins,bin_counts,freqk,intk,probek,rangetobin)
        for j in range(len(bin_counts)):
            if bin_counts[j]==0:
              bin_counts[j]=-1
        freq_bins=freq_bins/bin_counts
        f=open(Leadstring+"{:05d}_to_".format(startindex)+"{:05d}".format(endindex)+"compiledSpectrum.txt",'w')
        if rezero==1:
            freq_bins=freq_bins-np.min(freq_bins)
        for j in range(len(frequencies)):
            f.write(" {:23.15e} {:23.15e}\n".format(frequencies[j],freq_bins[j]))
        f.close()
    if overwrite != -1:
     for entry in masklines:
        entrylo=entry-maskwidth
        entryhi=entry+maskwidth
        idxlo=int((entrylo-frequencies[0])/binsize)
        idxhi=int((entryhi-frequencies[0])/binsize)+1
        frequenciesloav=np.average(freq_bins[idxlo-4:idxlo-1])
        frequencieshiav=np.average(freq_bins[idxhi+1:idxhi+4])
        freq_bins[idxlo:idxhi+1]=0.5*(frequenciesloav+frequencieshiav)
     f=open(Leadstring+"{:05d}_to_".format(startindex)+"{:05d}".format(endindex)+"compiledSpectrum_MASKED.txt",'w')
     for j in range(len(frequencies)):
            f.write(" {:23.15e} {:23.15e}\n".format(frequencies[j],freq_bins[j]))
     f.close()
    return frequencies, freq_bins



#### Here is the runcommand.
#### Start with inital string, then starting index, then end index.
#### Last value is zero if you would like to check for an already processed file, 1 if you want to force a recalculation
maskwidth=abs(maskwidth)
if endidx==startidx or JUGDATA==True:
    Leadstring=label
    if JUGDATA==True:
        print("Jugdata currently only implemented for single file. Trying to open file with filename "+Leadstring )
    else:
        print("Endindex is Startindex, reading only single file")

    k=startidx

    if JUGDATA==True:
            informations, numbers = obm.convertfileJUGDATA(Leadstring, outputnameview="view.txt",
                                        outputnamefid="fid{:d}.txt".format(k % 2), writefid=0,
                                        writeview=0)
    else:
            informations, numbers = obm.convertfile(Leadstring + "{:05d}.dat".format(k), outputnameview="view.txt",
                                        outputnamefid="fid{:d}.txt".format(k % 2), writefid=0,
                                        writeview=0)  # i alternate between 0 and 1 here, because sometimes, computer is a bit slow with read/write. For speedup, avoid save and directly pass to next function can be implemented...
    #print(informations[3])
    if Kaiser != 0 or Kaiser_Tail != 0 or exponentialwindow != 0 or Kaiser_Head!=0:
     if JUGDATA==True:
         timestep=informations[0]
     else:
        timestep=informations[3]*1000
     numbers = window(numbers, Kaiser_Tail, Kaiser_Head, timestep, Kaiser, exponentialwindow)
    if JUGDATA==True:
            freqk, intk = obm.createspectrumJUGDATA(fidfile=[informations, numbers], paddingfactor=3, outname="Freq_spectrum",
                               writespectrum=0, writefid=0)
            probe = informations[1]
    else:
        freqk, intk = obm.createspectrum(fidtype=0, fidfile=[informations, numbers], paddingfactor=3, outname="Freq_spectrum",
                                   writespectrum=0, writefid=0)
        probe = informations[6]
    idx1=0
    idx2=-1
    for j in range(0,len(freqk)):
        if freqk[j]>probe-rangeMHz:
            idx1=j
            break
    for j in range(0,len(freqk)):
        if freqk[j]>probe+rangeMHz:
            idx2=j-1
            break
    frequencies,freq_bins=freqk[idx1:idx2],intk[idx1:idx2]

    if rebin == 1:
        binsize = rebinsize
        rangetobin=rangeMHz

        frequencies = np.arange(probe - 2 * rangetobin, probe+ 2 * rangetobin, binsize)
        freq_bins = np.copy(frequencies) * 0
        bin_counts = np.copy(frequencies) * 0
        freq_bins, bin_counts = addtobins(frequencies, freq_bins, bin_counts, freqk, intk, probe, rangetobin)


    f=open(Leadstring+"{:05d}_to_".format(startidx)+"{:05d}".format(endidx)+"compiledSpectrum.txt",'w')
    if rezero==1:
        freq_bins=freq_bins-np.min(freq_bins)
    for j in range(len(frequencies)):
         f.write(" {:23.15e} {:23.15e}\n".format(frequencies[j],freq_bins[j]))
    f.close()
    if rebin != 1:
        binsize=frequencies[1]-frequencies[0]
    for entry in masklines:
        entrylo=entry-maskwidth
        entryhi=entry+maskwidth
        idxlo=int((entrylo-frequencies[0])/binsize)
        idxhi=int((entryhi-frequencies[0])/binsize)+1
        frequenciesloav=np.average(freq_bins[idxlo-4:idxlo-1])
        frequencieshiav=np.average(freq_bins[idxhi+1:idxhi+4])
        freq_bins[idxlo:idxhi+1]=0.5*(frequenciesloav+frequencieshiav)
    f=open(Leadstring+"{:05d}_to_".format(startidx)+"{:05d}".format(endidx)+"compiledSpectrum_MASKED.txt",'w')
    for j in range(len(frequencies)):
            f.write(" {:23.15e} {:23.15e}\n".format(frequencies[j],freq_bins[j]))
    f.close()
    x, y = frequencies, freq_bins
    if PowerSpectrum==True:
        y=y**2
else:
    x,y=evaluate_scan(label,startidx,endidx,forcecalc)



##plt.plot(x,y)
##plt.show()
print("Survey merging finished. Continue with plotting for peak finding.")
#### From here on just plotting with buttons and peakfinder etc.
def copypaste(number):
    def copier(event):
        if number==0:
            print("Error in copypaste")
        elif number==1:
            value = (T1.label.get_text())
        elif number==2:
            value = (T2.label.get_text())
        elif number==3:
            value = (T3.label.get_text())
        elif number==4:
            value = (T4.label.get_text())
        elif number==5:
            value = (T5.label.get_text())
        else:
            print("Error in copypaste 2")
        r = Tk()
        r.withdraw()
        r.clipboard_clear()
        r.clipboard_append(value)
        r.update()  # now it stays on the clipboard after the window is closed
        r.destroy()
    return copier

def line_at_cursor(event):
    if event.dblclick == True:
        cur_ylim=ax.get_ylim()
        cur_xlim=ax.get_xlim()

        xdata = event.xdata  # get event x location
        ydata = event.ydata  # get event y location
        #print(ydata)
        #print(event.button)
        mina = np.min(x)
        maxa = np.max(x)
        counts = len(x)
        step = (maxa - mina) / (counts)
        idxguess = (xdata - mina) / step
        if idxguess > 10:
            idxguess = idxguess - 10
        idxguess = int(idxguess)
        #print(idxguess)
        #print(a[:,0][idxguess])
        #print(xdata)
        for j in range(10000):
            if xdata < x[idxguess + j]:
                centreidx = idxguess + j
                break
        #print(x[centreidx])
        #x[centreidx - 2:centreidx + 3]
        #y[centreidx - 2:centreidx + 3]
        counts=0
        while counts<20:
            counts+=1
            checkarray=y[centreidx - 2:centreidx + 3]
            #print(len(checkarray))
            checkmax=np.argmax(checkarray)
            if checkmax==0:
                centreidx-=1
            elif checkmax==4:
                centreidx+=1
            else:
                break
        idxmax =centreidx - 2+ np.argmax(y[centreidx - 2:centreidx + 3])
        threepointinterpol0=tpi0=x[idxmax-1:idxmax+2]
        threepointinterpol0=tpi1=y[idxmax-1:idxmax+2]
        resort=np.argsort(tpi1)
        tpi1=tpi1[resort]
        tpi0=tpi0[resort]
        #print(tpi1)
        #print(tpi1[0]/tpi1[2])
        if tpi1[0]/tpi1[2] < 0.5:
            #print("Using linear interpolation")
            #linear interpolation if datapoint are sparse, works better than quadratic at far distances.
            weightedx=((tpi1[1]-tpi1[0])*tpi0[1]+(tpi1[2]-tpi1[0])*tpi0[2])/((tpi1[1]-tpi1[0])+(tpi1[2]-tpi1[0]))
            weightedy=((tpi1[1]-tpi1[0])**2+(tpi1[2]-tpi1[0])**2)/((tpi1[1]-tpi1[0])+(tpi1[2]-tpi1[0]))+tpi1[0]
        else:
            #print("Using quadratic interpolation")
            # quadratic interpolation if datapoints are close, the quadratic interpolation matches very well e.g. a lorentzian shape around its peak
            weightedx=((tpi1[1]-tpi1[0])**2*tpi0[1]+(tpi1[2]-tpi1[0])**2*tpi0[2])/((tpi1[1]-tpi1[0])**2+(tpi1[2]-tpi1[0])**2)
            weightedy=((tpi1[1]-tpi1[0])**3+(tpi1[2]-tpi1[0])**3)/((tpi1[1]-tpi1[0])**2+(tpi1[2]-tpi1[0])**2)+tpi1[0]


        #print(len(tpi0))

        valuex = weightedx
        valuey = weightedy

        if event.button==1:
            marker1.set_ydata([valuey])
            marker1.set_xdata([valuex])
            marker1b.set_ydata([0,valuey])
            marker1b.set_xdata([valuex,valuex])
            T1.label.set_text("{:11.5f}".format(valuex))
            CheckMarkers[0]=1
        elif event.button==3:
            marker2.set_ydata([valuey])
            marker2.set_xdata([valuex])
            marker2b.set_ydata([0,valuey])
            marker2b.set_xdata([valuex,valuex])
            T2.label.set_text("{:11.5f}".format(valuex))
            CheckMarkers[1]=1
        if CheckMarkers[0]==1 and CheckMarkers[1]==1:
            valuea=marker1.get_xdata()[0]
            valueb=marker2.get_xdata()[0]
            inti1=marker1.get_ydata()[0]
            inti2=marker2.get_ydata()[0]
            intaverage=(inti1+inti2)/2
            average=(valuea+valueb)/2
            doppler=abs((valuea-valueb)/2)
            T3.label.set_text("{:11.5f}".format(doppler))
            T4.label.set_text("{:11.5f}".format(average))
            T5.label.set_text("{:8.3e}".format(intaverage))
        plt.draw() # force re-draw
fig, ax = plt.subplots(figsize=(18,9))
line0, = ax.plot(x,y, color='k',lw=1)
if stickspectrum==1:
    av=np.average(y)
    thresh=av*stickthreshold
    s=np.copy(x)*0
    vx=[]
    vy=[]
    #print("hello")
    #print(y[10-stickmaxtest:10+stickmaxtest+1])
    for j in range(stickmaxtest,len(x)-stickmaxtest-1):
        if y[j]>thresh:
          if np.min(y[j-stickmaxtest:j+stickmaxtest+1]!=0):
            if np.max(y[j-stickmaxtest:j+stickmaxtest+1])==y[j]:
                if j>500:
                    idx1=j-500
                else:
                    idx1=0
                if j <len(x)-500:
                    idx2=j+500
                else:
                    idx2=-1
                localaverage=np.average(y[idx1:idx2])
                t=2.5*localaverage
                if y[j]>t and y[j+1]>t and y[j-1]>t: # peakintensity must be at least t
                    tpi0 = np.array([x[j - 1], x[j], x[j + 1]])
                    tpi1 = np.array([y[j - 1], y[j], y[j + 1]])
                    resort = np.argsort(tpi1)
                    tpi1 = tpi1[resort]
                    tpi0 = tpi0[resort]
                    # print(tpi1[resort])
                    if tpi1[0] / tpi1[2] < 0.5:
                        # print("Using linear interpolation")
                        # linear interpolation if datapoint are sparse, works better than quadratic at far distances.
                        weightedx = ((tpi1[1] - tpi1[0]) * tpi0[1] + (tpi1[2] - tpi1[0]) * tpi0[2]) / (
                                    (tpi1[1] - tpi1[0]) + (tpi1[2] - tpi1[0]))
                        weightedy = ((tpi1[1] - tpi1[0]) ** 2 + (tpi1[2] - tpi1[0]) ** 2) / (
                                    (tpi1[1] - tpi1[0]) + (tpi1[2] - tpi1[0])) + tpi1[0]
                    else:
                        # print("Using quadratic interpolation")
                        # quadratic interpolation if datapoints are close, the quadratic interpolation matches very well e.g. a lorentzian shape around its peak
                        weightedx = ((tpi1[1] - tpi1[0]) ** 2 * tpi0[1] + (tpi1[2] - tpi1[0]) ** 2 * tpi0[2]) / (
                                    (tpi1[1] - tpi1[0]) ** 2 + (tpi1[2] - tpi1[0]) ** 2)
                        weightedy = ((tpi1[1] - tpi1[0]) ** 3 + (tpi1[2] - tpi1[0]) ** 3) / (
                                    (tpi1[1] - tpi1[0]) ** 2 + (tpi1[2] - tpi1[0]) ** 2) + tpi1[0]
                    # print(len(tpi0))

                    vx+= [weightedx]
                    vy+= [weightedy]
                    #s[j]=y[j]
                else:
                    pass
    vx=np.array(vx)
    vy=np.array(vy)
    if stick_dopplerread == False:
        ax.errorbar(vx,vy,yerr=[vy,vy*0],ls="",color='red',capsize=0)
    else:
        #ax.errorbar(vx, vy, yerr=[vy, vy * 0], ls="", color='red', capsize=0)
        pass
    g=open(label+"{:05d}_to_".format(startidx)+"{:05d}".format(endidx)+"STICKS_single.txt",'w')
    for j in range(len(vx)):
        g.write(" {:23.15e} {:23.15e}\n".format(vx[j],vy[j]))
    g.close()
    if stick_dopplerread==True:
        dopplersplitting = stick_exampledoppler[0][1] / stick_exampledoppler[0][0] + stick_exampledoppler[1][1] / \
                           stick_exampledoppler[1][0]
        print('The velocity of the jet is assumed as {:10.0f}'.format(dopplersplitting/2*const.c))
        Doppler_linesx=[]
        Doppler_splits=[]
        Doppler_linesy=[]
        for j in range(len(vx)):
            for k in range(j,len(vx)):
                #print("try")
                if abs(vx[k]-vx[j]) > dopplersplitting*vx[j]+2*stick_dopplerthresh:
                    break
                if abs(abs(vx[k]-vx[j])-dopplersplitting*vx[j])<stick_dopplerthresh:
                    #print("hit")
                    Doppler_linesx += [(vx[k]+vx[j])/2]
                    Doppler_splits += [(vx[k] - vx[j]) / 2]
                    Doppler_linesy += [(vy[k] + vy[j]) / 2]
        Doppler_linesx=np.array(Doppler_linesx)
        Doppler_splits=np.array(Doppler_splits)
        Doppler_linesy=np.array(Doppler_linesy)
        ax.errorbar(Doppler_linesx-Doppler_splits, Doppler_linesy, yerr=[Doppler_linesy, Doppler_linesy * 0], ls="", color='green', capsize=0)
        ax.errorbar(Doppler_linesx + Doppler_splits, Doppler_linesy, yerr=[Doppler_linesy, Doppler_linesy * 0], ls="", color='green', capsize=0)
        g=open(label+"{:05d}_to_".format(startidx)+"{:05d}".format(endidx)+"STICKS_DOPPLER.txt",'w')
        for j in range(len(Doppler_linesx)):
            g.write(" {:23.5f} {:13.5e} {:13.5f}\n".format(Doppler_linesx[j],Doppler_linesy[j],Doppler_splits[j]))
        g.close()

marker1, = ax.plot([x[0]],[y[0]], color='red',lw=1, marker='o')
marker1b, = ax.plot([x[0]],[y[0]], color='red',lw=1)
marker2, = ax.plot([x[0]],[y[0]], color='blue',lw=1, marker='o')
marker2b, = ax.plot([x[0]],[y[0]], color='blue',lw=1)
#plt.plot(a[:,0],a[:,1])
TextBoxAxes1 = fig.add_axes([0.83, 0.84, 0.07, 0.04])
T1=Button(TextBoxAxes1,"LeftDoppler")
TextBoxAxes2 = fig.add_axes([0.83, 0.79, 0.07, 0.04])
T2=Button(TextBoxAxes2,"RightDoppler")
TextBoxAxes3 = fig.add_axes([0.83, 0.74, 0.07, 0.04])
T3=Button(TextBoxAxes3,"DopplerSplit")
TextBoxAxes4 = fig.add_axes([0.83, 0.69, 0.07, 0.04])
T4=Button(TextBoxAxes4,"CenterFreq")
TextBoxAxes5 = fig.add_axes([0.74, 0.84, 0.07, 0.04])
T5=Button(TextBoxAxes5,"Intensity")
fig.canvas.mpl_connect('button_press_event', line_at_cursor)
CheckMarkers=[0,0]
intis=[0,0]
T1.on_clicked(copypaste(1))
T2.on_clicked(copypaste(2))
T3.on_clicked(copypaste(3))
T4.on_clicked(copypaste(4))
T5.on_clicked(copypaste(5))
ax.set_xlabel("Frequency / MHz")
ax.set_ylabel("Intensity / mVrms")
plt.show()
#print("{:05d}".format(1))
