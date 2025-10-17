import numpy as np
import struct
from scipy.fft import fft, fftfreq
# Specify the size of each chunk to read
#### Lines in "view" file - format dddddddd:
#line 22: 0,0, span/µs, Interval/µs, Bandwidth/ MHz, Resolution/kHz, Probe /MHz,0
#line 75: Protection, recovery delay, recovery width, transfer delay, transfer width, excitation delay, excitation width
#line 48: 0,0,0,0, Tune: Span/Mode Separation, Tune/Span *s^-1, level, Crossover
def convertfile(filename,outputnameview,outputnamefid,writefid,writeview):
    'Function to convert Jens binary file to readable text. Saves complete file as well as extracted FID. Returns some metadata.'
    fn=filename
    ov=outputnameview
    of=outputnamefid
    chunk_size = 64
    start=0
    file = open(fn, "rb")
    if writeview==1:
        viewport=open(ov,'w')
    HeaderList=[]
    Numberlist=[]
#labelblock=["string","integer"
# Using while loop to iterate the file data
    while True:
        offset1=16
        offset2=0
        if start==192:
            chunk = file.read(chunk_size-offset1)
        elif start==393:
            chunk_size = 64
            chunk = file.read(chunk_size - offset2)
        else:
            chunk = file.read(chunk_size)

        if not chunk:
            break
        # Processing the chunk of binary data
        if start > -1:
            try:
              a=f"{chunk.decode('utf8')}"
            except:
                    a=f"{chunk}"
            if writeview==1:
                viewport.write(a+"\n")
            if start == 20:
                samples = struct.unpack('iiiiiiiiiiiiiiii', chunk)[3]
            if start == 21:
                informations=list(struct.unpack('dddddddd',chunk))
                informations+=[samples]
                HeaderList+=["Span:       {:10.3f} µs\n".format(informations[2])]
                HeaderList+=["Interval:   {:10.3f} ns\n".format(informations[3]*1000)]
                HeaderList+=["Bandwidth:  {:10.3f} MHz\n".format(informations[4])]
                HeaderList+=["Resolution: {:10.4f} kHz\n".format(informations[5]*1000)]
                HeaderList+=["Probe:      {:10.3f} MHz\n".format(informations[6])]
                HeaderList+=["Conversion: {:10.3f} MHz\n".format(informations[7])]
                HeaderList+=["Samples: {:13.0f}    \n".format(informations[8])]
            if start > 393: #and start < 399:
               try:
                values=(struct.unpack('dddddddd',chunk))
               except:
                print("chunk not matching -adjusting ")
                try:
                 values=(struct.unpack('ddddddd',chunk))
                 print("56")
                except:
                 try:
                  values=(struct.unpack('dddddd',chunk))
                  print("48")
                 except:
                  try:
                   values=(struct.unpack('ddddd',chunk))
                   print("40")
                  except:
                   try:
                     values=(struct.unpack('dddd',chunk))
                     print("32")
                   except:
                    try:
                     values=(struct.unpack('ddd',chunk))
                     print("24")
                    except:
                     try:
                      values=(struct.unpack('dd',chunk))
                      print("16")
                     except:
                      try:
                       values=(struct.unpack('d',chunk))
                       print("8")
                      except:
                         print("Error in Loading FID, chunksize problem. Passing for now.")

               for m in range(len(values)):
                   Numberlist += [values[m]]
            start+=1
    if writeview==1:
        viewport.close()
    if writefid==1:
        fid=open(of,'w')
        for element in HeaderList:
            fid.write(element)
        for n in Numberlist:
            fid.write(" {:23.15e}\n".format(n))
        fid.close()
    return informations, Numberlist
def convertfileJUGDATA(filename,outputnameview,outputnamefid,writefid,writeview):
    'Function to convert Jens binary file to readable text. Saves complete file as well as extracted FID. Returns some metadata.'
    fn=filename
    ov=outputnameview
    of=outputnamefid
    chunk_size = 32
    start=0
    file = open(fn, "rb")
    if writeview==1:
        viewport=open(ov,'w')
    Numberlist=[]
    informations=[]
# Using while loop to iterate the file data
    while True:
        chunk = file.read(chunk_size)
        if not chunk:
            break
        # Processing the chunk of binary data
        if start > -1:
            try:
              a=f"{chunk.decode('utf8')}"
            except:
              a=f"{chunk}"
            if writeview==1:
                viewport.write(a+"\n")
            if start==8:
                step=struct.unpack('llllllll',chunk)[0]
                informations+=[step]
            if start==10:
                centrefreq=struct.unpack('dddd',chunk)[0]
                informations+=[centrefreq]
            if start > 15: #and start < 399:
               try:
                     values=(struct.unpack('llllllll',chunk))
               except:
                    try:
                     values=(struct.unpack('llllll',chunk))
                     print("24")
                    except:
                     try:
                      values=(struct.unpack('llll',chunk))
                      print("16")
                     except:
                      try:
                       values=(struct.unpack('ll',chunk))
                       print("8")
                      except:
                         print("Error in Loading FID, chunksize problem. Pasing for now.")

               for m in range(len(values)):
                        Numberlist+=[values[m]]


            start+=1
    if writeview==1:
        viewport.close()
    if writefid==1:
        fid=open(of,'w')
        for n in range(len(Numberlist)):
            fid.write(" {:23.15e}\n".format(n))
        fid.close()
    return informations, Numberlist

def createspectrum(fidtype=1, fidfile="fid.txt",paddingfactor=3,outname="Freq_spectrum",writespectrum=0,writefid=0):
    'Function to read the fid created by convertfile, allows for zero padding for extra datapoints returns spectrum. fidtype=0 - list of [informations, Numberlist], fidtype=1 - string representing fid.txt file name'
    pf=paddingfactor  # 0 = no padding; else extend time-domain array by (pf+1) factor
    ic=intensityconversion=1/np.sqrt(2)# if 1.0 - peak voltage, if 1/sqrt2 - rms voltage
    step = 1.0
    FrequencyOffset = 0
    if fidtype==1:

        headlines=7
        f=open(fidfile,'r')

        for k in range(headlines):
            line=f.readline()
            splitty=line.split()
            if "Interval:" in line:
                step=np.float64(splitty[-2])
            if "Probe:" in line:
                FrequencyOffset=np.float64(splitty[-2])
            if "Conversion:" in line:
                FrequencyOffset+=-np.float64(splitty[-2])
            if "Samples:" in line:
                samples = np.float64(splitty[-1])
        a=np.loadtxt(fidfile,skiprows=headlines,dtype=np.complex128)
    elif fidtype==0:
        [ifo,nums]=fidfile
        FrequencyOffset=Probe=ifo[6]
        Conversion=ifo[7]
        FrequencyOffset+=-Conversion
        step=ifo[3]*1000
        samples=ifo[8]
        a=np.array(nums,dtype=np.complex128)
    else:
        print("Error in Choice of fidtype - use 0 for [informations, Numberlist] as input. Use 1 for string as inputfilename.")
    a=a/samples*1000
    counts=len(a)//2
    side1=a[0:counts]
    side2=a[counts:2*counts]

    combined=side1+1.0j*side2
    combined=np.pad(combined, (0, len(combined)*pf), 'constant')

    x=np.linspace(0,len(combined)-1,len(combined))
    x=x*step/1000

    xf = fftfreq(len(combined), step*10**-9)
    xf=xf/10**6+FrequencyOffset

    yfc=fft(combined)/counts

    yfc=np.abs(yfc)
    yfc=yfc*ic

    idx=np.argsort(xf)
    xf=xf[idx]
    yfc=yfc[idx]

    if writespectrum==1:
        output=open(outname+".txt",'w')
        for j in range(0,len(xf)):
            output.write("{:23.15e} {:23.15e}\n".format(xf[j],yfc[j]))
        output.close()

    if writefid==1:
     output2=open(outname+"_fid.txt",'w')
     xshort=x[:len(x)//(pf+1)]
     for j in range(0,len(xshort)):
        output2.write("{:23.15e} {:23.15e} {:23.15e}\n".format(xshort[j],np.real(side1[j]),np.real(side2[j])))
     output2.close()

    return xf,yfc # returns frequency and intensity

def createspectrumJUGDATA(fidfile="fid.txt",paddingfactor=3,outname="Freq_spectrum",writespectrum=0,writefid=0):
    'Function to read the fid created by convertfile, allows for zero padding for extra datapoints returns spectrum. fidtype=0 - list of [informations, Numberlist], fidtype=1 - string representing fid.txt file name'
    pf=paddingfactor  # 0 = no padding. else increase site of time domain signal array by factor of (pf+1)
    ic=intensityconversion=1/np.sqrt(2)# if 1.0 - peak voltage, if 1/sqrt2 - rms voltage
    [ifo,nums]=fidfile
    FrequencyOffset=Probe=ifo[1]
    Conversion=2.5
    FrequencyOffset+=-Conversion
    step=ifo[0]#*1000 - already in ns , no factor of 1000 needed compared to ftmw++ fids
    samples=10**6 #just be bring JUGData on a similar scale as FTMW dat file.
    a=np.array(nums,dtype=np.float64)
    a=a/samples*1000#conversion to mV
    counts=len(a)//2
    side1=a
    side2=a*0

    combined=side1+1.0j*side2
    combined=np.pad(combined, (0, len(combined)*pf), 'constant')

    x=np.linspace(0,len(combined)-1,len(combined))
    x=x*step/1000

    xf = fftfreq(len(combined), step*10**-9)#[1:-1]
    xf=xf/10**6+FrequencyOffset

    yfc=fft(combined)/counts

    yfc=np.abs(yfc)
    yfc=yfc*ic

    idx=np.argsort(xf)
    xf=xf[idx]
    yfc=yfc[idx]

    if writespectrum==1:
        output=open(outname+".txt",'w')
        for j in range(0,len(xf)):
            output.write("{:23.15e} {:23.15e}\n".format(xf[j],yfc[j]))
        output.close()

    if writefid==1:
     output2=open(outname+"_fid.txt",'w')
     xshort=x[:len(x)//(pf+1)]
     for j in range(0,len(xshort)):
        output2.write("{:23.15e} {:23.15e} {:23.15e}\n".format(xshort[j],np.real(side1[j]),np.real(side2[j])))
     output2.close()

    return xf,yfc # returns frequency and intensity
