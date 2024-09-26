import subprocess
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import json
from pathlib import Path

calc_msdmj = True #if true: msd_hdf5.py needed
only_calc_msdmj = False #if true: msd_hdf5.py needed, maxtimes not needed
refit = True #if true: msdMJ_fit.py msdMJ_inp.json needed
plot = True
show_plots = False

endtimes = [50] #in ns
starttime = 0
maxtimes = [600] #in ps
dt = 0.5 #ps

with open("../out/npt_7.out", "r") as f:
    second_last_line = f.readlines()[-2]
    boxl = round(float(second_last_line.split()[-1]),2)
print(f"Boxlength used is {boxl} A")

fitfunctions = {}
cond_dict = {}

if only_calc_msdmj:
    for time in endtimes:
        print(f"At time {time}")
        os.system(f"python msd_hdf5.py {str(starttime)} {str(time*1000)}")
    quit()

for time in endtimes:
    if calc_msdmj:
        os.system(f"python msd_hdf5.py {str(starttime)} {str(time*1000)}")
    fname = "msdmj/msdMJ_1-100_5_total.dat"
    if not os.path.isfile(fname):
        fname = f"msdmj/msdMJ_{str(starttime)}-{str(int(time*1000))}_total.dat"
    if not os.path.isfile(fname):
        print(f"Input file {fname} not found!")
        sys.exit()
    orig_data = np.loadtxt(fname)
    fitfunctions[time] = [orig_data]
    cond_dict[time] = []
    for maxtime in maxtimes:
        newfile = f"msdMJ_inp_{str(starttime)}-{str(int(time*1000))}_{str(maxtime)}.json"
        fitname = f"msdmj/msdMJ_fit_{str(starttime)}-{str(time*1000)}_total.dat"
        #fit again
        if refit:
            #copy input json
            os.system(f"cp msdMJ_inp.json {newfile}")
            # replace input file name
            replace_string = f's~"data":.*~"data":      "{fname}",~'
            os.system(f"sed -i '{replace_string}' {newfile}")
            # replace fit filename
            replace_string = f's~"residuals":.*~"residuals": "{fitname}",~'
            os.system(f"sed -i '{replace_string}' {newfile}")
            # replace boxl
            replace_string = f's~"boxl":.*~"boxl": "{boxl}",~'
            os.system(f"sed -i '{replace_string}' {newfile}")
            # replace maxtime
            replace_string = f's~"maxtime":.*~"maxtime": {maxtime},~'
            os.system(f"sed -i '{replace_string}' {newfile}")

            #os.system(f"python msdMJ_fit.py {newfile}")
            output = subprocess.check_output(f"python msdMJ_fit.py {newfile}", shell=True, text=True)
            #output = subprocess.check_output(f"python msdMJ_fit.py {newfile}", shell=True, text=True)
            cond = output.split("Conductivity is: ")[1].split(" ")[0]
        # do not fit again, use data from earlier cretaed fit files, basically just plot
        else:
            with open(newfile) as infile:
                data = json.load(infile)
                slope = float(data["line"][0]["slope"]["value"])
                boxl = float(data["boxl"])
            temp = 300.
            cond = ( 1.602176634 * slope ) / (6 * 8.617332478 * 10**-8 * temp * boxl**3 )
            
        cond_dict[time].append(round(float(cond)*10,2))
        print(f"{maxtime = } ps, {time = } ns")
        print("Conductivity is: [S/m]", cond)

        fit = np.loadtxt(fitname[:-4] + f"_{maxtime}.dat")
        fitfunctions[time].append(fit)


with open("conductivities.txt", "w") as f:
    f.write(f"     {'   '.join(str(i) for i in endtimes)}\n")
    for idx, maxtime in enumerate(maxtimes):
        c = [str(value[idx]) for value in cond_dict.values()]
        f.write(f"{maxtime} {' '.join(c)}\n")


if plot:
    savedir = "./plots"
    Path(savedir).mkdir(parents=True, exist_ok=True)
    # different cond plots
    evenly_spaced_intervals = np.linspace(0,1, len(endtimes))
    colors = [matplotlib.cm.rainbow(x) for x in evenly_spaced_intervals]
    for pos, time in enumerate(endtimes):
        color = colors[pos]
        plt.plot(maxtimes, cond_dict[time][:], color=color, marker="o", label=f"{time} ns")
    plt.axhline(y=4., color="r", linestyle="--")
    plt.legend(loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True)
    plt.ylabel("Conductivity (mS/cm)")
    plt.xlabel("Fittimes (ps)")
    fig = plt.gcf()
    if show_plots:
        plt.show()
    fig.savefig(f"{savedir}/cond_different_maxtimes.pdf")
    plt.close()


    #comparison plot for total times
    plt.axis([0,np.max(endtimes)*1000, 0, np.max(fitfunctions[int(np.max(endtimes))][0][:,1])])
    for time in endtimes:
        plt.plot(fitfunctions[time][0][:,0],fitfunctions[time][0][:,1], label=f"{time}ns,orig", linewidth=1.0)
    plt.ylabel("msdmj")
    plt.xlabel("time ps")
    plt.legend()
    plt.axes([0.55,0.16,0.34,0.34])
    plt.yticks([])
    for time in endtimes:
        plt.plot(fitfunctions[time][0][0:2000*int(1/dt),0],fitfunctions[time][0][0:2000*int(1/dt),1], label=f"{time},orig", linewidth=1.0)

    fig = plt.gcf()
    if show_plots:
        plt.show()
    fig.savefig(f"{savedir}/comparison_endtimes.pdf")
    plt.close()
    
    
    # one endtime, different fits
    for time in endtimes:
        plt.plot(fitfunctions[time][0][0:maxtimes[-1]*int(1/dt),0],fitfunctions[time][0][0:maxtimes[-1]*int(1/dt),1], "r", label=f"{time}ns,orig", linewidth=1.0)
        plt.axis([0,maxtimes[-1]*1.1, 0, np.max(fitfunctions[time][0][0:maxtimes[-1]*int(1/dt),1])*1.1])
        for idx, maxtime in enumerate(maxtimes,1):
            plt.plot(fitfunctions[time][0][0:maxtime*int(1/dt),0],fitfunctions[time][idx][0:maxtime*int(1/dt),2], label=f"{time}ns,{maxtime}ps,{cond_dict[time][idx-1]}mS/cm", linewidth=0.5)
        plt.ylabel("msdmj")
        plt.xlabel("time ps")
        plt.legend()
    
        a = plt.axes([0.55,0.15,0.35,0.35])
        plt.yticks([])
        plt.plot(fitfunctions[time][0][0:100*int(1/dt),0],fitfunctions[time][0][0:100*int(1/dt),1], "r", label=f"{time},orig", linewidth=1.0)
        for idx, maxtime in enumerate(maxtimes,1):
            plt.plot(fitfunctions[time][0][0:100*int(1/dt),0],fitfunctions[time][idx][0:100*int(1/dt),2], label=f"{time},{maxtime}", linewidth=0.5)
        fig = plt.gcf()
        if show_plots:
            plt.show()
        fig.savefig(f"{savedir}/comparison_endtime{time}.pdf")
        plt.close()
    
    
