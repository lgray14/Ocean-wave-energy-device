import matplotlib.pyplot as plt
import numpy as np
import re

pattern = r" Mean calculated power: ([0-9.]+)W"
amps = []
avgs = []

amp15 = []
with open("power_trials/powertrial1-5.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp15.append((float(l.split()[3])))
        avg15 = 0
    
    last = str(lines[-1])
    avg15 = re.findall(r'\d+\.{0,1}\d*', last)
    avg15 = float(avg15[0])
    amps.append(1.5)
    avgs.append(avg15)
# print(amp15, avg15)

amp16 = []
with open("power_trials/powertrial1-6A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp16.append((float(l.split()[3])))
    
    last = str(lines[-1])
    avg16 = re.findall(r'\d+\.{0,1}\d*', last)
    avg16 = float(avg16[0])
    amps.append(1.6)
    avgs.append(avg16)

amp17 = []
with open("power_trials/powertrial1-7A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp17.append((float(l.split()[3])))

    last = str(lines[-1])
    avg17 = re.findall(r'\d+\.{0,1}\d*', last)
    avg17 = float(avg17[0])
    amps.append(1.7)
    avgs.append(avg17)

amp18 = []
with open("power_trials/powertrial1-8.txt", "r") as f:
    lines = f.readlines()
    for l in lines[2:-2]:
        amp18.append((float(l.split()[3])))
    
    last = str(lines[-1])
    avg18 = re.findall(r'\d+\.{0,1}\d*', last)
    avg18 = float(avg18[0])
    amps.append(1.8)
    avgs.append(avg18)

amp19 = []
with open("power_trials/powertrial1-9A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp19.append((float(l.split()[3])))
    last = str(lines[-1])
    avg19 = re.findall(r'\d+\.{0,1}\d*', last)
    avg19 = float(avg19[0])
    amps.append(1.9)
    avgs.append(avg19)

amp20 = []
with open("power_trials/powertrial2-0A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp20.append((float(l.split()[3])))
    last = str(lines[-1])
    avg20 = re.findall(r'\d+\.{0,1}\d*', last)
    avg20 = float(avg20[0])
    amps.append(2.0)
    avgs.append(avg20)

amp21 = []
with open("power_trials/powertrial2-1A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp21.append((float(l.split()[3])))
    last = str(lines[-1])
    avg21 = re.findall(r'\d+\.{0,1}\d*', last)
    avg21 = float(avg21[0])
    amps.append(2.1)
    avgs.append(avg21)

amp22 = []
with open("power_trials/powertrial2-2A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp22.append((float(l.split()[3])))
    last = str(lines[-1])
    avg22 = re.findall(r'\d+\.{0,1}\d*', last)
    avg22 = float(avg22[0])
    amps.append(2.2)
    avgs.append(avg22)

amp23 = []
with open("power_trials/powertrial2-3A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp23.append((float(l.split()[3])))
    last = str(lines[-1])
    avg23 = re.findall(r'\d+\.{0,1}\d*', last)
    avg23 = float(avg23[0])
    amps.append(2.3)
    avgs.append(avg23)

amp24 = []
with open("power_trials/powertrial2-4A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp24.append((float(l.split()[3])))
    last = str(lines[-1])
    avg24 = re.findall(r'\d+\.{0,1}\d*', last)
    avg24 = float(avg24[0])
    amps.append(2.4)
    avgs.append(avg24)

amp25 = []
with open("power_trials/powertrial2-5A.txt", "r") as f:
    lines = f.readlines()
    for l in lines[1:-2]:
        amp25.append((float(l.split()[3])))
    last = str(lines[-1])
    avg25 = re.findall(r'\d+\.{0,1}\d*', last)
    avg25 = float(avg25[0])
    amps.append(2.5)
    avgs.append(avg25)

plt.scatter(1.5*(np.ones(len(amp15))), amp15, c='#b22222')
plt.scatter(1.6*(np.ones(len(amp16))), amp16, c='#b22222')
plt.scatter(1.7*(np.ones(len(amp17))), amp17, c='#b22222')
# print(amp17)
plt.scatter(1.8*(np.ones(len(amp18))), amp18, c='#b22222')
plt.scatter(1.9*(np.ones(len(amp19))), amp19, c='#b22222')
plt.scatter(2.0*(np.ones(len(amp20))), amp20, c='#b22222')
plt.scatter(2.1*(np.ones(len(amp21))), amp21, c='#b22222')
plt.scatter(2.2*(np.ones(len(amp22))), amp22, c='#b22222')
plt.scatter(2.3*(np.ones(len(amp23))), amp23, c='#b22222')
plt.scatter(2.4*(np.ones(len(amp24))), amp24, c='#b22222')
plt.scatter(2.5*(np.ones(len(amp25))), amp25, c='#b22222')

plt.plot(amps, avgs, "-*", label='Mean power')
plt.xlabel("Amperes [A]")
plt.ylabel('Power output of device [W]')
plt.legend()
plt.title('Current vs Power for .9Hz Wave Trial')

plt.show()