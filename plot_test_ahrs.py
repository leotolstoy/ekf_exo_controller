import matplotlib.pyplot as plt
import numpy as np
import csv
from attitudeEstimatorEKFFuncs import extractEulerAngles_new

with open('test_ahrs.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    Rs = []
    R_foots = []
    ts=[]
    accXs=[]
    accYs=[]
    accZs=[]

    rolls=[]
    pitchs=[]
    yaws=[]

    for i, row in enumerate(reader):
        if i==0:
            print (", ".join(row))
        else:
            

            X = np.array([float(x) for x in row])
            # print(X)
            ts.append(X[0])
            R = np.array(X)[1:10].reshape((3,3))
            Rs.append(R)

            R_foot = np.array(X)[10:19].reshape((3,3))
            R_foots.append(R_foot)
            roll, pitch, yaw = extractEulerAngles_new(R_foot)
            accXs.append(X[19])
            accYs.append(X[20])
            accZs.append(X[21])
            rolls.append(roll*180/np.pi)
            pitchs.append(pitch*180/np.pi)
            yaws.append(yaw*180/np.pi)

ts = ts-ts[0]

accXs=np.array(accXs)
accYs=np.array(accYs)
accZs=np.array(accZs)


accNorm = np.sqrt(accXs**2 + accYs**2 + accZs**2)
Rs=np.array(Rs)
R_foots=np.array(R_foots)

plt.plot(ts, R_foots[:,0,0], label="R11 foot")
plt.plot(ts, R_foots[:,1,1], label="R22")
plt.plot(ts, R_foots[:,2,2], label="R33")

plt.plot(ts, R_foots[:,0,1], label="R12")
plt.plot(ts, R_foots[:,0,2], label="R13")
plt.plot(ts, R_foots[:,1,2], label="R23")
plt.legend()


fig, axs = plt.subplots()
axs.plot(ts,accXs,label='accXs')
axs.plot(ts,accYs,label='accYs')
axs.plot(ts,accZs,label='accZs')
# axs.plot(ts,accNorm,'k',label='accNorm')
plt.legend()


fig, axs = plt.subplots()
axs.plot(ts,rolls,label='rolls')
axs.plot(ts,pitchs,label='pitchs')
axs.plot(ts,yaws,label='yaws')
# axs.plot(ts,accNorm,'k',label='accNorm')
plt.legend()











plt.show()
