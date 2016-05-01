#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import pylab

t = np.arange(0.0, 2.0, 0.01)
s1 = np.sin(2*np.pi*t)
s2 = np.sin(4*np.pi*t)

plt.figure(1)
plt.subplot(211)
plt.plot(t, s1)
plt.subplot(212)
plt.plot(t, 2*s1)

plt.figure(2)
plt.plot(t, s2)

# now switch back to figure 1 and make some changes
plt.figure(1)
plt.subplot(211)
plt.plot(t, s2, 'gs')
ax = plt.gca()
ax.set_xticklabels([])


with open("a.dat") as f:
    data = f.read()

    data = data.split()
    x1 =[]
    y1 =[]
    x1.append(float(data[0]))
    y1.append(float(data[1]))
#print y1

print "OK"
xv = np.array(x1)
yv = np.array(y1)
print xv


    
datalist = [ ( pylab.loadtxt("a.dat"), label) for filename, label in label ]
    
for data, label in datalist:
    pylab.plot( data[:,0], data[:,1], label=label )
    
    pylab.legend()
    pylab.title("Title of Plot")
    pylab.xlabel("X Axis Label")
    pylab.ylabel("Y Axis Label")

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Plot title...")
ax1.set_xlabel('your x label..')
ax1.set_ylabel('your y label...')
pylab.xlim([1,5])
pylab.ylim([0.2,1.0])

ax1.plot(xv,yv, c='r', label='the data')

leg = ax1.legend()

plt.show()