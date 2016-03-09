import pickle
from numpy import exp, array, pi, sqrt, log, unique, mean, std
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


with open('data/nondims.pickle', 'rb') as f:
    data = pickle.load(f)


def xAxis(nondimsArray, aRe, aSt, aWe, aCa, C1, C2):
    Re = array(nondimsArray['Re'])
    St = array(nondimsArray['St'])
    We = array(nondimsArray['We'])
    Ca = array(nondimsArray['Ca'])

    #return Re
    # c1
    #return Re / Ca

    # c4
    #return Re ** aRe / We ** aWe / Ca ** aCa 

    # c3
    #return We ** aWe / Re ** aRe

    # c2
    #return We * Re ** aRe / Ca ** aCa * St ** aSt

def dependency(nondimsArray, A, B, aRe, aSt, aWe, aCa, C1, C2):
    x = xAxis(nondimsArray, aRe, aSt, aWe, aCa, C1, C2)
    return A * x + B

target = 'c1'
Ci = array(data[target])
#Ci = array(data['c4']) / 1e12
popt, pcov = curve_fit(dependency, data, Ci)
print(popt)

x = xAxis(data, *popt[2:])

splitData = dict()
keys = data['key']
uniqueKeys = unique(keys)

for uK in uniqueKeys:
    splitData[uK] = dict()
    for key in data.keys():
        splitData[uK][key] = []

for i in range(len(keys)):
    name = data['key'][i]
    for key in data.keys():
        splitData[name][key].append(data[key][i])

fig = plt.figure()


#plt.plot(x, Ci, 'x')
plt.plot(x, dependency(data, *popt), '.')
for key, value in splitData.iteritems():
    x = xAxis(value, *popt[2:])
    plt.plot(x, value[target], 'o', label=key)
#plt.yscale('log')
#plt.xscale('log')

plt.legend(loc='best')
plt.show()

means = dict()
stds = dict()

for i in range(4):
    key = 'c' + repr(i+1)
    means[key] = mean(data[key])
    stds[key] = std(data[key])

print means[target], stds[target]
