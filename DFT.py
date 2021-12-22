import numpy as np
import math
import copy
import timeit
import random
import matplotlib.pyplot as plt
#import cv2

def fft1(src, dst = None):
	
	l = len(src)
	n = int(math.log(l,2))
 
	bfsize = np.zeros((l), dtype = "complex")
 
	for i in range(n + 1):
		if i == 0:
			for j in range(l):	
				bfsize[j] = src[Dec2Bin_Inverse2Dec(j, n)]
		else:
			tmp = copy.copy(bfsize)
			for j in range(l):
				pos = j%(pow(2,i))
				if pos < pow(2, i - 1):
					bfsize[j] = tmp[j] + tmp[j + pow(2, i - 1)] * np.exp(complex(0, -2*np.pi*pos/pow(2,i)))
					bfsize[j + pow(2, i - 1)] = tmp[j] - tmp[j + pow(2, i - 1)] * np.exp(complex(0, -2*np.pi*pos/(pow(2,i))))
	return bfsize
 
def ifft1(src):
 
	for i in range(len(src)):
		src[i] = complex(src[i].real, -src[i].imag)
 
	res = fft1(src)
 
	for i in range(len(res)):
		res[i] = complex(res[i].real, -res[i].imag)
 
	return res/len(res)
 
def Dec2Bin_Inverse2Dec(n, m):
	
	b = bin(n)[2:]
	if len(b) != m:
		b = "0"*(m-len(b)) + b
	b = b[::-1]
	return int(b,2)

def fft1shift(x):
    return np.roll(x, int(len(x) / 2))

def cov_in_time(a, b):
    c = np.tile(b, 3)
    ma = max(len(a), len(b))
    cov = np.convolve(a, c)
    r = cov[ma:2 * ma]
    return r

def get_lowpass_mask(l, k):
    a = np.zeros_like(l, np.double)
    c = int(len(l) / 2)
    a[c - k:c + k] = 1
    return a

def filt_by_frequency(s, mask):
    f = np.fft.fft(s)
    r = f * mask
    return np.fft.ifft(r)

def filt_by_timedomain(b, mask):
    ift = np.fft.ifft(mask)
    return cov_in_time(b, ift)

def filter():
    l = 50
    k = 10
    rang = np.arange(l)
    s = np.sin(rang)+2*np.cos(rang)
    mask = get_lowpass_mask(s, k)
    res_byf = filt_by_frequency(s, mask)
    res_byt = filt_by_timedomain(s, mask)
    fig, subs = plt.subplots(2, 2)
    subs[0][0].plot(s)
    subs[0][1].plot(mask)
    subs[1][0].plot(res_byf)
    subs[1][1].plot(res_byt)
    plt.show()


print('任务一：手工实现fft和ifft并进行比对')
src = np.random.randint(20,size=8)
print('自己实现的fft:')
print(fft1(src))
print('自己实现的ifft:')
print(ifft1(fft1(src)))
print('np里面的fft:')
hhh=np.fft.fft(src)
print(hhh)
print('np里面的ifft:')
print(np.fft.ifft(hhh))
print()



print('任务二：手写傅里叶变换和快速傅里叶变换的效率比较')
data = []
for i in [2**x for x in range(14)]:
	data.append((i, timeit.timeit('fft1(x)','from __main__ import fft1;import random; x=list(range(%d)); random.shuffle(x);' %i, number=10)))
plt.figure()
plt.scatter(*zip(*data))
plt.xlabel('n')
plt.ylabel('time')
print('手写傅里叶变换数据规模和耗时的关系如下：')
plt.show()

data = []
for i in [2**x for x in range(14)]:
	data.append((i, timeit.timeit('np.fft.fft(x)','import numpy as np;import random; x=list(range(%d)); random.shuffle(x);' %i, number=10)))
plt.figure()
plt.scatter(*zip(*data))
plt.xlabel('n')
plt.ylabel('time')
print('快速傅里叶变换数据规模和耗时的关系如下：')
plt.show()

data = []
for i in [2**x for x in range(14)]:
	data.append((i, timeit.timeit('ifft1(x)','from __main__ import ifft1;import random; x=list(range(%d)); random.shuffle(x);' %i, number=10)))
plt.figure()
plt.scatter(*zip(*data))
plt.xlabel('n')
plt.ylabel('time')
print('手写傅里叶逆变换数据规模和耗时的关系如下：')
plt.show()

data = []
for i in [2**x for x in range(14)]:
	data.append((i, timeit.timeit('np.fft.ifft(x)','import numpy as np;import random; x=list(range(%d)); random.shuffle(x);' %i, number=10)))
plt.figure()
plt.scatter(*zip(*data))
plt.xlabel('n')
plt.ylabel('time')
print('快速傅里叶逆变换数据规模和耗时的关系如下：')
plt.show()
print()
#


print('任务三')
print('经过fft：')
print(hhh)
fs1=fft1shift(hhh)
fs2=np.fft.fftshift(hhh)
print('再经过我的fftshift：')
print(fs1)
print('或者经过np内置的fftshift：')
print(fs2)
plt.plot(hhh,"ob")
plt.show()
plt.plot(fs1,"ob")
plt.show()
plt.plot(fs2,"ob")
plt.show()
print()
#

'''
print('任务四')
print('验证线性')
x=np.array([1,2,3,4])
y=np.array([2,3,4,5])
xy=np.array(x+y)
print(x)
print(y)
print(xy)
fx=np.fft.fft(x)
fy=np.fft.fft(y)
fxy=np.fft.fft(xy)
print(fx)
print(fy)
print(fxy)
print(np.fft.ifft(fx))
print(np.fft.ifft(fy))
print(np.fft.ifft(fxy))
print('验证对称性')
x=np.array([1,2,3,4])
print(x)
fx=np.fft.fft(x)
fy=4*np.fft.ifft(x)
print(fx)
print(fy)
print()
#
'''
print('任务四：滤波器')
filter()
