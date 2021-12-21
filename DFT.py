# DFT
# DSP Lab 2
# 1911406

import math, time
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.font_manager import fontManager


def plot_dft(x,y,z,Fs,N,t):
    #    plt.rcParams['font.sans-serif']=['SimHei']
    #    plt.rcParams['axes.unicode_minus']=False
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(12, 18))
    line1, = axes[0].plot(x, 'b')  # the source signal
    axes[0].set_xlabel(u'sample dots')
    axes[0].set_ylabel(u'amplitude')
    axes[0].set_title(u'the source signal')
    line2, = axes[1].plot([Fs/N*i for i in range(int(N/2))], [abs(com_y) for com_y in y[0:int(N/2)]],'r')
    axes[1].set_xlabel(u'frequency/Hz')
    axes[1].set_ylabel(u'amplitude')
    axes[1].set_title(u'the DFT signal amplitude')
    line3, = axes[2].plot([Fs/N*i for i in range(int(N/2))],
                          [-math.atan(com_y.imag/com_y.real) if abs(com_y)>0.1 else 0 for com_y in y[0:int(N/2)]],'r')
    axes[2].set_xlabel(u'frequency/Hz')
    axes[2].set_ylabel(u'phase')
    axes[2].set_title(u'the DFT signal phase')
    line4, = axes[3].plot(z,'b')  # the source signal
    axes[3].set_xlabel(u'sample dots')
    axes[3].set_ylabel(u'amplitude')
    axes[3].set_title(u'the recovered signal by IDFT')
    file_name = 'DFT_simulation_' + str(t) + '.png'
    fig.savefig(file_name,dpi=500,bbox_inches='tight')
    plt.close()


def plot_time(dft_x, my_t, sys_t):
    l = len(dft_x)
    plt.figure()
    plt.title('Time Complexity Comparison')
    plt.show()


def plot_signal(x, xn, lbx, lby):  # given coord x when plot a signal
    plt.figure()
    plt.xlabel(lbx)
    plt.ylabel(lby)
    plt.plot(x,xn)
    plt.show()

if __name__ == "__main__":
    Fs = 8000  # sample rate
    Nli = [100, 200, 300, 500, 700, 800]
    dft_x = []  # coord x for my dft
    dftshift_x = [i for i in range(-400, 400)]  # coord x for my fftshift
    #print(dftshift_x)
    my_dft_tli = []  # running times of my DFT
    sys_dft_tli = []  # running times of sys DFT
    for p in range(len(Nli)):
        # in condition that Fs>=2*fi
        f1 = 1000  # 1st component in signal
        x1 = np.array([np.cos(2*math.pi*f1/Fs*n+math.pi/8) for n in range(Nli[p])])
        f2 = 2000  # 2nd component in signal
        x2 = np.array([np.cos(2*math.pi*f2/Fs*n-math.pi/4)*2 for n in range(Nli[p])])
        f3 = 3000  # 3rd component in signal
        x3 = np.array([np.cos(2*math.pi*f3/Fs*n+math.pi/2)*3 for n in range(Nli[p])])
        f4 = 3000
        x4 = np.array([np.cos(2*math.pi*f4/Fs*n)+0j for n in range(Nli[p])])  # even signal
        x5 = np.array([np.sin(2*math.pi*f4/Fs*n)*1j for n in range(Nli[p])])  # odd signal
        x = np.array([x1[i]+x2[i]+x3[i] for i in range(Nli[p])])  # source signal for prob1&2

        # prob1, 2&3
        y = []  # the result of (my)DFT
        ts_my_dft = time.time()  # timing
        for k in range(Nli[p]):
            basis = [complex(math.cos(2*math.pi/Nli[p]*k*n),math.sin(2*math.pi/Nli[p]*k*n)) for n in range(Nli[p])]
            y.append(np.dot(x,np.transpose(basis)))
        te_my_dft = time.time()
        
        #test
        #print(y)

        dft_x.append(Nli[p])  # numerical magnitudes
        my_dft_tli.append(te_my_dft-ts_my_dft)  # running time in seconds
        # contain the basis exp(j*2*pi/N*k*n) and the projection weight

        ts_sys_dft = time.time()
        sys_y = np.fft.fft(x)  # fft func in numpy
        te_sys_dft = time.time()

        sys_dft_tli.append(te_sys_dft-ts_sys_dft)

        z = []  # the result of IDFT
        for k in range(Nli[p]):
            basis = [complex(math.cos(2*math.pi/Nli[p]*k*n), -math.sin(2*math.pi/Nli[p]*k*n)) for n in range(Nli[p])]
            z.append(np.dot(y, np.transpose(basis))/Nli[p])

        # when input 800 points
        # show fftshift
        #if Nli[p] == 800:
        if Nli[p] == 800:
            plt.figure()
            plt.xlim(-400,400)
            plt.plot(range(-400,400), sys_y)
            plt.xlim(-Nli[p]/2,Nli[p]/2)
            plt.title('myfftshift')
            plt.show()

            # prob4
            # linearity
            alpha1 = 2+0j
            alpha2 = 3+0j
            x6 = alpha1*x4+alpha2*x5
            # plot_signal([_ for _ in range(len(x6))], x6, 't', 'amp')  # show x6
            fft_x4 = np.fft.fft(x4)
            fft_x5 = np.fft.fft(x5)
            fft_sum = alpha1*fft_x4+alpha2*fft_x5
            fft_x6 = np.fft.fft(x6)

            res4 = fft_sum-fft_x6
            flag = True
            for i in range(len(res4)):
                if abs(res4[i]) >= 1+0j:
                    flag = False
                    break
            if flag:
                print('Linearity Proved')  # all approx to 0

            # Parity
            x7 = np.array([np.sin(2*math.pi*f4/Fs*n) for n in range(Nli[p])])  # odd
            fft_1 = np.fft.fft(x7)  # shall be imag
            fft_2 = np.fft.fft(x4)  # shall be real
            flag = True
            for i in range(len(fft_2)):
                if abs(fft_2[i].imag) >= 1:
                    flag = False
                    break
            # print(len(fft_2), len(fft_1))
            for i in range(len(fft_1)):
                if abs(fft_1[i].real) >= 1:
                    flag = False
                    break
            if flag:
                print('Parity Proved')  # all approx to 0
            
            print("problem5: mean filter\n\n")
            # prob5: mean filter
            # time domain
            h = np.array([0.25,0.25,0.25,0.25])
            # print(np.fft.fft(h))
            x = np.array([0]*3+list(x4)+[0]*3)  # zero padding
            xi = np.array([0]*3+list(x4)+[0]*3)
            l = len(x)-3  # 0-
            for i in range(l):
                x[i+3] = x[i+3]*h[0]+x[i+2]*h[1]+x[i+1]*h[2]+x[i]*h[3]  # convolution

            # freq domain
            hf = [0,0,0]+[1,0,0,0]*200+[0,0,0]
            # hfi = np.array(np.fft.ifft(hf))
            xf = np.fft.fft(xi)*hf
            xf = np.fft.ifft(xf)
            # xf = np.array([0]*3+list(xf)+[0]*3)

            plt.figure()
            plt.xlim(0,800)
            plt.plot(range(len(x)), x, label="x")
            plt.plot(range(len(x)), xf, label="xf")
            plt.legend()
            plt.show()
        
    plt.figure()
    plt.plot(dft_x, my_dft_tli, label="my dft")
    plt.plot(dft_x, sys_dft_tli, label="sys fft")
    plt.legend()
    plt.title("Time Complexity Comparison")
    plt.show()



# 0 1 2 3 4 5 6 7 8 9 10, len=5