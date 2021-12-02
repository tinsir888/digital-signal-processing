import copy

class Sequence:
    start = 0
    tail = 0
    seq = []

    # 初始化信号，可指定下标范围
    def __init__(self, s, t, x):
        self.start = s
        self.tail = t
        self.seq = x

    # 索引
    def __getitem__(self, k):
        if k < self.start or k > self.tail:
            return 0
        else:
            return self.seq[k-self.start]

    # 补零
    def extension(self, s, t):
        if s < self.start:
            self.seq = [0]*(self.start - s) + self.seq
            self.start = s
        if t > self.tail:
            self.seq = self.seq + [0]*(t-self.tail)
            self.tail = t
        return self.seq

    # 延迟和提前，作变换y[n]=x[n+t]，t>0为提前
    def delay(self, t):
        # self.extension(self.start+t, self.tail+t
        self.start += t
        self.tail += t
        return self.seq

    # 反折
    def rev(self):
        self.seq.reverse()
        return self.seq

    # 拉伸/压缩（上/下采样）
    def stretch(self, k, samp='u'):
        if k == 0:
            return self.seq
        li = []
        if samp == 'u':  # 上采样，y[n]=x[n/k](n=mk)，序列总长不变
            if k < 0:
                k = -k
                self.rev()
            for i in range(self.start, self.tail+1):
                if i%k:
                    li += [0]
                else:
                    x = int(i/k)
                    li += [self[x]]
        elif samp == 'd':  # 下采样，y[n]=x[nk]，序列缩短
            if k < 0:
                k = -k
                self.rev()
            for i in range(self.start, self.tail + 1):
                if self.start <= i*k <= self.tail:
                    li += [self[i*k]]
        self.seq = li
        return self.seq

    # 差分累加
    def diff(self, tp='d'):
        if tp == 'd':
            #self.seq = [(x1[i + 1] - x1[i]) for i in range(self.start, self.tail)]
            for i in range(len(self.seq) - 1):
                self.seq[i] = self.seq[i + 1] - self.seq[i]
            self.seq.pop()
            self.tail = self.tail - 1
            return self.seq
        elif tp == 'a':
            sum = 0
            for i in range(len(self.seq)):
                sum += self.seq[i]
                self.seq[i] = sum
            return self.seq

    # 加法
    def __add__(self, other):
        lindex = min(self.start, other.start)
        rindex = max(self.tail, other.tail)
        self.extension(lindex, rindex)
        other.extension(lindex, rindex)
        self.seq = list(map(lambda y: y[0]+y[1], zip(self.seq, other.seq)))
        return self.seq

    # 乘法
    def __mul__(self, other):
        # if type(other) == 'list':  # 调变
            lindex = min(self.start, other.start)
            rindex = max(self.tail, other.tail)
            self.extension(lindex, rindex)
            other.extension(lindex, rindex)
            self.seq = list(map(lambda y: y[0]*y[1], zip(self.seq, other.seq)))
        # elif type(other) == ('int' or 'float' or 'double' or 'long'):  # 数乘
        #     self.seq = [i*other for i in self.seq]
            return self.seq

    # 卷积
    def convolution(self, kern, mode='L'):
        kern.reverse()
        l = len(kern)
        if mode == 'L':
            self.extension(self.start-(l-1), self.tail+(l-1))
            # print(self.seq)
            li = copy.deepcopy(self.seq)
            for i in range(self.start, self.tail-(l-1)+1):
                sum = 0
                for j in range(l):
                    sum += kern[j] * self[i+j]
                print(sum)
                li[i-self.start+l-1]=sum
            self.start += (l-1)
            self.seq = copy.deepcopy(li[1:])
        elif mode == 'C':
            self.seq = self.seq[-l+1:]+self.seq+self.seq[:self.start+l-1]
            self.start -= l-1
            self.tail += l-1
            li = copy.deepcopy(self.seq)
            for i in range(self.start,self.tail-(l-1)+1):
                sum = 0
                for j in range(l):
                    sum += kern[j] * self[i + j]
                li[i - self.start + l - 1] = sum
            self.start+=l-1
            self.seq = copy.deepcopy(li[1:])
        return self.seq

    # 相似性比对
    def comp(self, k, mode='N'):
        li = []
        if mode == 'W':
            sum = 0
            if k == '1':
                for i in range(-4,5):
                    x2.extension(x1.start-i, x1.tail+i)
                    l=len(x2.seq)
                    sum = 0
                    for j in range(self.start, self.tail + 1):
                        sum += x2[j-i] * self[j]
                    li.append(sum)
            if k == '2':
                for i in range(-4,5):
                    x1.extension(x2.start-i, x2.tail+i)
                    l=len(x1.seq)
                    sum = 0
                    for j in range(self.start, self.tail + 1):
                        sum += x1[j-i] * self[j]
                    li.append(sum)
        elif mode == 'N':
            sum = 0
            if k == '1':
                for i in range(-4, 5):
                    x1.norm()
                    x2.norm()
                    x2.extension(x1.start - i, x1.tail + i)
                    l = len(x2.seq)
                    sum = 0
                    for j in range(self.start, self.tail + 1):
                        sum += x2[j - i] * self[j]
                    li.append(sum)
            if k == '2':
                for i in range(-4, 5):
                    x1.norm()
                    x2.norm()
                    x1.extension(x2.start - i, x2.tail + i)
                    l = len(x1.seq)
                    sum = 0
                    for j in range(self.start, self.tail + 1):
                        sum += x1[j - i] * self[j]
                    li.append(sum)
        return li

    # 零均值与归一化
    def norm(self):
        mean=0
        for i in range(len(self.seq)):
            mean += self.seq[i]
        mean /= len(self.seq)
        for i in range(len(self.seq)):
            self.seq[i] -= mean  # 零均值
        mean=0
        for i in range(len(self.seq)):
            mean += (self.seq[i] * self.seq[i])  # 归一化
        import math
        mean = math.sqrt(mean)
        for i in range(len(self.seq)):
            self.seq[i] = self.seq[i] / mean
if __name__ == '__main__':
    print("--------序列计算机v1.0！--------")
    print(">> 说明：本计算机对原始序列x进行变换，但操作需要时也可输入新的序列")
    print(">> 准备就绪，输入'q'退出使用")

    # 初始化x1和x2两个序列
    xs1, xt1 = map(int, input(">> 初始化序列x1，请输入起止下标（以空格分隔）：\n").split())
    l = [int(x) for x in input(">> 请输入序列内容(以空格分隔)：\n").split()]
    x1 = Sequence(xs1, xt1, l)

    xs2, xt2 = map(int, input(">> 初始化序列x2，请输入起止下标（以空格分隔）：\n").split())
    l = [int(x) for x in input(">> 请输入序列内容(以空格分隔)：\n").split()]
    x2 = Sequence(xs2, xt2, l)

    # 选择运算
    c = input(">> 输入以下数字之一选择操作：\n\t0：查看某序列内容；1：补零；2：延迟/提前；3：反转；4：拉伸/压缩；5：一阶差分/累加；6：序列间加法；7：序列间乘法；8：卷积；9：相似性比对\n")
    while c != 'q':
        if c == '0':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                n = input(">> 输入下标（输入'w'退出）：")
                while n != 'w':
                    print(x1[int(n)])
                    n = input(">> 输入下标（输入'w'退出）：")
            elif t == '2':
                n = input(">> 输入下标（输入'w'退出）：")
                while n != 'w':
                    print(x2[int(n)])
                    n = input(">> 输入下标（输入'w'退出）：")
        elif c == '1':
            t = input(">> 输入待补零序列编号（1或2）：")
            if t == '1':
                print(">> 将x1相对于序列x2补零，结果为：")
                print(x1.extension(x2.start, x2.tail))
            elif t == '2':
                print(">> 将x2相对于序列x1补零，结果为：")
                print(x2.extension(x1.start, x1.tail))
        elif c == '2':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                k = int(input(">> 输入延迟时长（正数时延迟）："))
                print(x1.delay(k))
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                k = int(input(">> 输入延迟时长（正数时延迟）："))
                print(x2.delay(k))
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        elif c == '3':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                print(x1.rev())
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                print(x2.rev())
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        elif c == '4':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                k = int(input(">> 输入采样间距："))
                s = input(">> 输入'u'/'d'选择上/下采样：")
                print(x1.stretch(k,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                k = int(input(">> 输入采样间距："))
                s = input(">> 输入'u'/'d'选择上/下采样：")
                print(x2.stretch(k,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        elif c == '5':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                s = input(">> 输入'd'/'a'选择差分/累加：")
                print(x1.diff(s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                s = input(">> 输入'd'/'a'选择差分/累加：")
                print(x2.diff(s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        elif c == '6':
            print(x1+x2)
            print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
        elif c == '7':
            print(x1*x2)
            print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
        elif c == '8':
            print(">> 初始化卷积核")
            li = [int(x) for x in input(">> 输入卷积核内容(以空格分隔)：\n").split()]
            print(li)

            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                s = input(">> 输入'L'/'C'选择线性/循环卷积：")
                print(x1.convolution(li,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                s = input(">> 输入'L'/'C'选择线性/循环卷积：")
                print(x2.convolution(li,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        elif c == '9':
            t = input(">> 输入序列编号（1或2）：")
            if t == '1':
                s = input(">> 输入'W'/'N'选择滑动窗/归一化：")
                print(x1.comp(t,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x1.start, "e": x1.tail})
            elif t == '2':
                s = input(">> 输入'W'/'N'选择滑动窗/归一化：")
                print(x2.comp(t,s))
                print("index: start=%(s)d, end=%(e)d" % {"s": x2.start, "e": x2.tail})
        c = input(
            ">> 输入以下数字之一选择操作：\n\t0：初始化序列以及查看该序列内容；1：补零；2：延迟/提前；3：反转；4：拉伸/压缩；5：一阶差分/累加；6：序列间加法；7：序列间乘法；8：卷积；9：相似性比对；q：退出使用\n")
    print("--------感谢使用！--------")
