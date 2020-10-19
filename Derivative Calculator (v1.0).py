#differential equations

def isnum(var):
    if str(type(var)) == "<class '__main__.f.const'>":
        return True
    else: return False

def parens(s):
    if isnum(s):
        if type(s.c) == int or type(s.c) == float:
            return str(s.c)
        else:
            return '(%s)' % s.c
    else:
        if str(type(s)) == "<class '__main__.f.var'>":
            return s.toStr()
        else:
            return '(%s)' % s.toStr()
        
def tofconst(var):
    if type(var) == int or type(var) == float or str(type(var)) == "<class 'fractions.Fraction'>":
        return f.const(var)
    else:
        return var


class f:
    class var:
        def __init__(self, varname='x'):
            self.name = varname

        def der(self):
            return 1

        def toStr(self):
            return self.name

        def evalu(self, val):
            return val

    class const:
        def __init__(self, c):
            self.c=c

        def der(self):
            return f.const(0)

        def toStr(self):
            return str(self.c)

        def evalu(self,val):
            return self.c
    
        
    class sum:
        def __init__(self, *args):
            if type(args[0]) == list: self.p = args[0]
            else: self.p = list(args)
            self.p[:] = [tofconst(x) for x in self.p]

        def der(self):
            ret = []
            for x in self.p:
                if not isnum(x):
                    ret.append(x.der())
            return f.sum(ret)

        def toStr(self):
            return '+'.join([x.toStr() for x in self.p])

        def evalu(self, val):
            return sum([x.evalu(val) for x in self.p])
            

    class prod:
        def __init__(self, *args):
            if type(args[0])== list: self.p = args[0]
            else: self.p = list(args)
            self.p[:] = [tofconst(x) for x in self.p]

        def der(self):
            funcs = []
            mult = 1
            for x in self.p:
                if isnum(x):
                    mult*=x.c
                else:
                    funcs.append(x)
            if len(funcs) == 0:
                return 0
            else:
                ret = []
                for n,func in enumerate(funcs):
                    ret.append(f.prod([mult,func.der()]+funcs[:n]+funcs[n+1:]))
                return f.sum(ret)

        def toStr(self):
            funcs = []
            mult = 1
            for x in self.p:
                if isnum(x):
                    mult*=x.c
                else:
                    funcs.append(parens(x))
            if mult == 1 and len(funcs)>0:
                return '*'.join(funcs)
            else:
                return '*'.join([str(mult)]+funcs)

        def evalu(self,val):
            ret = 1
            for x in self.p:
                ret *= x.evalu(val)
            return ret

    class exp:
        def __init__(self, power):
            self.p = tofconst(power)

        def der(self):
            return f.prod(self.p.der(),f.exp(self.p))

        def toStr(self):
            return 'e^%s' % parens(self.p)

        def evalu(self,val):
            return exp(self.p.evalu(val))

    class ln:
        def __init__(self, fval):
            self.v = tofconst(fval)

        def der(self):
            return f.fract(self.v.der(),self.v)

        def toStr(self):
            return 'ln(%s)' % self.v.toStr()

        def evalu(self,val):
            return log(self.v.evalu(val))

    class logb:
        def __init__(self,base,fval):
            self.v = tofconst(fval)
            self.b = tofconst(base)

        def der(self):
            if not isnum(self.b): raise Exception('Variable logarithmic bases not supported')
            return f.fract(self.v.der(),f.prod(f.ln(self.b),self.v))

        def toStr(self):
            return 'log_%s (%s)' % (parens(self.b),self.v.toStr())

        def evalu(self,val):
            return log(self.v.evalu(val))/log(self.b.evalu(val))


    class pow:
        def __init__(self, *args):
            if len(args) == 1 and type(args[0]) == list:
                [self.base,self.power] = args[0]
            elif len(args) == 2:
                [self.base,self.power] = args
            else: raise Exception('Power input wrong')
            self.base = tofconst(self.base)
            self.power = tofconst(self.power)

        def der(self):
            if isnum(self.power) and isnum(self.base):
                return f.const(0)
            elif isnum(self.power):
                return f.prod(self.power,self.base.der(),f.pow(self.base,f.const(self.power.c-1)))
            elif isnum(self.base):
                if self.base.c<0: raise Exception('negative exponential base')
                return f.prod(f.ln(self.base),self.power.der(),f.pow(self.base,self.power))
            else:
                return f.exp(f.prod(self.power,f.ln(self.base))).der()

        def toStr(self):
            return '%s^%s' % (parens(self.base),parens(self.power))

        def evalu(self,val):
            return (self.base.evalu(val))**(self.power.evalu(val))

    class fract:
        def __init__(self, numerator, denominator):
            self.n = tofconst(numerator)
            self.d = tofconst(denominator)

        def der(self):
            return f.fract(f.sum(f.prod(self.n.der(),self.d),f.prod(-1,self.d.der(),self.n)),f.pow(self.d,2))

        def toStr(self):
            return '%s/%s' % (parens(self.n),parens(self.d))

        def evalu(self,val):
            return self.n.evalu(val)/self.d.evalu(val)

    class sqrt:
        def __init__(self, root):
            self.r = tofconst(root)
        def der(self):
            return f.pow(self.r,frac(1,2)).der()
        def toStr(self):
            return '√%s' % parens(self.r)
        def evalu(self,val):
            return (self.r.evalu(val))**0.5

    class abs:
        def __init__(self, fval):
            self.v = tofconst(fval)
        def der(self):
            return f.prod(self.v.der(),f.fract(self.v,f.abs(self.v)))
        def toStr(self):
            return '|%s|' % self.v.toStr()
        def evalu(self,val):
            return abs(self.v.evalu(val))
        
    #Trig functions
    class sin:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(self.v.der(),f.cos(self.v))
        def toStr(self): return 'sin(%s)' % self.v.toStr()
        def evalu(self,val): return sin(self.v.evalu(val))
    class cos:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(-1,self.v.der(),f.sin(self.v))
        def toStr(self): return 'cos(%s)' % self.v.toStr()
        def evalu(self,val): return cos(self.v.evalu(val))
    class tan:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(self.v.der(),f.pow(f.sec(self.v),2))
        def toStr(self): return 'tan(%s)' % self.v.toStr()
        def evalu(self,val): return tan(self.v.evalu(val))
    class csc:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(-1,self.v.der(),f.csc(self.v),f.cot(self.v))
        def toStr(self): return 'csc(%s)' % self.v.toStr()
        def evalu(self,val): return 1/sin(self.v.evalu(val))
    class sec:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(self.v.der(),f.sec(self.v),f.tan(self.v))
        def toStr(self): return 'sec(%s)' % self.v.toStr()
        def evalu(self,val): return 1/cos(self.v.evalu(val))
    class cot:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.prod(-1,self.v.der(),f.pow(f.csc(self.v),2))
        def toStr(self): return 'cot(%s)' % self.v.toStr()
        def evalu(self,val): return cos(self.v.evalu(val))/sin(self.v.evalu(val))
    #inverse trig
    class arcsin:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(self.v.der(),f.sqrt(f.sum(1,f.prod(-1,f.pow(self.v,2)))))
        def toStr(self): return 'arcsin(%s)' % self.v.toStr()
        def evalu(self,val): return asin(self.v.evalu(val))
    class arccos:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(f.prod(-1,self.v.der()),f.sqrt(f.sum(1,f.prod(-1,f.pow(self.v,2)))))
        def toStr(self): return 'arccos(%s)' % self.v.toStr()
        def evalu(self,val): return acos(self.v.evalu(val))
    class arctan:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(self.v.der(),f.sum(1,f.pow(self.v,2)))
        def toStr(self): return 'arctan(%s)' % self.v.toStr()
        def evalu(self,val): return atan(self.v.evalu(val))
    class arccsc:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(self.v.der(),f.prod(f.abs(self.v),f.sqrt(f.sum(f.pow(self.v,2),-1))))
        def toStr(self): return 'arccsc(%s)' % self.v.toStr()
        def evalu(self,val): return asin(1/self.v.evalu(val))
    class arcsec:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(f.prod(-1,self.v.der()),f.prod(f.abs(self.v),f.sqrt(f.sum(f.pow(self.v,2),-1))))
        def toStr(self): return 'arcsec(%s)' % self.v.toStr()
        def evalu(self,val): return acos(1/self.v.evalu(val))
    class arccot:
        def __init__(self,fval): self.v = tofconst(fval)
        def der(self): return f.fract(f.prod(-1,self.v.der()),f.sum(1,f.pow(self.v,2)))
        def toStr(self): return 'arccot(%s)' % self.v.toStr()
        def evalu(self,val): return pi/2-atan(self.v.evalu(val))

splitters = '+ - ^ * / e exp ln √ log_ sin cos tan csc sec cot arcsin arccos arctan arccsc arcsec arccot abs'.split(' ')

def findLi(s,fds):
    ret = []
    for fd in fds:
        index = 0
        while index<len(s):
            index = s.find(fd,index)
            if index == -1:
                break
            ret.append([index,len(fd)])
            index+= len(fd)
    ret.sort()
    return ret
        

def fparseRe(li,var):
    newLi = []
    for n,el in enumerate(li):
        if type(el) == list:
            newLi.append(fparseRe(el,var))
        else:
            if el != '':
                spts = findLi(el,splitters+[var])
                ind = 0
                for i in spts:
                    if el[ind:i[0]]!= '':
                        newLi.append(el[ind:i[0]])
                    if el[i[0]:sum(i)] == '-' and (n!=0 or i[0]!=0):
                        newLi.append('+')
                        newLi.append('-')
                    else:
                        newLi.append(el[i[0]:sum(i)])
                    ind = sum(i)
                try:
                    if el[sum(spts[-1]):] != '':
                        newLi.append(el[sum(spts[-1]):])
                except: newLi.append(el)
    return newLi

fnames = '√ ln exp sin cos tan csc sec cot arcsin arccos arctan arccsc arcsec arccot abs'.split(' ')

def fparseNegFix(li):
    newLi = []
    noAdd = False
    for n,el in enumerate(li):
        if el == '-':
            try:
                newLi.append(str(-1*int(li[n+1])))
                noAdd = True
            except:
                newLi.append('-1');
        elif not noAdd:
            if type(el) == list:
                newLi.append(fparseNegFix(el))
            else:
                newLi.append(el)
        else: noAdd = False
    return newLi

def fgroup(li):
    newLi = []
    noAdd = 0
    for n,el in enumerate(li):
        for fname in fnames:
            if el == fname:
                if type(li[n+1]) == list:
                    newLi.append([el,fgroup(li[n+1])])
                else:
                    newLi.append([el,li[n+1]])
                noAdd = 2
        if noAdd == 0:
            if type(el) == list:
                newLi.append(fgroup(el))
            else:
                newLi.append(el)
        else: noAdd -= 1
    return newLi

def opgroup(li,st):
    if type(li)!= list:
        return li
    newLi = []
    noAdd = 0
    for n,el in enumerate(li):
        if el == st:
            del newLi[-1]
            if st == '^' and li[n-1] == 'e':
                newLi.append(['exp',opgroup(li[n+1],st)])
            else:
                newLi.append([el,opgroup(li[n-1],st),opgroup(li[n+1],st)])
            noAdd = 2
        if noAdd == 0:
            newLi.append(opgroup(el,st))
        else: noAdd -= 1
    return newLi
            
def sumprod(li):
    if type(li)!= list:
        return li
    newLi = []
    for el in li:
        newLi.append(sumprod(el))
    li = list(newLi)
    newLi = []
    index = 0
    for n,el in enumerate(li):
        if el == '+':
            if n-index == 1:
                newLi.append(li[index])
            else:
                newLi.append(['*']+li[index:n])
            index = n+1
    if len(li)-index == 1:
        newLi.append(li[index])
    else:
        newLi.append(['*']+li[index:])
    if len(newLi) == 1:
        newLi = newLi[0]
    else:
        newLi = ['+']+newLi
    if newLi[0] == '*':
        for spt in splitters:
            if newLi[1] == spt:
                del newLi[0]
                break
    return newLi

fnameDict = {'*': f.prod, '+': f.sum, 'exp': f.exp, 'ln': f.ln, '^': f.pow, '/': f.fract, '√':f.sqrt, 'sqrt':f.sqrt,
'abs': f.abs, 'sin':f.sin, 'cos':f.cos, 'tan':f.tan, 'csc':f.csc, 'sec':f.sec,
'cot':f.cot, 'arcsin':f.arcsin, 'arccos':f.arccos, 'arctan':f.arctan, 'arccsc':f.arccsc, 'arcsec':f.arcsec, 'arccot':f.arccot}

def fnreplace(li,var='x'):
    try:
        return f.const(int(li))
    except:
        try:
            return f.const(float(li))
        except:
            if li == var:
                return f.var(var)
    res = []
    for el in li[1:]:
        res.append(fnreplace(el))
    fname = fnameDict[li[0]]
    if li[0] == '*' or li[0] == '+':
        return fname(res)
    elif li[0] == '^':
        return fname(res[0],res[1])
    elif li[0] == '/':
        try: return f.const(frac(int(li[1]),int(li[2])))
        except: return fname(res[0],res[1])
    else:
        return fname(res[0])

def fparse(s, var='x'):
    sLi = eval("['"+s.replace('(',"',['").replace(')',"'],'")+"']")
    negFix = fparseNegFix(fparseRe(sLi,var))
    fnList = sumprod(opgroup(opgroup(fgroup(negFix),'^'),'/'))
    return fnreplace(fnList,var)
    

def derWT():
    print('Derivative Calculator in Python by Nathan Mihm with ~500 lines of code\nMake sure to include parenthesis to prevent the parser from getting confused\n')
    while True:
        fstring = input('f(x) = ')
        fparsed = fparse(fstring)
        print('\nResult (check to make sure your function was interpreted correctly):')
        
        print('f (x) = %s' % fparsed.toStr() )
        print('f\'(x) = %s' % fparsed.der().toStr() )
        if input('Press enter for another, or e to exit loop ') == 'e':
            break
            
derWT()
