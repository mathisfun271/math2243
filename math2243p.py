#math2243
#programmed by Nathan Mihm

#although numpy has many things preprogrammed, I am only using:
#array: creates a numpy array (matrix). Simplifies row opperations significantly from generic Python.
#random: Used to test these functions against random matricies. 
from numpy import array,random,zeros #only basic things imported from numpy
from math import gcd
from fractions import Fraction as frac

def mgcd(*args): #gcd for any amount of values
    #can interpret a several arguments or a single list.
    if len(args)==1 and type(args[0])==list:
            args=args[0]
    if len(args)==2:
        return gcd(args[0],args[1])
    else:
        ar = []#cuts number of terms in half
        for i in range(0, len(args), 2):
            try:ar.append(gcd(args[i],args[i+1]))
            except:ar.append(args[i])
        return mgcd(ar)#then nests back, until there are only 2 terms.


def swap(matrix,rowswap1,rowswap2):#swap rows
    matrix[[rowswap1,rowswap2]]=matrix[[rowswap2,rowswap1]]

def ref(matrix):#row echelon form
    a=array(matrix)#sets 'a' to a new array to not disrupt the old one.
    yc=0#y check in matrix.
    for x in range(0,len(a[0])):
        nozero = True
        if a[yc,x] == 0:#makes sure there in not a zero in the position we check,
            for y in range(yc+1,len(a)+1):
                if y==len(a):nozero = False#if the entire column after that is zero, we increment the x
                elif a[y,x]!=0:
                    swap(a,yc,y)#if it is not, we swap the non-zero row with the checking row.
                    break
        if nozero: #only runs if there is no zero. If there is, the loop increments x and tries again.
            for y in range(yc+1,len(a)):
                if a[y,x] != 0:#if this is zero, we are already done
                    a[y]=(a[yc]*a[y,x]-a[y]*a[yc,x])//gcd(a[yc,x],a[y,x]) #makes this position 0.
                    #works by multiplying each row to the lcm, then subtracing.
            yc+=1
        if yc>=len(a):break #end loop once all rows are reached
    for y in a:#checks through each row
        for x in y:
            if x!=0:#at the leading term:
                y //= mgcd(list(y))#divide by gcd of row
                if x<0:y*=-1#make the leading term positive by negating row.
                break
    return a
        
def rref(matrix,dtype='frac',prec=5): #reduced row echelon form
    a = ref(matrix).astype('object')#allow adding fractions/decimals into matrix
    #These are only needed for the final step, making the leading coefficient one.
    if dtype!='frac' and dtype!='dec': raise Exception("Use either 'frac' or 'dec' for dtype")
    for yc in range(0,len(a))[::-1]:#iterates rows bottom to top.
        allzeros = False
        for x in range(0,len(a[0])+1):#finds leading coefficients
            if x==len(a[0]):allzeros = True;break
            if a[yc,x] != 0:xc=x;break
        if (not allzeros):
            for y in range(0,yc):#iterates through each value in the column coresponding to the leading coefficients
                if a[y,xc] != 0:#very similar to method used earlier in ref() to clear values.
                    a[y]=(a[y]*a[yc,xc]-a[yc]*a[y,xc])//gcd(a[yc,xc],a[y,xc])
            for x in range(xc+1,len(a[0])):#iterates through rest of row after leading coefficient
                if a[yc,x] % a[yc,xc] == 0:
                    a[yc,x] //= a[yc,xc]# if it is divisible, it returns the integer division.
                else:#if it is not it does either of these
                    if dtype=='frac':#returns a fraction
                        a[yc,x] = frac(a[yc,x],a[yc,xc])
                    elif dtype=='dec':#returns a decimal
                        a[yc,x] = round(a[yc,x]/a[yc,xc],prec)
            a[yc,xc] = 1#sets leading coefficient to one, now that we have divided the rest of the row.
    return a



#variables, separated by commas in a string, (spaces/no spaces, or as tuple)
#(vars muust be single letter, with an optional number following)
#list of equations (as strings)
def sysSolve(variables, *equations):
    vs = variables.replace('(','').replace(')','').replace(' ','').split(',')
    vdict  = dict((v,n) for n,v in enumerate(vs))
    vdict1 = dict((n,v) for n,v in enumerate(vs))
    mat = zeros((len(equations),len(vs)+1),dtype=int)
    for n,eq in enumerate(equations):
        eq = eq.replace(' ','')
        const = int(eq.split('=')[1])
        mons = eq.split('=')[0].replace('-','+-').split('+')
        try: mons.remove('')
        except: pass
        for mon in mons:
            for i in range(0,len(mon)):
                for v in vs:
                    if mon[-i:]==v:
                        try: mat[n,vdict[v]]=int(mon[:i])
                        except:
                            if mon[:i] == '': mat[n,vdict[v]]=1
                            elif mon[:i] == '-': mat[n,vdict[v]]=-1
        mat[n,-1]=const

    sts = []
    mat=rref(mat)
    for row in list(mat):
        for n,val in enumerate(row):
            if val == 1:
                sts.append('%s=%s' % (vdict1[n],row[-1]))
                break
                
    return ', '.join(sts)

    
