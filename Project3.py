def create(m,n):
    l=[]
    for i in range(m):
        z=[]
        for j in range(n):
            x=input("Enter the a%d%d element of matrix:"%(i,j))
            if float(x)==int(x):
                y=int(x)
            else:
                y=float(x)
            z.append(y)
        l.append(z)
    return l
def disp(a):
    for i in a:
        print(i)
    print('\n')
def zero(m,n):
    O=[]
    for i in range(m):
        z=[]
        for j in range(n):
            z.append(0)
        O.append(z)
    return O
def isrow(a):
    return len(a)==1
def iscolumn(a):
    for i in range(len(a)):
        if len(a[i])!=1:
            return False
    else:
        return True
def issquare(a):
    for i in range(len(a)):
        if len(a[i])!=len(a):
            return False
    else:
        return True
def issymmetric(a):
    return transpose(a)==a
def isskewsymmetric(a):
    return transpose(a)==scale(a,-1)
def isdiag(a):
    f=issquare(a)
    for i in range(len(a)):
        for j in range(len(a[0])):
            if i!=j and a[i][j]!=0:
                return False
    return  f
def isscalar(a):
    f=isdiag(a)
    for i in range(len(a)):
        for j in range(len(a[0])):
            if i==j:
                if a[i][j]!=a[0][0]:
                    return False
    return f    
def isidentity(a):
    f=isscalar(a)
    for i in range(len(a)):
        for j in range(len(a[0])):
            if i==j:
                if a[i][j]!=1:
                    return False
    return f
def iszero(a):
    f=True
    for i in range(len(a)):
        for j in range(len(a[0])):
            if a[i][j]!=0:
                return False
    return f
def isnilpotent(a):
    return iszero(mpow(a,2))
def isinvolutory(a):
    return isidentity(mpow(a,2))
def isidempotent(a):
    return (mpow(a,2)==a)
def identity(m):
    I=[]
    for i in range(m):
        z=[] 
        for j in range(m):
            if i==j:
                z.append(1)
            else:
                z.append(0)
        I.append(z)
    return I
def msum(a,b):
    if len(a)!=len(b) or len(a[0])!=len(b[0]):
        raise Exception('Matrices have to be of the same order to be added')
    c=zero(len(a),len(b[0]))
    for i in range(len(a)):
        for j in range(len(b[i])):
            c[i][j]=a[i][j]+b[i][j]
    return c
def mdiff(a,b):
    if len(a)!=len(b) or len(a[0])!=len(b[0]):
        raise Exception('Matrices have to be of the same order to be subtracted')
    c=zero(len(a),len(b[0]))
    for i in range(len(a)):
        for j in range(len(b[i])):
            c[i][j]=a[i][j]-b[i][j]
    return c
def mprod(a,b):
    if len(a[0])!=len(b):
        raise Exception('Matrices cannot be multiplied')
    c=zero(len(a),len(b[0]))
    for i in range(len(a)):
        for j in range(len(b[i])):
            for k in range(len(b)):
                c[i][j]+=a[i][k]*b[k][j]
    return c
def transpose(a):
    b=[]
    for i in range(len(a[0])):
        z=[] 
        for j in range(len(a)):
            z.append(None)
        b.append(z)
    for i in range(len(a)):
        for j in range(len(a[i])):
            b[j][i]=a[i][j]
    return b
def scale(a,k):
    l=[]
    for r in a:
        x=list(r)
        l.append(x)    
    for i in range(len(l)):
        for j in range(len(l[i])):
            l[i][j]*=k
    return l  
def mpow(a,x):
    b=a
    for i in range(x-1):
        b=mprod(a,b)
    return b
def fact(n):
    f=1
    if n<0:
        for i in range (n,0):
            f=i*f
    else:
        for i in range (1,n+1):
            f=i*f
    return f            
def mexp(a):
    c=identity(len(a))
    for x in range(1,100):
        c=msum(c,scale(mpow(a,x),1/fact(x)))
    return c
def delete(a,i,j):
    l=[]
    for r in a:
        x=list(r)
        l.append(x)
    for k in range(1,len(a)+1):
        del l[k-1][j-1]
    del l[i-1]
    return l 
def det(a):
    if not issquare(a):
        raise Exception('Determinant of non-square matrices is not defined')
    elif len(a)==1:
        return a[0][0]
    else:
        d=0
        for j in range(1,len(a)+1):
            e=a[0][j-1]
            b=delete(a,1,j)
            c=det(b)
            d+=((-1)**(j-1))*c*e
        return d
def minor(a,i,j):
    return det(delete(a,i,j))
def cofactor(a,i,j):
    return (-1)**(i+j)*(minor(a,i,j))
def adj(a):
    b=[]
    for i in range(1,len(a)+1):
        r=[]
        for j in range(1,len(a[i-1])+1):
            x=cofactor(a,i,j)
            r.append(x)
        b.append(r)
    return b
def inverse(a):
    return scale(transpose(adj(a)),(1/det(a)))
def tr(a):
    s=0
    for i in range(len(a)):
        for j in range (len(a[i])):
            if i==j:
                s+=a[i][j]
    return s
def balance():
    #Program balances chemical reactions with elements 1 less than the number of compounds
    a=input('Enter the Equation to be balanced: ')
    #Identification of elements in the Equation
    element=[]
    for i in range(len(a)): 
        if a[i].isupper():
            z=a[i]
            try:
                if a[i+1].islower():
                    z+=a[i+1]
            except IndexError:
                if z not in element:
                    element.append(z)
        if z not in element:
            element.append(z)        
    #Identification of Reactants and Products and conversion of chemical reaction to Data in matrix form
    c=list(a.partition('->'))
    react=c[0].split('+')
    Eq={}
    for x in element:
        l=[]
        for i in react:
            t=1
            if '(' in i :
                lig=i[i.index('(')+1:i.index(')')]
                if x in lig:
                    t=int(i[i.index(')')+1])
            if x not in i:
                n=0
            elif i.index(x)+len(x)<len(i): 
                n=1
                if i[i.index(x)+len(x)].isdigit():
                    m=0
                    subs='0'
                    while i.index(x)+len(x)+m<len(i):
                        if i[i.index(x)+len(x)+m].isdigit():
                            subs+=i[i.index(x)+len(x)+m]
                            m+=1
                        else:
                            break
                    n=int(subs)
            else:
                n=1
            l.append(n*t)
        Eq[x]=l
    prod=c[2].split('+')
    for x in element:
        l=[]
        for i in prod:
            t=1
            if '(' in i :
                lig=i[i.index('(')+1:i.index(')')]
                if x in lig:
                    t=int(i[i.index(')')+1])
            if x not in i:
                n=0
            elif i.index(x)+len(x)<len(i):
                n=1
                if i[i.index(x)+len(x)].isdigit():
                    m=0
                    subs='0'
                    while i.index(x)+len(x)+m<len(i):
                        if i[i.index(x)+len(x)+m].isdigit():
                            subs+=i[i.index(x)+len(x)+m]
                            m+=1
                        else:
                            break
                    n=int(subs)
            else:
                n=1
            n=-n
            l.append(n*t  )
        Eq[x]+=l
    b=input("Enter the Compound whose stoichometric co-efficient is one: ")
    #Solving the system using matrix method and finding the stoichometric co-efficients for balancing the equation
    comp=react+prod
    const=comp.index(b)
    out=[]
    for i in Eq:
        out.append([-Eq[i].pop(const)])
    chem=[]
    for i in Eq:
        chem.append(Eq[i])
    co_eff=(mprod(inverse(chem),out))
    co_eff.insert(const,[1.0])
    for i in range(len(comp)):
        comp[i]=str(co_eff[i][0])+' '+comp[i]
    #String manipulation for output of balanced Equation
    balanced=''
    for i in react:
        balanced+=comp[0]+' +'
        del comp[0]
    balanced=balanced[0:-1]
    balanced+=' -> '
    for i in prod:
        balanced+=comp[0]+' +'
        del comp[0]
    balanced=balanced[0:-1]
    print(balanced)
#Binary file manipulation for storage and accessing of matrices
import pickle
try:
    f=open('C:\\Users\\tanus\\OneDrive\\Desktop\\IIT\\comp\\matrix.dat','rb')
except FileNotFoundError:
    f=open('C:\\Users\\tanus\\OneDrive\\Desktop\\IIT\\comp\\matrix.dat','r')
try:
    d=pickle.load(f)
except EOFError:
    d={}
def save(m):
    hl=input("Would you like to save the matrix?(y/n):")
    if hl=='y':
        name=input("Enter the name of the matrix: ")
        d[name]=m
        f=open('C:\\Users\\tanus\\OneDrive\\Desktop\\IIT\\comp\\matrix.dat','wb')
        pickle.dump(d,f)
        f.close()
def access(m):
    try:
        matrix=eval(m)
        return matrix
    except NameError:
        try:
            matrix=d[m]
            return matrix
        except KeyError:
            return("Matrix not found")
while True:
    print(""" \n 1.Create a matrix \n 2.Access A matrix \n 3.Identify the type of matrix \n 4.Matrix calculator \n 5.Balance chemical equations \n 6.Convert a matrix into a table \n 7.End Program """)
    try:
        Job = int(input("what would you like to do?"))
    except ValueError:
        print("Invalid input")  
    if Job==1:
        order=input("Enter the order of the matrix(mxn)")
        matrix=create(int(order[0]),int(order[2]))
        disp(matrix)
        save(matrix)
    elif Job==2:
        name=input("Enter the name of the matrix to be accessed:")
        if type(access(name))==list:
            disp(access(name))
        else:
            print(access(name))
    elif Job==3:
        m=input("Enter the matrix")
        mtype=''
        m=access(m)
        if isnilpotent(m):
            mtype+='Nilpotent, '
        if isinvolutory(m):
            mtype+='involutory, '
        if iszero(m):
            mtype+='Zero, '
        elif issquare(m):
            if isdiag(m):
                if isscalar(m):
                    if isidentity(m):
                        mtype+='Identity, '
                    else:
                        mtype+='Scalar, '
                else:
                    mtype+='Diagonal, '
            elif issymmetric(m):
                mtype+='Symmetric, '
            elif isskewsymmetric(m):
                mtype+='Skew-symmetric, '
            else:
                mtype+='Square, '
        elif isrow(m):
            mtype+='Row, '
        elif iscolumn(m):
            mtype+='Column, '
        mtype=mtype[0:-2]
        print('The given matrix is a/an '+mtype+' matrix')
    elif Job==4:
        while True:
            print(f"""
                    |{'Welcome to Matrix Calculator'.center(80,'~')}|
                    |{'What would you like to do?:'.center(80,'~')}|
                    |1. Addition of matrices{' '*(80-23)}|
                    |2. Subtraction of matrices{' '*(80-26)}|
                    |3. Multiplication of matrices{' '*(80-29)}|
                    |4. Scalar multiplication{' '*(80-24)}|
                    |5. Transpose of matrix{' '*(80-22)}|
                    |6. Adjoint of a Matrix{' '*(80-22)}|
                    |7. Inverse of a Matrix{' '*(80-22)}|
                    |8. Determinant of a matrix{' '*(80-26)}|
                    |9. nth Power of a matrix{' '*(80-24)}|
                    |10. Exponential of a matrix{' '*(80-27)}|
                    |11. Trace of the matrix{' '*(80-23)}|
                    |12. Exit Matrix Calculator{' '*(80-26)}|
                    |{'-'*80}|
                    """)
            n=int(input("Enter your Choice: "))
            if n == 1:
                m1=input("Enter the first matrix")
                m1=access(m1)
                m2=input("Enter the second matrix")
                m2=access(m2)
                disp(msum(m1,m2))
                save(msum(m1,m2))
            elif n==2:
                m1=input("Enter the first matrix")
                m1=access(m1)
                m2=input("Enter the second matrix")
                m2=access(m2)
                disp(mdiff(m1,m2))
                save(mdiff(m1,m2))
            elif n==3:
                m1=input("Enter the first matrix")
                m1=access(m1)
                m2=input("Enter the second matrix")
                m2=access(m2)
                disp(mprod(m1,m2))
                save(mprod(m1,m2))
            elif n==4:
                m1=input("Enter the matrix")
                m1=access(m1)
                k=int(input("Enter the scalar to be multiplied"))
                disp(scale(m1,k))
                save(scale(m1,k))
            elif n==5:
                m1=input("Enter the matrix")
                m1=access(m1)
                disp(transpose(m1))
                save(transpose(m1))
            elif n==6:
                m1=input("Enter the matrix")
                m1=access(m1)
                disp(adj(m1))
                save(adj(m1))
            elif n==7:
                m1=input("Enter the matrix")
                m1=access(m1)
                disp(inverse(m1))
                save(inverse(m1))
            elif n==8:
                m1=input("Enter the matrix")
                m1=access(m1)
                print(det(m1))
            elif n==9:
                m1=input("Enter the matrix")
                m1=access(m1)
                p=int(input("Enter the power to which the matrix is to be raised"))
                print(mpow(m1,p))
                save(mpow(m1,p))
            elif n==10:
                m1=input("Enter the matrix")
                m1=access(m1)
                print(mexp(m1))
                save(mexp(m1))
            elif n==11:
                m1=input("Enter the matrix")
                m1=access(m1)
                print(tr(m1))
            elif n==12:
                break
            else:
                print("Invalid response")
    elif Job==5:
        balance()
    elif Job==6:
        import mysql.connector
        conn = mysql.connector.connect(host='localhost', password='mathanos', user='root')
        if conn.is_connected():
            print("Connection established...")
        c=conn.cursor()
        c.execute('CREATE DATABASE IF NOT EXISTS Matrices')
        c.execute('Use Matrices')
        matrix=eval(input("Enter the matrix to be converted into a table"))
        table_name=input("Enter the name of the table")
        com='Create Table '+table_name +'('
        for i in matrix[0]:
            n=input("Enter name of the column")
            com+=n+' int(9), '
        com=com[:-2]
        com+=')'
        c.execute(com)
        for i in matrix:
            com='Insert Into '+table_name+' Values('
            row=''
            for j in i:
                row+=str(j)+', '
            com+=row[:-2]+')'
            c.execute(com)
            conn.commit()
        com='Select * from'+table_name
        c.execute(com)
    elif Job==7:
        break
    else:
        print("Invalid response")
