import random
import json
class CFG:
    scale=100
    num=10
    n_fq=5

def fun1(id):
    Data={'T':[],'D':[]}
    for i in range(CFG.num):
        T=random.random()*CFG.scale
        Data['T'].append(T)
        Data['D'].append(random.random()*T/CFG.num)
    with open('data/task1/{0}.json'.format(id),'w',encoding='utf-8') as f:
        json.dump(Data,f)

    with open('data/task1/{0}.txt'.format(id),'w',encoding='utf-8') as f:
        f.write('{0} {1:.4f}\n'.format(CFG.num,random.random()*CFG.scale/(CFG.num*5)))
        for i in range(CFG.num):
            f.write("{0:.4f} {1:.4f}\n".format(Data['D'][i],Data['T'][i]))

def fun2(id):
    Data={'phi':[],'D':[],'T':[],'f':[]}
    for i in range(CFG.num):
        T=random.random()*CFG.scale
        Data['T'].append(T)
        Data['D'].append(random.random()*T/CFG.num)
        Data['f'].append(random.randint(0,CFG.n_fq-1))
        Data['phi'].append(random.random()*T)
    
    with open('data/task2/{0}.txt'.format(id),'w',encoding='utf-8') as f:
        f.write('{0} {1} {2:.4f}\n'.format(CFG.num,CFG.scale*5,random.random()*CFG.scale/(CFG.num*5)))
        for i in range(CFG.num):
            f.write("{0:.4f} {1:.4f} {2:.4f} {3}\n".format(Data['phi'][i],Data['D'][i],Data['T'][i],Data['f'][i]))

def fun3_1(id):
    Data={'T':[],'D':[]}
    for i in range(CFG.num):
        T=random.random()*CFG.scale
        Data['T'].append(T)
        Data['D'].append(random.random()*T/CFG.num)
    pT=random.random()*CFG.scale
    pD=random.random()*pT/CFG.num
    with open('data/task3/1_{0}.txt'.format(id),'w',encoding='utf-8') as f:
        f.write('{0} {1:.4f}\n'.format(CFG.num,random.random()*CFG.scale/(CFG.num*10)))
        for i in range(CFG.num):
            f.write("{0:.4f} {1:.4f}\n".format(Data['D'][i],Data['T'][i]))
        f.write("\n{0:.4f} {1:.4f}\n".format(pD,pT))
        
def fun3_2(id):
    Data={'phi':[],'D':[],'T':[],'f':[]}
    for i in range(CFG.num):
        T=random.random()*CFG.scale
        Data['T'].append(T)
        Data['D'].append(random.random()*T/CFG.num)
        Data['f'].append(random.randint(0,CFG.n_fq-1))
        Data['phi'].append(random.random()*T)
    
    pT=random.random()*CFG.scale
    pD=random.random()*pT/CFG.num
    pf=random.randint(0,CFG.n_fq-1)
    pphi=random.random()*pT

    with open('data/task3/2_{0}.txt'.format(id),'w',encoding='utf-8') as f:
        f.write('{0} {1} {2:.4f}\n'.format(CFG.num,CFG.scale*5,random.random()*CFG.scale/(CFG.num*5)))
        for i in range(CFG.num):
            f.write("{0:.4f} {1:.4f} {2:.4f} {3}\n".format(Data['phi'][i],Data['D'][i],Data['T'][i],Data['f'][i]))
        f.write("\n{0:.4f} {1:.4f} {2:.4f} {3}\n".format(pphi,pD,pT,pf))

if __name__=='__main__':
    for i in range(100):
        fun3_1(i)
