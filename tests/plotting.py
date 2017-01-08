import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import glob
import re

#for filename in ['vector_addition','daxpy','finite_difference','multiquadric','multiquadric_scaling','linear_spring']:
for filename in ['linear_spring']:
    if filename in ['linear_spring']:
        files = glob.glob('tests/'+filename+'*.csv')
        datas = []
        radius_div_h = []
        for f in files:
            radius_div_h.append(float(re.search('tests\/'+filename+'(.*)\.csv',f).group(1)))
            print 'reading ',f,' with r = ',radius_div_h[-1]
            datas.append(np.loadtxt(f))
    else:
        data = np.loadtxt('tests/'+filename+'.csv')
    plt.figure(figsize= (6,5))
    if filename in ['vector_addition','daxpy']:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria_level1')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_level2')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
    elif filename in ['multiquadric']:
        filename_matlab = 'matlab_matrix_multiply_size'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
        plt.semilogx(data[:,0],data[:,4]/1e6,label='std::vector')
        plt.semilogx(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')
    elif filename in ['multiquadric_scaling']:
        filename_matlab = 'matlab_matrix_multiply_scaling'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.plot(data[:,0],data[:,1]/1e6,label='aboria')
        plt.plot(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.plot(data[:,0],data[:,3]/1e6,label='eigen')
        plt.plot(data[:,0],data[:,4]/1e6,label='std::vector')
        plt.plot(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')
    elif filename in ['linear_spring']:
        for data,radius in zip(datas,radius_div_h):
            plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria %f'%radius)
            plt.semilogx(data[:,0],data[:,2]/1e6,linestyle='--',label='gromacs %f'%radius)
        print 'dummy'
    else:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')


    if filename in ['multiquadric_scaling']:
        leg = plt.legend(fancybox=True,loc='upper right')
    else:
        leg = plt.legend(fancybox=True,loc='upper left')
    leg.get_frame().set_alpha(0.3)
    plt.xlabel('N')
    plt.ylabel('MFLOPS')
    plt.tight_layout()
    plt.savefig(filename+'.pdf')


