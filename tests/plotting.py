import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



for filename in ['vector_addition','daxpy','finite_difference','multiquadric','multiquadric_scaling','linear_spring']:
#for filename in ['linear_spring']:
    if filename in ['linear_spring']:
        files = glob.glob('tests/'+filename+'*.csv')
        datas = []
        radius_div_h = []
        for f in sorted(files):
            radius_div_h.append(float(re.search('tests\/'+filename+'(.*)\.csv',f).group(1)))
            print 'reading ',f,' with r = ',radius_div_h[-1]
            datas.append(np.loadtxt(f))
    else:
        data = np.loadtxt('tests/'+filename+'.csv')
    plt.figure(figsize= (6,5))

    handles = []
    if filename in ['vector_addition']:
        plt.title(r'$a_i = b_i + c_i$')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria (Level 1)')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria (Level 2)')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')
    elif filename in ['daxpy']:
        plt.title(r'$a_i = a_i + 0.1\, b_i $')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria (Level 1)')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria (Level 2)')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')
    elif filename in ['multiquadric']:
        plt.title(r'$s_i = s_i+ \sum_j^N  s_j \sqrt{|\mathbf{dx}_{ij}|+c_j^2}$')
        filename_matlab = 'matlab_matrix_multiply_size'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria')
        #plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria-eigen')
        #plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
        plt.semilogx(data[:,0],data[:,4]/1e6,label='std::vector')
        #plt.semilogx(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')
    elif filename in ['multiquadric_scaling']:
        plt.title(r'$s_i = s_i+ \sum_j^N  s_j \sqrt{|\mathbf{dx}_{ij}|+c_j^2}$')
        filename_matlab = 'matlab_matrix_multiply_scaling'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.plot(data[:,0],data[:,1]/1e6,label='Aboria')
        #plt.plot(data[:,0],data[:,2]/1e6,label='Aboria-Eigen')
        #plt.plot(data[:,0],data[:,3]/1e6,label='Eigen')
        plt.plot(data[:,0],data[:,4]/1e6,label='std::vector')
        plt.plot(data_matlab[:,0],data_matlab[:,1]/1e6,label='Matlab')
    elif filename in ['linear_spring']:
        plt.title(r'$dr_i = \sum_{|\mathbf{dx}_{ij}|<2r} \frac{(2r-|\mathbf{dx}_{ij}|}{|\mathbf{dx}_{ij}|}\mathbf{dx}_{ij}$')
        dashed_line = mpl.lines.Line2D([], [], color='black',linestyle='--',label='std::vector + Gromacs search')
        solid_line = mpl.lines.Line2D([], [], color='black',linestyle='-',label='Aboria')
        handles.append(dashed_line)
        handles.append(solid_line)
        for data,radius in zip(datas,radius_div_h):
            handle, = plt.loglog(data[:,0],data[:,1]/1e6,label='radius = %2.1f'%radius)
            handles.append(handle)
            plt.loglog(data[:,0],data[:,2]/1e6,color=handle.get_color(),linestyle='--',label='gromacs %f'%radius)
        plt.ylim(10,1e4)
    else:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria-Eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')


    if filename in ['multiquadric_scaling']:
        leg = plt.legend(fancybox=True,loc='upper right')
    elif filename in ['linear_spring']:
        leg = plt.legend(handles=handles,fancybox=True,loc='upper left')
    elif filename in ['daxpy']:
        leg = plt.legend(fancybox=True,loc='upper right')
    else:
        leg = plt.legend(fancybox=True,loc='upper left')

    leg.get_frame().set_alpha(0.3)
    plt.xlabel(r'$N$')
    if filename in ['vector_addition', 'daxpy']:
        plt.ylabel(r'$N / T_e$ ($\times 10^6$)')
    elif filename in ['linear_spring']:
        plt.ylabel(r'$N^3 / T_e$ ($\times 10^6$)')
    else:
        plt.ylabel(r'$N^2 / T_e$ ($\times 10^6$)')
    plt.tight_layout()
    plt.savefig(filename+'.pdf')


