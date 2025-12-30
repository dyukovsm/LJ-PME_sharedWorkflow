import signac
import numpy as np
import os
current_directory = os.getcwd()
current_directory_name = os.path.basename(current_directory)
project = signac.init_project() #current_directory)

r_cut = [6.0, 4.0, 3.0, 1.4] # [6.0, 4.0, 3.0, 2.0, 1.4, 1.0] 
cut_type = ['Cut-off','PME'] #, 138.1, 72.53] # 343, 310.3, 241.7, 138.1 and 72.53
temperature = [135.0, 140.0, 145.0, 165.0]
replicas = [0]

# filter the list of dictionaries
total_statepoints = list()
legend = open('legend.txt','w')
legend.write('job \t sp \n')
print('job \t sp ')

for i in range(len(r_cut)):
    if r_cut[i] == 1.4 : # and r_cut[i] > 0.99 :
        for j in range(len(cut_type)):
            for k in range(len(replicas)):
                for l in range(len(temperature)):
                    statepoint = {
                        "r_cut": r_cut[i],
                        "cut_type": cut_type[j],
                        "replicas": replicas[k],
                        "temperature": temperature[l]
                    }
                    total_statepoints.append(statepoint)
    else:
        #for j in range(len(cut_type)):
            for k in range(len(replicas)):
                for l in range(len(temperature)):
                    statepoint = {
                        "r_cut": r_cut[i],
                        "cut_type": cut_type[0],
                        "replicas": replicas[k],
                        "temperature": temperature[l]
                    }
                    total_statepoints.append(statepoint)
        
for sp in total_statepoints:
    job=project.open_job(
        statepoint=sp,
    ).init()
    legend.write(f'{job} \t\t {sp}\n')
    print(f'{job} \t\t {sp}')
    
legend.close()

