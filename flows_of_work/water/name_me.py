

#naming_file = open('aggregate_general_Data.txt','w')
#
#job_str = ['']*5
#
#job_str[0] = "jobid"
#job_str[1] = "rcut"
#job_str[2] = "cuttype"
#job_str[3] = "replicas"
#job_str[4] = "temp. [K]"
#
#properties_of_interest = ['']*8
#
#properties_of_interest[0] = 'Potential'
#properties_of_interest[1] = 'LJ (SR)'
#properties_of_interest[2] = 'Coulomb (SR)'
#properties_of_interest[3] = 'Coul. recip.'
#properties_of_interest[4] = 'Total Energy'
#properties_of_interest[5] = 'Vir-ZZ'
#properties_of_interest[6] = 'Pres-ZZ'
#properties_of_interest[7] = '#Surf*SurfTen'
#
#naming_file.write(f"{job_str[0]:<42} {job_str[1]:<8} {job_str[2]:<8} {job_str[3]:<8}   {job_str[4]:<8}"
#                                   f" {properties_of_interest[0]:<42} " 
#                                   f" {properties_of_interest[1]:<42} " 
#                                   f" {properties_of_interest[2]:<42} " 
#                                   f" {properties_of_interest[3]:<42} "
#                                   f" {properties_of_interest[4]:<42} "
#                                   f" {properties_of_interest[5]:<42} "
#                                   f" {properties_of_interest[6]:<42} "
#                                   f" {properties_of_interest[7]:<42} "
#                                   "\n\n")
#

naming_file = open('aggregate_dens_Data.txt','w')
job_str = ['']*5
job_str[0] = "jobid"
job_str[1] = "rcut"
job_str[2] = "cuttype"
job_str[3] = "replicas"
job_str[4] = "temp. [K]"
job_str[5] = 'liq dens'
job_str[6] = 'gas dens'

naming_file.write(f"{job_str[0]:<42} {job_str[1]:<8} {job_str[2]:<8} {job_str[3]:<8}   {job_str[4]:<8}"
                  f" \t\t {job_str[5]:<9}"
                  f" \t\t {job_str[6]:<9}"
                  "\n\n")

