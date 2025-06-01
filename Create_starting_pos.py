#Creating starting Position for a 1 CD3 Complex with all chains 
import pandas as pd
import numpy as np
#calc starting position for each chain 
#position_complex.txt gives information about the cyroem studies starting position
zeta_a= [(103.46000, 121.81600, 158.16900),
(102.79800, 120.80800, 158.98800),
(103.38400, 120.76900, 160.39300),
(102.64200, 120.67800, 161.37800),
(102.13900, 118.34500, 159.05000),
(100.67400, 118.35800, 158.65000),
(99.84000,  118.98300, 159.67000), 
(98.59700,  119.39800, 159.47100),
(98.00800,  119.27400, 158.29300),
(97.92900,  119.95100, 160.47900)]
Zeta_b=[
(115.46700, 124.45800, 161.71100), 
(115.40600, 124.97000, 163.07700),
(114.10100, 124.59800, 163.77700),
(114.10900, 124.24400, 164.9630),
(115.59900, 126.47900, 163.10200),
(116.29500, 126.88500, 164.35500),
(117.64700, 126.24500, 164.42800),
(118.18000, 126.38700, 165.77100),
(118.79800, 127.47400, 166.19200),
(118.90500, 128.53300, 165.41700),
(119.29800, 127.50300, 167.42300)]
Delta_a =[
(125.89900, 115.81200, 165.12600),
(126.01600, 115.88700, 166.57200),
(125.18800, 116.99600, 167.19000),
(124.67700, 116.84900, 168.30300)]
Epsilon_a = [
(128.62300, 142.98100, 164.31300),
(127.66700, 142.13300, 165.01700),
(128.21000, 141.67600, 166.36300),
(127.47600, 141.65300, 167.35700),
(127.30100, 140.92700, 164.15800),
(126.40100, 140.08300, 164.85000)]
Epsilon_b=[
(126.02200, 108.61600, 161.89500),
(126.56000, 109.56300, 162.86200),
(125.56800, 109.95800, 163.94700),
(125.98700, 110.50300, 164.97400),
(127.04600, 110.82700, 162.14900),
(128.12500, 110.53800, 161.27900)]
Gamma_a= [
(127.73600, 135.18600, 164.41900),
(128.66200, 136.05000, 165.13600),
(128.47000, 135.89400, 166.63800),
(127.34500, 135.97500, 167.14200),
(128.46400, 137.50700, 164.72500)]

def get_avr_pos(data):
    #get average position
    x_sum=0
    y_sum=0
    z_sum=0
    for i in range(len(data)):
        x_sum += data[i][0]
        y_sum += data[i][1]
        z_sum += data[i][2]
    return x_sum/len(data), y_sum/len(data), z_sum/len(data)

#get average position of an amino acid from the coordinates given by cryo-em study
Zeta_a_avr = get_avr_pos(zeta_a)
Zeta_b_avr = get_avr_pos(Zeta_b)
Delta_a_avr = get_avr_pos(Delta_a)
Epsilon_a_avr = get_avr_pos(Epsilon_a)
Epsilon_b_avr = get_avr_pos(Epsilon_b)
Gamma_a_avr = get_avr_pos(Gamma_a)

Zeta_b_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Zeta_b_avr))
Delta_a_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Delta_a_avr))
Epsilon_a_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Epsilon_a_avr))
Epsilon_b_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Epsilon_b_avr))
Gamma_a_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Gamma_a_avr))
Zeta_a_start = tuple(map(lambda i, j: i - j, Zeta_a_avr, Zeta_a_avr))

#define length of the chains
N_zeta= 164-52
N_epilon= 208-153
N_delta= 171-127
N_gamma= 182-137
N_gesamt = 2*N_zeta+2*N_epilon+N_gamma+N_delta

#define simulation box size
xlo=-115*6.4
xhi=115*6.4
ylo=-115*6.4
yhi=115*6.4
zlo=-115*6.4
zhi=50
#filepath and name of the in putfile
filepath = "example"
with open(filepath,'w') as fdata:
    
    #headers
    fdata.write('polymere chain with 150 points\n\n')
    fdata.write('{} atoms\n'.format(int((N_gesamt))))
    fdata.write('{} bonds\n'.format(int((N_gesamt-6))))
    fdata.write('{} atom types\n'.format(3))
    fdata.write('{} bond types\n'.format(1))
    fdata.write('{} {} xlo xhi\n'.format(xlo,xhi))
    fdata.write('{} {} ylo yhi\n'.format(ylo,yhi))
    fdata.write('{} {} zlo zhi\n\n'.format(zlo,zhi))
    fdata.write('Masses\n\n')
    fdata.write('{} {}\n'.format(1,0.00015))
    fdata.write('{} {}\n'.format(2,0.00015))
    fdata.write('{} {}\n\n'.format(3,0.00015))
    fdata.write('Atoms\n\n')
    i=1
    g=1

#starting coordinates for each unit
#first CD3 complex
    #first zeta chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Zeta_a_start[0],Zeta_a_start[1],Zeta_a_start[2],0,0,0))
    g+=1
    while i < N_zeta:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,1,Zeta_a_start[0],Zeta_a_start[1],float(Zeta_a_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    #second zeta chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Zeta_b_start[0],Zeta_b_start[1],Zeta_b_start[2],0,0,0))
    g+=1
    while i < N_zeta:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,1,Zeta_b_start[0],Zeta_b_start[1],float(Zeta_b_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    #First Delta chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Delta_a_start[0],Delta_a_start[1],Delta_a_start[2],0,0,0))
    g+=1
    while i < N_delta:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,2,Delta_a_start[0],Delta_a_start[1],float(Delta_a_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    #First Epsilon chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Epsilon_a_start[0],Epsilon_a_start[1],Epsilon_a_start[2],0,0,0))
    g+=1
    while i < N_epilon:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,1,Epsilon_a_start[0],Epsilon_a_start[1],float(Epsilon_a_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    #Second Epsilon chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Epsilon_b_start[0],Epsilon_b_start[1],Epsilon_b_start[2],0,0,0))
    g+=1
    while i < N_epilon:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,1,Epsilon_b_start[0],Epsilon_b_start[1],float(Epsilon_b_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    #First Gamma chain
    fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,3,Gamma_a_start[0],Gamma_a_start[1],Gamma_a_start[2],0,0,0))
    g+=1
    while i < N_gamma:
        fdata.write('{} {} {} {} {} {} {} {} {}\n'.format(g,1,2,Gamma_a_start[0],Gamma_a_start[1],float(Gamma_a_start[2]-6.2*(i)),0,0,0))
        i+=1
        g+=1
        u = 1
    i=1
    K=0
    
    #stargin velocities
    fdata.write('\n')
    fdata.write('Velocities\n\n')
    i=1
    while i < g:
        fdata.write('{} {} {} {}\n'.format(i,0,0,0))
        i+=1
    fdata.write('\n')
    
    #define bonds 
    fdata.write('Bonds\n\n')
    Bondnumber = 1
    i=1
    while i < N_zeta:
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
    i+=1
    while i < 2*N_zeta:
        
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
    i+=1
    while i < N_delta+2*N_zeta:
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
    i+=1
    while i < N_epilon+N_delta+2*N_zeta:
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
    i+=1
    while i < 2*N_epilon+N_delta+2*N_zeta:
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
    i+=1
    while i < N_gamma+2*N_epilon+N_delta+2*N_zeta:
        fdata.write('{} {} {} {}\n'.format(Bondnumber,1,i,i+1))
        i+=1
        Bondnumber += 1
#check if bondnumer is correct
print(Bondnumber,N_gesamt-5)        